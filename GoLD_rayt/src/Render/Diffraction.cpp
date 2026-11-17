// Diffraction at the aperture (--diffraction); see Render/Diffraction.hpp.
#include "pch.h"
#include "Render/Diffraction.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "Color/ChromaticAdaptation.hpp"
#include "Color/ColorSpace.hpp"
#include "Color/Illuminants.hpp"
#include "Render/Film.hpp"

namespace rayt {

    namespace {

        constexpr double kPi = 3.14159265358979323846;

        /// A deterministic number in [0, 1) from integers (the kernel does not depend on the run's random seed).
        double Hash01(std::uint64_t a, std::uint64_t b, std::uint64_t c) {
            std::uint64_t h = a * 0x9E3779B97F4A7C15ull ^ (b + 0x632BE59BD9B4E019ull) * 0xC2B2AE3D27D4EB4Full ^ c * 0x165667B19E3779F9ull;
            h ^= h >> 33; h *= 0xFF51AFD7ED558CCDull; h ^= h >> 33; h *= 0xC4CEB9FE1A85EC53ull; h ^= h >> 33;
            return double(h >> 11) * (1.0 / 9007199254740992.0);
        }

        /// In-place radix-2 FFT of n (a power of two) values; inverse without the 1/n.
        void Fft(std::complex<double>* a, int n, bool inverse) {
            for (int i = 1, j = 0; i < n; ++i) {   // bit reversal
                int bit = n >> 1;
                for (; j & bit; bit >>= 1) j ^= bit;
                j ^= bit;
                if (i < j) std::swap(a[i], a[j]);
            }
            for (int len = 2; len <= n; len <<= 1) {
                const double angle = 2 * kPi / len * (inverse ? 1 : -1);
                const std::complex<double> wl(std::cos(angle), std::sin(angle));
                for (int i = 0; i < n; i += len) {
                    std::complex<double> w(1.0);
                    for (int k = 0; k < len / 2; ++k) {
                        const std::complex<double> u = a[i + k], v = a[i + k + len / 2] * w;
                        a[i + k] = u + v;
                        a[i + k + len / 2] = u - v;
                        w *= wl;
                    }
                }
            }
        }

        /// 2D FFT of an nx x ny array (row-major), rows then columns.
        void Fft2(std::vector<std::complex<double>>& a, int nx, int ny, bool inverse) {
#pragma omp parallel for schedule(static)
            for (int y = 0; y < ny; ++y) Fft(&a[static_cast<size_t>(y) * nx], nx, inverse);
#pragma omp parallel
            {
                std::vector<std::complex<double>> col(static_cast<size_t>(ny));
#pragma omp for schedule(static)
                for (int x = 0; x < nx; ++x) {
                    for (int y = 0; y < ny; ++y) col[static_cast<size_t>(y)] = a[static_cast<size_t>(y) * nx + x];
                    Fft(col.data(), ny, inverse);
                    for (int y = 0; y < ny; ++y) a[static_cast<size_t>(y) * nx + x] = col[static_cast<size_t>(y)];
                }
            }
        }

        int NextPow2(int n) { int p = 1; while (p < n) p <<= 1; return p; }

        /// CIE X, Y, Z times D65 at a wavelength (nm): the weights of the three kernels.
        void XyzWeights(double nm, double w[3]) {
            const color::CIE_XYZ c = color::evaluateCIE_XYZ_1931(Real(nm));
            const double d65 = color::evaluateD65(Real(nm));
            w[0] = double(c.X) * d65;
            w[1] = double(c.Y) * d65;
            w[2] = double(c.Z) * d65;
        }

        /// The aperture of unit diameter: a disc, or a regular polygon of the disc's area (circumradius sqrt(pi / (2 n sin(2 pi / n))), corners
        /// at rotation + 2 pi k / n, counter-clockwise), with its edges prepared once.
        struct Aperture {
            struct Edge { double len, tx, ty, nx, ny, cx, cy; };
            int blades = 0;
            std::vector<Edge> edges;

            Aperture(int n, double rotation) : blades(n >= 3 ? n : 0) {
                if (blades == 0) return;
                const double rc = std::sqrt(kPi / (2.0 * n * std::sin(2 * kPi / n)));
                for (int k = 0; k < n; ++k) {
                    const double a0 = rotation + 2 * kPi * k / n, a1 = rotation + 2 * kPi * (k + 1) / n;
                    const double x0 = rc * std::cos(a0), y0 = rc * std::sin(a0), x1 = rc * std::cos(a1), y1 = rc * std::sin(a1);
                    const double ex = x1 - x0, ey = y1 - y0, len = std::sqrt(ex * ex + ey * ey);
                    edges.push_back({ len, ex / len, ey / len, ey / len, -ex / len, 0.5 * (x0 + x1), 0.5 * (y0 + y1) });
                }
            }

            /// F(q) = integral over the aperture of exp(-i q.r) d^2 r.
            std::complex<double> transform(double qx, double qy) const {
                const double q2 = qx * qx + qy * qy, area = kPi / 4;
                if (blades == 0) {   // a disc of radius a = 1/2: 2 pi a^2 J1(q a) / (q a)
                    const double x = std::sqrt(q2) * 0.5;
                    return x < 1e-8 ? area : 2 * area * std::cyl_bessel_j(1.0, x) / x;
                }
                if (q2 < 1e-16) return area;
                // By the divergence theorem: F(q) = (i / |q|^2) sum over the edges of (q.n) L exp(-i q.c) sinc(q.t L / 2), with the outward
                // normal n, the direction t, the length L and the midpoint c of each edge.
                std::complex<double> sum(0.0);
                for (const Edge& e : edges) {
                    const double qt = (qx * e.tx + qy * e.ty) * e.len * 0.5;
                    const double sinc = std::fabs(qt) < 1e-8 ? 1.0 : std::sin(qt) / qt;
                    const double phase = -(qx * e.cx + qy * e.cy);
                    sum += (qx * e.nx + qy * e.ny) * e.len * sinc * std::complex<double>(std::cos(phase), std::sin(phase));
                }
                return std::complex<double>(0.0, 1.0) * sum / q2;
            }
        };

    } // namespace

    std::complex<double> ApertureTransform(int blades, double rotation, double qx, double qy) {
        return Aperture(blades, rotation).transform(qx, qy);
    }

    DiffractionKernels ComputeDiffractionKernels(const DiffractionSetup& s) {
        DiffractionKernels out;
        const int r = std::max(s.radius, 1), size = 2 * r + 1;
        out.radius = r;
        for (auto& k : out.xyz) k.assign(static_cast<size_t>(size) * size, 0.0);
        const double pitch = double(s.pixelPitch) / std::max(double(s.scale), 1e-6);   // the pixel the pattern is measured in (metres)
        const double fNumber = double(s.fNumber), area = kPi / 4;
        const Aperture aperture(s.blades, double(s.bladeRotation));

        // Each pixel: n x n stratified points, each with its own wavelength (a Latin hypercube over 380..780 nm), the pattern of that
        // wavelength at the point weighted with X, Y, Z x D65. The pixels round the centre, where the core lies, get 64 x 64 points.
#pragma omp parallel for schedule(dynamic, 4)
        for (int iy = -r; iy <= r; ++iy) {
            for (int ix = -r; ix <= r; ++ix) {
                // At the real size the core (1.22 lambda N, a fraction of a pixel) puts most of the light into the centre pixel: dense there
                const int ring = std::max(std::abs(ix), std::abs(iy));
                const int count = ring <= 1 ? 65536 : ring <= 3 ? 4096 : 256;
                double acc[3] = { 0, 0, 0 }, wsum[3] = { 0, 0, 0 };
                const std::uint64_t pix = static_cast<std::uint64_t>((iy + r) * size + (ix + r));
                // The points (position in the pixel, wavelength) of the R3 low-discrepancy sequence (Roberts 2018), shifted per pixel: spread
                // evenly in all three at once, so the narrow core is not sampled at wavelengths of its own (a lattice pairing of position
                // and wavelength made the kernel hold 102 % of the light at the real size).
                const double s0 = Hash01(pix, 0, 0), s1 = Hash01(pix, 0, 1), s2 = Hash01(pix, 0, 2);
                for (int i = 0; i < count; ++i) {
                    const double p0 = std::fmod(s0 + 0.8191725133961645 * (i + 1), 1.0);
                    const double p1 = std::fmod(s1 + 0.6710436067037893 * (i + 1), 1.0);
                    const double p2 = std::fmod(s2 + 0.5497004779019703 * (i + 1), 1.0);
                    const double x = (ix - 0.5 + p0) * pitch;    // to the right
                    const double y = -(iy - 0.5 + p1) * pitch;   // up (the kernel's rows go down)
                    const double nm = 380.0 + 400.0 * p2;
                    const double lambdaN = nm * 1e-9 * fNumber;
                    const std::complex<double> f = aperture.transform(2 * kPi * x / lambdaN, 2 * kPi * y / lambdaN);
                    const double density = std::norm(f) / area / (lambdaN * lambdaN);   // per square metre; integrates to 1
                    double w[3];
                    XyzWeights(nm, w);
                    for (int c = 0; c < 3; ++c) { acc[c] += w[c] * density; wsum[c] += w[c]; }
                }
                for (int c = 0; c < 3; ++c)
                    out.xyz[c][static_cast<size_t>(iy + r) * size + (ix + r)] = wsum[c] > 0 ? acc[c] / wsum[c] * pitch * pitch : 0.0;
            }
        }
        for (int c = 0; c < 3; ++c) {
            double sum = 0;
            for (const double v : out.xyz[c]) sum += v;
            if (c == 1) out.captured = sum;
            if (sum > 0) for (double& v : out.xyz[c]) v /= sum;
        }
        return out;
    }

    void ApplyDiffraction(Film& film, const DiffractionKernels& k) {
        const int w = film.width(), h = film.height(), r = k.radius, size = 2 * r + 1;
        const int px = NextPow2(w + 2 * r), py = NextPow2(h + 2 * r);   // no wrap-around into the picture
        // linear sRGB <-> XYZ, the exact inverse of the renderer's XYZ_to_sRGB
        color::Matrix3 toRgb;
        for (int c = 0; c < 3; ++c) {
            color::CIE_XYZ e; e.X = c == 0; e.Y = c == 1; e.Z = c == 2;
            toRgb[c] = color::XYZ_to_sRGB(e);
        }
        const color::Matrix3 toXyz = glm::inverse(toRgb);
        const Vector3* src = film.getData();
        std::vector<Vector3> result(static_cast<size_t>(w) * h, Vector3(0.0));
        std::vector<std::complex<double>> img, ker;
        for (int c = 0; c < 3; ++c) {
            img.assign(static_cast<size_t>(px) * py, 0.0);
            ker.assign(static_cast<size_t>(px) * py, 0.0);
            // the picture in XYZ, shifted by r, with its border pixels repeated r pixels outward
            for (int y = 0; y < h + 2 * r; ++y)
                for (int x = 0; x < w + 2 * r; ++x) {
                    const int sx = std::clamp(x - r, 0, w - 1), sy = std::clamp(y - r, 0, h - 1);
                    img[static_cast<size_t>(y) * px + x] = (toXyz * src[static_cast<size_t>(sy) * w + sx])[c];
                }
            // the kernel with its centre at (0, 0), wrapped
            for (int y = -r; y <= r; ++y)
                for (int x = -r; x <= r; ++x)
                    ker[static_cast<size_t>((y + py) % py) * px + (x + px) % px] = k.xyz[c][static_cast<size_t>(y + r) * size + (x + r)];
            Fft2(img, px, py, false);
            Fft2(ker, px, py, false);
            for (size_t i = 0; i < img.size(); ++i) img[i] *= ker[i];
            Fft2(img, px, py, true);
            const double norm = 1.0 / (double(px) * double(py));
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x)
                    result[static_cast<size_t>(y) * w + x][c] = Real(img[static_cast<size_t>(y + r) * px + (x + r)].real() * norm);
        }
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) film.setPixel(x, y, toRgb * result[static_cast<size_t>(y) * w + x]);
    }

} // namespace rayt
