// Polarization checks (--scene polarization): Fresnel amplitudes against closed forms, the chain of Mueller matrices and basis rotations
// against a direct propagation of the electric field in 3D, and the Stokes vector of a path tracer ray reflected by smooth glass.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/SceneParts.hpp"
#include "Color/SpectralConversion.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Polarization.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Sphere.hpp"
#include "Lights/AreaLight.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/RoughConductor.hpp"
#include "Materials/ThinFilm.hpp"
#include "Render/Integrator.hpp"

#include <functional>

namespace rayt::cli {

namespace {

using pol::Complex;
constexpr Real kPi = Real(3.14159265358979323846);

struct CVec3 { Complex v[3]; };   // an electric field: three complex components in world space

Complex dotC(const CVec3& e, const Vector3& a) { return e.v[0] * a.x + e.v[1] * a.y + e.v[2] * a.z; }
CVec3 scaled(const Vector3& a, Complex c) { return { { c * a.x, c * a.y, c * a.z } }; }
CVec3 sum(const CVec3& a, const CVec3& b) { return { { a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2] } }; }

/// The Stokes vector of a field in the basis (x, y = k x x): the definitions of Core/Polarization.hpp.
void stokesOf(const CVec3& e, const Vector3& x, const Vector3& k, Real out[4]) {
    const Complex ex = dotC(e, x), ey = dotC(e, glm::cross(k, x));
    const Complex w = ex * std::conj(ey);
    out[0] = std::norm(ex) + std::norm(ey);
    out[1] = std::norm(ex) - std::norm(ey);
    out[2] = Real(2) * std::real(w);
    out[3] = -Real(2) * std::imag(w);
}

Vector3 randomUnit() {
    for (;;) {
        const Vector3 v(Real(2) * sampling::Random() - 1, Real(2) * sampling::Random() - 1, Real(2) * sampling::Random() - 1);
        const Real l = glm::length(v);
        if (l > Real(0.1) && l <= Real(1)) return v / l;
    }
}

Vector3 perpendicularTo(const Vector3& k) {
    const Vector3 r = randomUnit();
    return glm::normalize(r - k * glm::dot(r, k));
}

Real mean3(const Vector3& v) { return (v.x + v.y + v.z) / Real(3); }

} // namespace

int RunPolarizationTests() {
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok, const std::string& detail = "") {
        std::cout << "[Polarization]   " << name << (detail.empty() ? "" : "  " + detail) << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    std::cout << "[Polarization] Stokes vectors, Mueller matrices and the Fresnel amplitudes\n";

    // (1) The amplitudes against the closed forms of a dielectric: normal incidence, Brewster's angle, the unpolarized reflectance,
    //     energy conservation, total internal reflection.
    {
        const Real n = 1.5;
        const pol::Amplitudes a0 = pol::fresnelAmplitudes(1.0, Complex(n));
        report("normal incidence: rs = (1 - n) / (1 + n) = -rp",
            std::abs(a0.rs - Complex((1 - n) / (1 + n))) < 1e-12 && std::abs(a0.rp + a0.rs) < 1e-12);

        const Real thetaB = std::atan(n);
        const pol::Amplitudes ab = pol::fresnelAmplitudes(std::cos(thetaB), Complex(n));
        report("Brewster's angle atan(n): rp = 0", std::abs(ab.rp) < 1e-12, "|rp| = " + std::to_string(std::abs(ab.rp)));

        Real worstF = 0, worstE = 0;
        for (Real deg = 0; deg < 89.5; deg += 3.5) {
            const Real ci = std::cos(deg * kPi / 180);
            const pol::Amplitudes a = pol::fresnelAmplitudes(ci, Complex(n));
            worstF = std::max(worstF, std::abs(Real(0.5) * (std::norm(a.rs) + std::norm(a.rp)) - fresnel::fresnelDielectric(ci, 1.0, n)));
            const Complex ct = std::sqrt(Complex(1) - Complex(1 - ci * ci) / Complex(n * n));
            const Real g = std::real(Complex(n) * ct) / ci;   // the energy factor of transmission
            worstE = std::max(worstE, std::abs(std::norm(a.rs) + g * std::norm(a.ts) - 1));
            worstE = std::max(worstE, std::abs(std::norm(a.rp) + g * std::norm(a.tp) - 1));
        }
        report("unpolarized reflectance = fresnelDielectric, 26 angles", worstF < 1e-12, "max difference " + std::to_string(worstF));
        report("R + T = 1 for s and p, 26 angles", worstE < 1e-12, "max difference " + std::to_string(worstE));

        // Inside the glass beyond the critical angle (41.8 deg for n = 1.5): all the light is reflected, with a phase difference.
        const pol::Amplitudes t = pol::fresnelAmplitudes(std::cos(55 * kPi / 180), Complex(1 / n));
        const Complex c = t.rs * std::conj(t.rp);
        report("total internal reflection (55 deg inside): |rs| = |rp| = 1, phase difference", std::abs(std::abs(t.rs) - 1) < 1e-12 && std::abs(std::abs(t.rp) - 1) < 1e-12 && std::abs(std::imag(c)) > 0.1,
            "arg(rs conj(rp)) = " + std::to_string(std::arg(c) * 180 / kPi) + " deg");

        // A conductor: the unpolarized reflectance of the amplitudes is the engine's conductor Fresnel.
        Real worstM = 0;
        for (Real deg = 0; deg < 89.5; deg += 7) {
            const Real ci = std::cos(deg * kPi / 180);
            const pol::Amplitudes a = pol::fresnelAmplitudes(ci, Complex(0.27, 2.95));
            worstM = std::max(worstM, std::abs(Real(0.5) * (std::norm(a.rs) + std::norm(a.rp)) - fresnel::fresnelConductor(ci, Real(0.27), Real(2.95))));
        }
        report("metal (0.27 + 2.95 i): unpolarized reflectance = fresnelConductor, 13 angles", worstM < 1e-9, "max difference " + std::to_string(worstM));
    }

    // (2) A chain of two reflections, the second on a metal, with random directions: the Stokes vector at the camera from the matrix chain
    //     (the integrator's two lines) against the electric field propagated with the Fresnel amplitudes, vector by vector.
    {
        Real worst = 0;
        int trials = 0;
        for (int trial = 0; trial < 400; ++trial) {
            Vector3 n1 = randomUnit(), n2 = randomUnit();
            const Vector3 d0 = randomUnit();
            const Vector3 wo1 = -d0;
            if (glm::dot(wo1, n1) < 0) n1 = -n1;
            const Vector3 wi1 = Real(2) * glm::dot(wo1, n1) * n1 - wo1;
            const Vector3 wo2 = -wi1;
            if (glm::dot(wo2, n2) < 0) n2 = -n2;
            const Vector3 wi2 = Real(2) * glm::dot(wo2, n2) * n2 - wo2;
            const Real c1 = glm::dot(wo1, n1), c2 = glm::dot(wo2, n2);
            if (c1 < 0.05 || c2 < 0.05) continue;   // a grazing bounce: the s axis is ill-conditioned
            ++trials;

            const pol::Amplitudes a1 = pol::fresnelAmplitudes(c1, Complex(1.5));
            const pol::Amplitudes a2 = pol::fresnelAmplitudes(c2, Complex(0.3, 2.8));
            const Vector3 s1 = pol::planeNormalS(wo1, n1), s2 = pol::planeNormalS(wo2, n2);

            // The light that arrives at the second surface, along -wi2, with a random (elliptical) polarization.
            const Vector3 kIn2 = -wi2;
            const Vector3 ea = perpendicularTo(kIn2), eb = glm::cross(kIn2, ea);
            const Complex ca(sampling::Random() - Real(0.5), sampling::Random() - Real(0.5));
            const Complex cb(sampling::Random() - Real(0.5), sampling::Random() - Real(0.5));
            const CVec3 e0 = sum(scaled(ea, ca), scaled(eb, cb));

            // The field through the two reflections: s and p components with their own signs (p = k x s on each side).
            const Vector3 pIn2 = glm::cross(kIn2, s2), pOut2 = glm::cross(wo2, s2);
            const CVec3 e1 = sum(scaled(s2, a2.rs * dotC(e0, s2)), scaled(pOut2, a2.rp * dotC(e0, pIn2)));
            const Vector3 kIn1 = -wi1;
            const Vector3 pIn1 = glm::cross(kIn1, s1), pOut1 = glm::cross(wo1, s1);
            const CVec3 e2 = sum(scaled(s1, a1.rs * dotC(e1, s1)), scaled(pOut1, a1.rp * dotC(e1, pIn1)));

            const Vector3 k0 = wo1;
            const Vector3 xCam = perpendicularTo(k0);
            Real sCam[4], sIn[4];
            stokesOf(e2, xCam, k0, sCam);
            stokesOf(e0, s2, kIn2, sIn);   // in the basis of the incoming light at the last event

            // The chain, as updatePolarization() builds it.
            pol::Mueller M1 = pol::interfaceMueller(a1.rs, a1.rp);
            pol::Mueller T = pol::Mueller::identity();
            Vector3 xPrev = xCam;
            T = T * pol::basisChange(s1, xPrev, wo1) * M1;
            xPrev = s1;
            // The second event is the metal: its matrix uses one amplitude pair for every component here.
            pol::Mueller M2 = pol::interfaceMueller(a2.rs, a2.rp);
            T = T * pol::basisChange(s2, xPrev, wo2) * M2;
            const Real gain = Real(0.25) * (std::norm(a1.rs) + std::norm(a1.rp)) * (std::norm(a2.rs) + std::norm(a2.rp));   // M00 of both
            Real got[4];
            for (int i = 0; i < 4; ++i) {
                Real v = 0;
                for (int j = 0; j < 4; ++j) v += T.e[i][j][0] * sIn[j];
                got[i] = v * gain;
            }
            for (int i = 0; i < 4; ++i) worst = std::max(worst, std::abs(got[i] - sCam[i]) / std::max(Real(1e-9), sIn[0]));
        }
        report("two reflections (glass, metal) with random directions: Mueller chain = field propagation, " + std::to_string(trials) + " cases",
            worst < 1e-9, "max difference " + std::to_string(worst) + " of the input intensity");
    }

    // (2b) A thin film: the stack amplitudes against the plain interface (d = 0, a dielectric and a metal substrate), and energy conservation
    //      for a film on glass (|r|^2 + g |t|^2 = 1 for s and p, g = n2 cos2 / (n0 cos0)).
    {
        Real worst0 = 0, worstE = 0;
        for (Real deg = 0; deg < 89.5; deg += 4.5) {
            const Real ci = std::cos(deg * kPi / 180);
            const Real n0 = 1.0, nF = 1.4, n2 = 1.5;
            const thinfilm::StackAmplitudes a0 = thinfilm::amplitudes(ci, n0, nF, 0.0, Complex(n2), 550.0);
            const pol::Amplitudes f = pol::fresnelAmplitudes(ci, Complex(n2));
            worst0 = std::max({ worst0, std::abs(a0.rs - f.rs), std::abs(a0.rp - f.rp), std::abs(a0.ts - f.ts), std::abs(a0.tp - f.tp) });
            const thinfilm::StackAmplitudes m0 = thinfilm::amplitudes(ci, n0, nF, 0.0, Complex(0.27, 2.95), 550.0);
            const pol::Amplitudes fm = pol::fresnelAmplitudes(ci, Complex(0.27, 2.95));
            worst0 = std::max({ worst0, std::abs(m0.rs - fm.rs), std::abs(m0.rp - fm.rp) });
            for (const Real d : { Real(120), Real(310), Real(777) }) {
                const thinfilm::StackAmplitudes a = thinfilm::amplitudes(ci, n0, nF, d, Complex(n2), 550.0);
                const Real sin2 = (1 - ci * ci) / (n2 * n2);
                const Real g = n2 * std::sqrt(1 - sin2) / (n0 * ci);
                worstE = std::max({ worstE, std::abs(std::norm(a.rs) + g * std::norm(a.ts) - 1), std::abs(std::norm(a.rp) + g * std::norm(a.tp) - 1) });
            }
        }
        report("film amplitudes with d = 0 = the plain interface (dielectric and metal substrate), 20 angles", worst0 < 1e-12, "max difference " + std::to_string(worst0));
        report("film on glass (120, 310, 777 nm): |r|^2 + g |t|^2 = 1 for s and p, 20 angles", worstE < 1e-12, "max difference " + std::to_string(worstE));
    }

    // (2c) The film material with d = 0 on a glass substrate is plain glass: the same Mueller matrices for reflection and transmission, seen
    //      from the outside and from the inside (the film material in both build modes: per wavelength, or integrated per color channel).
    {
        FilmThickness none;
        none.base = 0;
        ThinFilm film(1.4, none, 1.5);
        Dielectric glass(1.5, 0.0);
        const SampledWavelengths lambda = [] {
#if defined(RAYT_SPECTRAL_MODE)
            return SampledWavelengths::sampleVisible(Real(0.37));
#else
            return SampledWavelengths();
#endif
        }();
        Real worst = 0;
        int cases = 0;
        for (const bool outside : { true, false })
            for (const Real deg : { Real(15), Real(35) }) {
                const Real th = deg * kPi / 180;
                const Real nWo = outside ? 1.0 : 1.5, nOther = outside ? 1.5 : 1.0;
                const Real thOther = std::asin(nWo / nOther * std::sin(th));
                SurfaceInteraction rec;
                rec.p = Point3(0.0);
                rec.n = rec.gn = Vector3(0, 1, 0);
                rec.frontFace = outside;
                const Vector3 wo(std::sin(th), std::cos(th), 0);
                for (const bool refl : { true, false }) {
                    const Vector3 wi = refl ? Vector3(-std::sin(th), std::cos(th), 0) : Vector3(-std::sin(thOther), -std::cos(thOther), 0);
                    pol::ScatterEvent a, b;
                    if (!film.scatterMueller(rec, wo, wi, lambda, a) || !glass.scatterMueller(rec, wo, wi, lambda, b)) { worst = 1; continue; }
                    ++cases;
                    for (int r = 0; r < 4; ++r)
                        for (int c = 0; c < 4; ++c)
                            for (int i = 0; i < pol::kN; ++i) worst = std::max(worst, std::abs(a.m.e[r][c][i] - b.m.e[r][c][i]));
                    worst = std::max(worst, glm::length(a.s - b.s));
                    Real w = 0;
                    for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) for (int i = 0; i < pol::kN; ++i) w = std::max(w, std::abs(a.m.e[r][c][i] - b.m.e[r][c][i]));
                    if (w > 1e-9) std::cout << "[Polarization]     (" << (outside ? "outside" : "inside") << ", " << deg << " deg, " << (refl ? "reflection" : "transmission") << ": difference " << w << ")" << std::endl;
                }
            }
#if defined(RAYT_SPECTRAL_MODE)
        const Real tol = Real(1e-9);
#else
        const Real tol = Real(2e-2);   // the RGB film is tabulated over cos(theta) (64 nodes): the error is largest near the critical angle (inside, 35 deg)
#endif
        report("film (d = 0) on glass = glass: Mueller matrices of reflection and transmission, " + std::to_string(cases) + " cases", worst < tol && cases == 8,
            "max difference " + std::to_string(worst));
    }

    // (3) The path tracer: a camera ray reflected by a smooth surface under a uniform sky. Only the reflection counts (depth 2: a refracted
    //     ray is not followed), so the Stokes vector is (|rs|^2 + |rp|^2, |rs|^2 - |rp|^2, 0, 0) / 2 in the basis (s, p) of the reflection:
    //     horizontal polarization for a horizontal surface seen from the side. The image reference direction is turned from s by 0, 30, 90 degrees.
    //     Glass is compared with the closed forms of the Fresnel equations, a metal (a narrow GGX lobe, light sampled by both strategies) with
    //     the amplitudes that part (1) and (2) have checked.
    auto reflectionCase = [&](const std::string& label, const std::shared_ptr<Material>& material, Real deg, int maxDepth, int N,
        const std::function<void(Real, Real&, Real&)>& reflectances, Real relTol) {
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0, material));
        auto scene = std::make_shared<Scene>(world);
        auto sky = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
        scene->finalize("spatial");

        const Real ti = deg * kPi / 180;
        Real rs2 = 0, rp2 = 0;
        reflectances(std::cos(ti), rs2, rp2);
        const Vector3 d(std::sin(ti), -std::cos(ti), 0);   // plane of incidence: x-y, so s is the z axis
        const Ray ray(Point3(0, 0, 0) - d * Real(10), d);
        const Vector3 k = -d;
        for (const Real turnDeg : { Real(0), Real(30), Real(90) }) {
            const Real turn = turnDeg * kPi / 180;
            const Vector3 sAxis(0, 0, 1), u = glm::cross(k, sAxis);
            const Vector3 xRef = std::cos(turn) * sAxis + std::sin(turn) * u;   // the reference direction of the image, turned by `turn` from s
            const Vector3 yRef = glm::cross(k, xRef);
            const Real alpha = std::atan2(glm::dot(sAxis, yRef), glm::dot(sAxis, xRef));   // angle of s in the image basis
            const Real e0 = Real(0.5) * (rs2 + rp2), e1 = Real(0.5) * (rs2 - rp2) * std::cos(2 * alpha), e2 = Real(0.5) * (rs2 - rp2) * std::sin(2 * alpha);

            PathIntegrator integrator(nullptr, maxDepth, 1);
            Real sum[4] = { 0, 0, 0, 0 }, sum2[4] = { 0, 0, 0, 0 };
            for (int i = 0; i < N; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
                SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
                SampledWavelengths lambda;
#endif
                PathStokes st;
                st.xRef = xRef;
                integrator.Li(ray, *scene, lambda, nullptr, &st);
                for (int c = 0; c < 4; ++c) {
                    const Real v = mean3(SpectrumToRGB(st.S[c], lambda));
                    sum[c] += v;
                    sum2[c] += v * v;
                }
            }
            const Real expected[4] = { e0, e1, e2, 0 };
            bool ok = true;
            std::string detail;
            for (int c = 0; c < 4; ++c) {
                const Real m = sum[c] / N, se = std::sqrt(std::max(Real(0), sum2[c] / N - m * m) / N);
                ok = ok && std::abs(m - expected[c]) <= 4 * se + relTol * e0;
                detail += (c ? ", " : "") + std::to_string(m) + " (" + std::to_string(expected[c]) + ")";
            }
            report(label + " at " + std::to_string(int(std::lround(deg))) + " deg, image axis turned " + std::to_string(int(turnDeg)) + " deg from s: S0..S3 (expected)", ok, detail);
        }
    };

    {
        const Real n = 1.5;
        auto glass = std::make_shared<Dielectric>(n, 0.0);
        for (const Real deg : { Real(56.3099), Real(40) })
            reflectionCase("glass", glass, deg, 2, 80000, [&](Real ci, Real& rs, Real& rp) {
                const Real ti = std::acos(ci), tt = std::asin(std::sin(ti) / n);
                rs = std::pow(std::sin(ti - tt) / std::sin(ti + tt), Real(2));
                rp = std::pow(std::tan(ti - tt) / std::tan(ti + tt), Real(2));
            }, Real(1e-3));
    }
    {
        auto table = std::make_shared<color::ComplexIorTable>();
        table->assign({ { 300.0, 0.27, 2.95 }, { 900.0, 0.27, 2.95 } });   // a gray metal (the same index at every wavelength)
        auto metal = std::make_shared<RoughConductor>(table, Real(0.02));
        for (const Real deg : { Real(70), Real(45) })
            reflectionCase("metal (narrow GGX lobe)", metal, deg, 2, 60000, [&](Real ci, Real& rs, Real& rp) {
                const pol::Amplitudes a = pol::fresnelAmplitudes(ci, Complex(0.27, 2.95));
                rs = std::norm(a.rs);
                rp = std::norm(a.rp);
            }, Real(0.03));   // the lobe has a width: the direction (and the Fresnel factor) is not exactly the mirror's
    }

    // (4) A polarized sky (EnvLight::setSkyPolarization, --sky-pol). The expected Stokes vectors come from the electric fields alone: the
    //     polarized part (degree P along e) and the unpolarized rest (two incoherent orthogonal fields of (1 - P) / 2 each) are reflected with
    //     the Fresnel amplitudes and projected on the image basis by stokesOf, without the Mueller matrices or basis rotations of the integrator.
    auto fieldStokes = [](const Vector3& e, const Vector3& kIn, Real P, const std::function<CVec3(const CVec3&)>& through,
        const Vector3& xOut, const Vector3& kOut, Real out[4]) {
        Real a[4], b[4];
        stokesOf(through(scaled(e, Complex(1))), xOut, kOut, a);
        stokesOf(through(scaled(glm::cross(kIn, e), Complex(1))), xOut, kOut, b);
        for (int i = 0; i < 4; ++i) out[i] = P * a[i] + Real(0.5) * (1 - P) * (a[i] + b[i]);
    };
    auto estimateStokes = [&](const Scene& scene, const Ray& ray, const Vector3& xRef, int maxDepth, int N, const Real expected[4], Real relTol,
        std::string& detail) {
        PathIntegrator integrator(nullptr, maxDepth, 1);
        Real sum[4] = { 0, 0, 0, 0 }, sum2[4] = { 0, 0, 0, 0 };
        for (int i = 0; i < N; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
            SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
            SampledWavelengths lambda;
#endif
            PathStokes st;
            st.xRef = xRef;
            integrator.Li(ray, scene, lambda, nullptr, &st);
            for (int c = 0; c < 4; ++c) {
                const Real v = mean3(SpectrumToRGB(st.S[c], lambda));
                sum[c] += v;
                sum2[c] += v * v;
            }
        }
        bool ok = true;
        for (int c = 0; c < 4; ++c) {
            const Real m = sum[c] / N, se = std::sqrt(std::max(Real(0), sum2[c] / N - m * m) / N);
            ok = ok && std::abs(m - expected[c]) <= 4 * se + relTol * expected[0];
            detail += (c ? ", " : "") + std::to_string(m) + " (" + std::to_string(expected[c]) + ")";
        }
        return ok;
    };
    {
        const Real P = 0.7;
        // The sky seen straight: its own Stokes vector, degree P (1 - c^2) / (1 + c^2) at angle acos(c) from the sun, in the image basis.
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0, std::make_shared<Dielectric>(1.5, 0.0)));   // below the rays
        auto scene = std::make_shared<Scene>(world);
        auto envLight = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
        scene->setEnvLight(envLight);
        scene->finalize("spatial");
        const Vector3 sun = glm::normalize(Vector3(1.0, 0.5, 0.2));
        envLight->setSkyPolarization(sun, P);
        for (const Vector3& dir : { Vector3(0.3, 0.8, -0.5), Vector3(-0.6, 0.3, 0.7), Vector3(0.1, 0.2, -0.9) }) {
            const Vector3 dv = glm::normalize(dir), kCam = -dv;
            const Vector3 xRef = glm::normalize(glm::cross(Vector3(0, 1, 0), kCam));
            const Real c = glm::dot(sun, dv), degree = P * (1 - c * c) / (1 + c * c);
            Real expected[4];
            fieldStokes(glm::normalize(glm::cross(sun, dv)), kCam, degree, [](const CVec3& e) { return e; }, xRef, kCam, expected);
            std::string detail;
            const bool ok = estimateStokes(*scene, Ray(Point3(0, 5, 0), dv), xRef, 2, 256, expected, Real(1e-6), detail);
            report("polarized sky seen straight, " + std::to_string(int(std::lround(std::acos(c) * 180 / kPi))) + " deg from the sun: S0..S3 (expected)", ok, detail);
        }
        // Nothing below the horizon is polarized, nor anything with P = 0.
        Vector3 x, s;
        const bool below = envLight->polarization(glm::normalize(Vector3(0.3, -0.2, 0.5)), x, s);
        envLight->setSkyPolarization(sun, 0.0);
        const bool off = envLight->polarization(glm::normalize(Vector3(0.3, 0.8, -0.5)), x, s);
        report("polarized sky: unpolarized below the horizon and with degree 0", !below && !off);
    }

    //     The reflection of the polarized sky: the sun at 90 deg from the mirror direction, so that the light is polarized along s, along p or
    //     at 45 deg between them. Glass reaches the sky by the escaping ray only; the metal's narrow lobe by both strategies (NEE and BSDF
    //     sampling with MIS), so a polarization factor that one of them missed would show.
    auto polarizedSkyCase = [&](const std::string& label, const std::shared_ptr<Material>& material, Real deg, int N,
        const std::function<pol::Amplitudes(Real)>& amplitudes, Real relTol) {
        const Real P = 0.7;
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0, material));
        auto scene = std::make_shared<Scene>(world);
        auto envLight = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
        scene->setEnvLight(envLight);
        scene->finalize("spatial");

        const Real ti = deg * kPi / 180;
        const Vector3 d(std::sin(ti), -std::cos(ti), 0), wo = -d;      // the camera ray, and the reflected light travelling to the camera
        const Vector3 r(std::sin(ti), std::cos(ti), 0), kIn = -r;      // the mirror direction (toward the sky) and the sky light along it
        const Vector3 sAxis(0, 0, 1), pIn = glm::cross(kIn, sAxis), pOut = glm::cross(wo, sAxis);
        const pol::Amplitudes a = amplitudes(std::cos(ti));
        auto reflect = [&](const CVec3& e) { return sum(scaled(sAxis, a.rs * dotC(e, sAxis)), scaled(pOut, a.rp * dotC(e, pIn))); };
        const Real turn = 30 * kPi / 180;
        const Vector3 xRef = std::cos(turn) * sAxis + std::sin(turn) * glm::cross(wo, sAxis);
        const Vector3 sunS(-std::cos(ti), std::sin(ti), 0), sunP(0, 0, 1);   // both at 90 deg from r: the field cross(sun, r) along s, along p
        const std::pair<const char*, Vector3> suns[] = { { "s", sunS }, { "p", sunP }, { "45 deg", glm::normalize(sunS + sunP) } };
        for (const auto& [name, sun] : suns) {
            envLight->setSkyPolarization(sun, P);
            Real expected[4];
            fieldStokes(glm::normalize(glm::cross(sun, r)), kIn, P, reflect, xRef, wo, expected);
            std::string detail;
            const bool ok = estimateStokes(*scene, Ray(Point3(0.0) - d * Real(10), d), xRef, 2, N, expected, relTol, detail);
            report(label + " at " + std::to_string(int(std::lround(deg))) + " deg under a sky polarized along " + name + ": S0..S3 (expected)", ok, detail);
        }
    };
    {
        auto glass = std::make_shared<Dielectric>(1.5, 0.0);
        for (const Real deg : { Real(56.3099), Real(40) })
            polarizedSkyCase("glass", glass, deg, 80000, [](Real ci) { return pol::fresnelAmplitudes(ci, Complex(1.5)); }, Real(1e-3));
        auto table = std::make_shared<color::ComplexIorTable>();
        table->assign({ { 300.0, 0.27, 2.95 }, { 900.0, 0.27, 2.95 } });
        auto metal = std::make_shared<RoughConductor>(table, Real(0.02));
        for (const Real deg : { Real(70), Real(45) })
            polarizedSkyCase("metal (narrow GGX lobe)", metal, deg, 60000, [](Real ci) { return pol::fresnelAmplitudes(ci, Complex(0.27, 2.95)); }, Real(0.03));
    }

    // (5) A polarized area light (AreaLight::setLinearPolarizer: a screen). Seen straight it has its own Stokes vector: along the normal the
    //     field is the transmission axis; off the normal, without refraction (film index 1), it is normal to the ray and to the absorbing axis.
    //     Then its reflection: a large quad square to the mirror direction, the axis along s, along p and at 45 deg, degree 0.7. Glass finds
    //     the light with the reflected ray only (addSurfaceEmission); the metal's narrow lobe by both strategies. Expected as in part (4).
    auto addScreen = [](Aggregate& world, Scene& scene, const Point3& center, const Vector3& normal, const Vector3& axis, Real half, Real degree, Real index) {
        const Vector3 u = glm::normalize(axis - glm::dot(axis, normal) * normal), v = glm::cross(normal, u);   // cross(u, v) = normal
        auto light = AddQuadLight(world, &scene, center - half * (u + v), Real(2) * half * u, Real(2) * half * v, EmissionSpectrum(RgbColor(1.0)), 1.0, true);
        static_cast<AreaLight&>(*light).setLinearPolarizer(axis, normal, degree, index);
        return light;
    };
    {
        const Vector3 n = glm::normalize(Vector3(0.2, 1.0, 0.1)), axis = glm::normalize(glm::cross(n, Vector3(0.3, 0.1, 1.0)));
        const Vector3 absorb = glm::cross(n, axis);
        for (const Real tiltDeg : { Real(0), Real(50) }) {
            auto world = std::make_shared<Aggregate>();
            auto scene = std::make_shared<Scene>(world);
            addScreen(*world, *scene, Point3(0, 0, 0), n, axis, 5.0, 0.7, 1.0);
            scene->finalize("spatial");
            const Real t = tiltDeg * kPi / 180;
            const Vector3 side = glm::normalize(axis + Real(0.6) * absorb);                 // tilt the view toward a direction between the two axes
            const Vector3 kCam = glm::normalize(std::cos(t) * n + std::sin(t) * side);      // the light's direction of travel, toward the camera
            const Vector3 xRef = glm::normalize(glm::cross(Vector3(0, 1, 0), kCam));
            const Vector3 e = tiltDeg == 0 ? axis : glm::normalize(glm::cross(kCam, absorb));
            Real expected[4];
            fieldStokes(e, kCam, 0.7, [](const CVec3& f) { return f; }, xRef, kCam, expected);
            std::string detail;
            const bool ok = estimateStokes(*scene, Ray(Point3(0.0) + kCam * Real(3), -kCam), xRef, 2, 256, expected, Real(1e-6), detail);
            report("polarized quad light seen " + std::to_string(int(tiltDeg)) + " deg off its normal: S0..S3 (expected)", ok, detail);
        }
        // With the sheet's index the field off the normal lies between the projected axis and the no-refraction one; degree 0 is unpolarized.
        AreaLight probe;
        Vector3 x, x1, xn, s;
        const Vector3 k = glm::normalize(n + Real(0.9) * glm::normalize(axis + absorb));
        probe.setLinearPolarizer(axis, n, 1.0, 1.0);
        probe.polarization(-k, x1, s);
        probe.setLinearPolarizer(axis, n, 1.0, 1.5);
        const bool on = probe.polarization(-k, xn, s);
        const Vector3 proj = glm::normalize(axis - glm::dot(axis, k) * k);
        const Real a1 = std::acos(std::min(Real(1), std::abs(glm::dot(x1, proj)))), an = std::acos(std::min(Real(1), std::abs(glm::dot(xn, proj))));
        probe.setLinearPolarizer(axis, n, 0.0);
        report("polarized quad light: the sheet's index turns the field less than no refraction; perpendicular to the ray; degree 0 is unpolarized",
            on && an > Real(1e-3) && an < a1 && std::abs(glm::dot(xn, k)) < Real(1e-9) && !probe.polarization(-k, x, s),
            "angle from the projected axis: " + std::to_string(an * 180 / kPi) + " deg (index 1.5), " + std::to_string(a1 * 180 / kPi) + " deg (index 1)");
    }
    auto polarizedScreenCase = [&](const std::string& label, const std::shared_ptr<Material>& material, Real deg, int N,
        const std::function<pol::Amplitudes(Real)>& amplitudes, Real relTol) {
        const Real P = 0.7;
        const Real ti = deg * kPi / 180;
        const Vector3 d(std::sin(ti), -std::cos(ti), 0), wo = -d;
        const Vector3 r(std::sin(ti), std::cos(ti), 0), kIn = -r;
        const Vector3 sAxis(0, 0, 1), pIn = glm::cross(kIn, sAxis), pOut = glm::cross(wo, sAxis);
        const pol::Amplitudes a = amplitudes(std::cos(ti));
        auto reflect = [&](const CVec3& e) { return sum(scaled(sAxis, a.rs * dotC(e, sAxis)), scaled(pOut, a.rp * dotC(e, pIn))); };
        const Real turn = 30 * kPi / 180;
        const Vector3 xRef = std::cos(turn) * sAxis + std::sin(turn) * glm::cross(wo, sAxis);
        const std::pair<const char*, Vector3> axes[] = { { "s", sAxis }, { "p", pIn }, { "45 deg", glm::normalize(sAxis + pIn) } };
        for (const auto& [name, axis] : axes) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0, material));
            auto scene = std::make_shared<Scene>(world);
            addScreen(*world, *scene, Point3(0.0) + r * Real(20), kIn, axis, 6.0, P, 1.5);   // faces the surface point, square to the mirror direction
            scene->finalize("spatial");
            Real expected[4];
            fieldStokes(axis, kIn, P, reflect, xRef, wo, expected);
            std::string detail;
            const bool ok = estimateStokes(*scene, Ray(Point3(0.0) - d * Real(10), d), xRef, 2, N, expected, relTol, detail);
            report(label + " at " + std::to_string(int(std::lround(deg))) + " deg under a quad light polarized along " + name + ": S0..S3 (expected)", ok, detail);
        }
    };
    {
        auto glass = std::make_shared<Dielectric>(1.5, 0.0);
        for (const Real deg : { Real(56.3099), Real(40) })
            polarizedScreenCase("glass", glass, deg, 80000, [](Real ci) { return pol::fresnelAmplitudes(ci, Complex(1.5)); }, Real(1e-3));
        auto table = std::make_shared<color::ComplexIorTable>();
        table->assign({ { 300.0, 0.27, 2.95 }, { 900.0, 0.27, 2.95 } });
        auto metal = std::make_shared<RoughConductor>(table, Real(0.02));
        for (const Real deg : { Real(70), Real(45) })
            polarizedScreenCase("metal (narrow GGX lobe)", metal, deg, 60000, [](Real ci) { return pol::fresnelAmplitudes(ci, Complex(0.27, 2.95)); }, Real(0.03));
    }

    std::cout << "[Polarization] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
