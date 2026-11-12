// Gaussian splats (Geometry/GaussianSplats.hpp): the 3DGS PLY reader, the particle BVH, the ray-Gaussian response and the compositing.
#include "pch.h"
#include "Geometry/GaussianSplats.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>

#include "Core/Constants.hpp"
#include "Color/ColorSpace.hpp"
#include "IO/PlyLoader.hpp"

namespace rayt {

    GaussianSplats::GaussianSplats(const std::vector<SplatSource>& sources, int degreeInFile, const SplatParams& params) : m_params(params) {
        m_degree = std::clamp(std::min(params.degree, degreeInFile), 0, 3);
        m_coeffs = (m_degree + 1) * (m_degree + 1);
        const int restInFile = (degreeInFile + 1) * (degreeInFile + 1) - 1;
        m_rotT = glm::transpose(params.rotation);

        m_particles.reserve(sources.size());
        m_sh.reserve(sources.size() * static_cast<size_t>(3 * m_coeffs));
        std::vector<AABB> boxes;
        boxes.reserve(sources.size());
        for (const SplatSource& s : sources) {
            if (!(s.opacity > params.alphaMin)) continue;   // never visible
            Particle p;
            p.mean = Point3(params.rotation * (Vector3(s.mean) * params.scale) + params.translate);
            if (params.hasKeep && !inside(params.keep, p.mean)) continue;
            if (params.hasRemove && inside(params.remove, p.mean)) continue;
            const Vector3 scale = s.scale * params.scale;
            if (!(scale.x > 0 && scale.y > 0 && scale.z > 0) || !std::isfinite(scale.x + scale.y + scale.z)) continue;
            const glm::dmat3 R = params.rotation * quatToMatrix(s.rot);   // columns = the particle's axes in the world
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) p.rotT[i][j] = static_cast<float>(R[i][j]);   // row i of R^T = column i of R (glm: R[col][row])
            p.invScale = glm::vec<3, float>(float(1 / scale.x), float(1 / scale.y), float(1 / scale.z));
            p.opacity = static_cast<float>(std::min(s.opacity, Real(0.999)));
            // the box: alpha reaches alphaMin at the Mahalanobis radius r = sqrt(2 ln(opacity / alphaMin)) (times each axis' scale)
            const Real r = std::sqrt(Real(2) * std::log(s.opacity / params.alphaMin));
            Vector3 half(0.0);
            for (int i = 0; i < 3; ++i) {
                Real sum = 0;
                for (int k = 0; k < 3; ++k) { const Real e = R[k][i] * scale[k] * r; sum += e * e; }   // glm: R[k] is column k = axis k
                half[i] = std::sqrt(sum);
            }
            boxes.emplace_back(Vector3(p.mean) - half, Vector3(p.mean) + half);
            m_bounds = AABB::unite(m_bounds, boxes.back());
            // the coefficients: degree 0 (f_dc), then for k = 1.. the three channels (f_rest is stored channel-major in the file)
            for (int c = 0; c < 3; ++c) m_sh.push_back(c < static_cast<int>(s.sh.size()) ? s.sh[static_cast<size_t>(c)] : 0.f);
            for (int k = 1; k < m_coeffs; ++k)
                for (int c = 0; c < 3; ++c) {
                    const size_t idx = 3 + static_cast<size_t>(c) * static_cast<size_t>(restInFile) + static_cast<size_t>(k - 1);
                    m_sh.push_back(idx < s.sh.size() ? s.sh[idx] : 0.f);
                }
            m_particles.push_back(p);
        }
        buildBvh(boxes);
    }

    void GaussianSplats::gather(const Ray& r, Real tMax, Vector3& radiance, Real& transmittance) const {
        thread_local std::vector<Hit> hits;
        hits.clear();
        collect(r, tMax, hits);
        radiance = Vector3(0.0);
        transmittance = 1;
        if (hits.empty()) return;
        std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.t < b.t; });
        const Vector3 dirLocal = m_rotT * r.d;   // the harmonics live in the cloud's frame
        for (const Hit& h : hits) {
            radiance += Vector3(transmittance * h.alpha) * color(h.index, dirLocal);
            transmittance *= Real(1) - h.alpha;
            if (transmittance < kStopT) { transmittance = 0; break; }
        }
    }

    Real GaussianSplats::transmittance(const Ray& r, Real tMax) const {
        thread_local std::vector<Hit> hits;
        hits.clear();
        collect(r, tMax, hits);
        Real T = 1;
        for (const Hit& h : hits) {
            T *= Real(1) - h.alpha;
            if (T < kStopT) return 0;
        }
        return T;
    }

    void GaussianSplats::collectBruteForce(const Ray& r, Real tMax, std::vector<std::pair<Real, Real>>& tAlpha) const {
        tAlpha.clear();
        for (size_t i = 0; i < m_particles.size(); ++i) {
            Real t, alpha;
            if (respond(m_particles[i], r, tMax, t, alpha)) tAlpha.emplace_back(t, alpha);
        }
    }

    Vector3 GaussianSplats::colorForDirection(size_t index, const Vector3& dirWorld) const { return color(static_cast<int>(index), m_rotT * dirWorld); }

    // ------------------------------------------------------------------------------------------------------------------------------
    // the file
    // ------------------------------------------------------------------------------------------------------------------------------
    namespace {
        double decodePly(io::plydetail::Type t, const unsigned char* buf) {
            using io::plydetail::Type;
            switch (t) {
            case Type::Char: { std::int8_t x; std::memcpy(&x, buf, 1); return x; }
            case Type::UChar: { std::uint8_t x; std::memcpy(&x, buf, 1); return x; }
            case Type::Short: { std::int16_t x; std::memcpy(&x, buf, 2); return x; }
            case Type::UShort: { std::uint16_t x; std::memcpy(&x, buf, 2); return x; }
            case Type::Int: { std::int32_t x; std::memcpy(&x, buf, 4); return x; }
            case Type::UInt: { std::uint32_t x; std::memcpy(&x, buf, 4); return x; }
            case Type::Float: { float x; std::memcpy(&x, buf, 4); return x; }
            case Type::Double: { double x; std::memcpy(&x, buf, 8); return x; }
            default: return 0;
            }
        }
    }

    bool GaussianSplats::loadPly(const std::string& filename, std::vector<SplatSource>& out, int& degreeOut, std::string* outError) {
        using namespace io::plydetail;
        auto fail = [&](const std::string& what) { if (outError) *outError = what + " (" + filename + ")"; return false; };
        std::ifstream f(filename, std::ios::binary);
        if (!f) return fail("cannot open");
        bool ascii = false, big = false;
        std::vector<Element> elements;
        std::string err;
        if (!parseHeader(f, ascii, big, elements, nullptr, err)) return fail(err);
        const Element* vertex = nullptr;
        for (const Element& e : elements) if (e.name == "vertex") { vertex = &e; break; }
        if (!vertex) return fail("no vertex element");
        if (!elements.empty() && &elements.front() != vertex) return fail("the vertex element must come first");
        for (const Property& p : vertex->props) if (p.list) return fail("a list property in the vertex element");

        // the columns
        int ix = -1, iy = -1, iz = -1, iop = -1, is[3] = { -1, -1, -1 }, ir[4] = { -1, -1, -1, -1 }, idc[3] = { -1, -1, -1 };
        std::vector<int> irest;
        for (size_t k = 0; k < vertex->props.size(); ++k) {
            const std::string& n = vertex->props[k].name;
            const int i = static_cast<int>(k);
            if (n == "x") ix = i; else if (n == "y") iy = i; else if (n == "z") iz = i; else if (n == "opacity") iop = i;
            else if (n.rfind("scale_", 0) == 0) { const int j = std::stoi(n.substr(6)); if (j >= 0 && j < 3) is[j] = i; }
            else if (n.rfind("rot_", 0) == 0) { const int j = std::stoi(n.substr(4)); if (j >= 0 && j < 4) ir[j] = i; }
            else if (n.rfind("f_dc_", 0) == 0) { const int j = std::stoi(n.substr(5)); if (j >= 0 && j < 3) idc[j] = i; }
            else if (n.rfind("f_rest_", 0) == 0) { const size_t j = static_cast<size_t>(std::stoi(n.substr(7))); if (irest.size() <= j) irest.resize(j + 1, -1); irest[j] = i; }
        }
        if (ix < 0 || iy < 0 || iz < 0) return fail("no x y z");
        if (iop < 0 || is[0] < 0 || is[1] < 0 || is[2] < 0 || ir[0] < 0 || ir[1] < 0 || ir[2] < 0 || ir[3] < 0 || idc[0] < 0 || idc[1] < 0 || idc[2] < 0)
            return fail("not a 3D Gaussian Splatting point cloud (needs opacity, scale_0..2, rot_0..3, f_dc_0..2)");
        for (int i : irest) if (i < 0) return fail("a gap in the f_rest_ properties");
        const size_t nRest = irest.size();
        int degree = 0;
        if (nRest == 0) degree = 0; else if (nRest == 9) degree = 1; else if (nRest == 24) degree = 2; else if (nRest == 45) degree = 3;
        else return fail("f_rest has " + std::to_string(nRest) + " values (expected 0, 9, 24 or 45)");
        degreeOut = degree;

        // the records: binary rows are read in blocks and decoded in memory
        const size_t nProps = vertex->props.size();
        std::vector<double> vals(nProps);
        out.clear();
        out.reserve(vertex->count);
        size_t recordSize = 0;
        for (const Property& p : vertex->props) recordSize += sizeOf(p.type);
        std::vector<unsigned char> block;
        const size_t perBlock = 4096;
        Reader reader(f, ascii, big);
        for (size_t i = 0; i < vertex->count; ++i) {
            if (ascii) {
                for (size_t k = 0; k < nProps; ++k) if (!reader.read(vertex->props[k].type, vals[k])) return fail("truncated data");
            }
            else {
                const size_t inBlock = i % perBlock;
                if (inBlock == 0) {
                    const size_t n = std::min(perBlock, vertex->count - i);
                    block.resize(n * recordSize);
                    f.read(reinterpret_cast<char*>(block.data()), static_cast<std::streamsize>(block.size()));
                    if (!f) return fail("truncated data");
                }
                const unsigned char* rec = block.data() + inBlock * recordSize;
                for (size_t k = 0; k < nProps; ++k) {
                    const Type t = vertex->props[k].type;
                    unsigned char buf[8];
                    std::memcpy(buf, rec, sizeOf(t));
                    rec += sizeOf(t);
                    if (big) std::reverse(buf, buf + sizeOf(t));
                    vals[k] = decodePly(t, buf);
                }
            }
            SplatSource s;
            s.mean = Point3(vals[static_cast<size_t>(ix)], vals[static_cast<size_t>(iy)], vals[static_cast<size_t>(iz)]);
            s.opacity = Real(1) / (Real(1) + std::exp(-vals[static_cast<size_t>(iop)]));
            for (int k = 0; k < 3; ++k) s.scale[k] = std::exp(vals[static_cast<size_t>(is[k])]);
            glm::dvec4 q(vals[static_cast<size_t>(ir[0])], vals[static_cast<size_t>(ir[1])], vals[static_cast<size_t>(ir[2])], vals[static_cast<size_t>(ir[3])]);
            const Real len = glm::length(q);
            s.rot = len > 0 ? q / len : glm::dvec4(1, 0, 0, 0);
            s.sh.resize(3 + nRest);
            for (int k = 0; k < 3; ++k) s.sh[static_cast<size_t>(k)] = static_cast<float>(vals[static_cast<size_t>(idc[k])]);
            for (size_t k = 0; k < nRest; ++k) s.sh[3 + k] = static_cast<float>(vals[static_cast<size_t>(irest[k])]);
            out.push_back(std::move(s));
        }
        return true;
    }

    // ------------------------------------------------------------------------------------------------------------------------------
    // mathematics of a particle
    // ------------------------------------------------------------------------------------------------------------------------------
    glm::dmat3 GaussianSplats::quatToMatrix(const glm::dvec4& q) {
        const Real w = q[0], x = q[1], y = q[2], z = q[3];
        glm::dmat3 R;
        // glm is column-major: R[col][row]
        R[0][0] = 1 - 2 * (y * y + z * z); R[1][0] = 2 * (x * y - w * z);     R[2][0] = 2 * (x * z + w * y);
        R[0][1] = 2 * (x * y + w * z);     R[1][1] = 1 - 2 * (x * x + z * z); R[2][1] = 2 * (y * z - w * x);
        R[0][2] = 2 * (x * z - w * y);     R[1][2] = 2 * (y * z + w * x);     R[2][2] = 1 - 2 * (x * x + y * y);
        return R;
    }

    Vector3 GaussianSplats::evalSh(const float* sh, int degree, const Vector3& d) {
        constexpr Real C0 = 0.28209479177387814;
        constexpr Real C1 = 0.4886025119029199;
        constexpr Real C2[5] = { 1.0925484305920792, -1.0925484305920792, 0.31539156525252005, -1.0925484305920792, 0.5462742152960396 };
        constexpr Real C3[7] = { -0.5900435899266435, 2.890611442640554, -0.4570457994644658, 0.3731763325901154, -0.4570457994644658, 1.445305721320277, -0.5900435899266435 };
        auto coef = [&](int k) { return Vector3(sh[3 * k], sh[3 * k + 1], sh[3 * k + 2]); };
        Vector3 c = C0 * coef(0);
        if (degree >= 1) {
            const Real x = d.x, y = d.y, z = d.z;
            c += -C1 * y * coef(1) + C1 * z * coef(2) - C1 * x * coef(3);
            if (degree >= 2) {
                const Real xx = x * x, yy = y * y, zz = z * z, xy = x * y, yz = y * z, xz = x * z;
                c += C2[0] * xy * coef(4) + C2[1] * yz * coef(5) + C2[2] * (2 * zz - xx - yy) * coef(6) + C2[3] * xz * coef(7) + C2[4] * (xx - yy) * coef(8);
                if (degree >= 3) {
                    c += C3[0] * y * (3 * xx - yy) * coef(9) + C3[1] * xy * z * coef(10) + C3[2] * y * (4 * zz - xx - yy) * coef(11)
                        + C3[3] * z * (2 * zz - 3 * xx - 3 * yy) * coef(12) + C3[4] * x * (4 * zz - xx - yy) * coef(13)
                        + C3[5] * z * (xx - yy) * coef(14) + C3[6] * x * (xx - 3 * yy) * coef(15);
                }
            }
        }
        return glm::max(c + Vector3(0.5), Vector3(0.0));
    }

    bool GaussianSplats::inside(const AABB& b, const Point3& p) {
        return p.x >= b.min.x && p.x <= b.max.x && p.y >= b.min.y && p.y <= b.max.y && p.z >= b.min.z && p.z <= b.max.z;
    }

    /// The ray-Gaussian response: the point of the ray closest to the center in the covariance metric, and alpha there.
    bool GaussianSplats::respond(const Particle& p, const Ray& r, Real tMax, Real& tOut, Real& alphaOut) const {
        const Vector3 o = Vector3(r.o) - Vector3(p.mean);
        Real ol[3], dl[3];
        for (int i = 0; i < 3; ++i) {
            ol[i] = (p.rotT[i][0] * o.x + p.rotT[i][1] * o.y + p.rotT[i][2] * o.z) * p.invScale[i];
            dl[i] = (p.rotT[i][0] * r.d.x + p.rotT[i][1] * r.d.y + p.rotT[i][2] * r.d.z) * p.invScale[i];
        }
        const Real dd = dl[0] * dl[0] + dl[1] * dl[1] + dl[2] * dl[2];
        if (!(dd > 0)) return false;
        const Real t = -(ol[0] * dl[0] + ol[1] * dl[1] + ol[2] * dl[2]) / dd;
        if (t <= r.tMin || t >= tMax) return false;
        Real q = 0;
        for (int i = 0; i < 3; ++i) { const Real x = ol[i] + t * dl[i]; q += x * x; }
        const Real alpha = p.opacity * std::exp(-Real(0.5) * q);
        if (alpha < m_params.alphaMin) return false;
        tOut = t;
        alphaOut = alpha;
        return true;
    }

    Vector3 GaussianSplats::color(int index, const Vector3& dirLocal) const {
        Vector3 c = evalSh(&m_sh[static_cast<size_t>(index) * static_cast<size_t>(3 * m_coeffs)], m_degree, dirLocal);
        if (m_params.linearize) c = color::sRGBToLinear(glm::min(c, Vector3(1.0)));
        return c * m_params.intensity;
    }

    // ------------------------------------------------------------------------------------------------------------------------------
    // the BVH
    // ------------------------------------------------------------------------------------------------------------------------------
    void GaussianSplats::collect(const Ray& r, Real tMax, std::vector<Hit>& hits) const {
        if (m_nodes.empty()) return;
        const Vector3 invDir(Real(1) / r.d.x, Real(1) / r.d.y, Real(1) / r.d.z);
        thread_local std::vector<int> stack;
        stack.clear();
        stack.push_back(0);
        while (!stack.empty()) {
            const Node& n = m_nodes[static_cast<size_t>(stack.back())];
            stack.pop_back();
            if (!slab(n.bounds, r, invDir, tMax)) continue;
            if (n.count > 0) {
                for (int k = 0; k < n.count; ++k) {
                    const int idx = m_order[static_cast<size_t>(n.start + k)];
                    Real t, alpha;
                    if (respond(m_particles[static_cast<size_t>(idx)], r, tMax, t, alpha)) hits.push_back({ t, alpha, idx });
                }
            }
            else {
                stack.push_back(n.left);
                stack.push_back(n.right);
            }
        }
    }

    bool GaussianSplats::slab(const AABB& b, const Ray& r, const Vector3& invDir, Real tMax) {
        Real t0 = r.tMin, t1 = tMax;
        for (int i = 0; i < 3; ++i) {
            Real near = (b.min[i] - r.o[i]) * invDir[i];
            Real far = (b.max[i] - r.o[i]) * invDir[i];
            if (near > far) std::swap(near, far);
            t0 = std::max(t0, near);
            t1 = std::min(t1, far);
            if (t0 > t1) return false;
        }
        return true;
    }

    void GaussianSplats::buildBvh(const std::vector<AABB>& boxes) {
        m_order.resize(boxes.size());
        for (size_t i = 0; i < boxes.size(); ++i) m_order[i] = static_cast<int>(i);
        m_nodes.clear();
        if (boxes.empty()) return;
        m_nodes.reserve(2 * boxes.size() / kLeafSize + 16);
        buildNode(boxes, 0, static_cast<int>(boxes.size()));
    }

    int GaussianSplats::buildNode(const std::vector<AABB>& boxes, int start, int end) {
        const int me = static_cast<int>(m_nodes.size());
        m_nodes.emplace_back();
        AABB bounds, centroids;
        for (int i = start; i < end; ++i) {
            const AABB& b = boxes[static_cast<size_t>(m_order[static_cast<size_t>(i)])];
            bounds = AABB::unite(bounds, b);
            const Vector3 c = (b.min + b.max) * Real(0.5);
            centroids = AABB::unite(centroids, AABB(c, c));
        }
        m_nodes[static_cast<size_t>(me)].bounds = bounds;
        const int n = end - start;
        const Vector3 ext = centroids.max - centroids.min;
        const int axis = ext.x >= ext.y && ext.x >= ext.z ? 0 : (ext.y >= ext.z ? 1 : 2);
        if (n <= kLeafSize || !(ext[axis] > 0)) {
            m_nodes[static_cast<size_t>(me)].start = start;
            m_nodes[static_cast<size_t>(me)].count = n;
            return me;
        }
        const int mid = start + n / 2;
        std::nth_element(m_order.begin() + start, m_order.begin() + mid, m_order.begin() + end, [&](int a, int b) {
            const AABB& ba = boxes[static_cast<size_t>(a)];
            const AABB& bb = boxes[static_cast<size_t>(b)];
            return ba.min[axis] + ba.max[axis] < bb.min[axis] + bb.max[axis];
        });
        const int left = buildNode(boxes, start, mid);
        const int right = buildNode(boxes, mid, end);
        m_nodes[static_cast<size_t>(me)].left = left;
        m_nodes[static_cast<size_t>(me)].right = right;
        return me;
    }

} // namespace rayt
