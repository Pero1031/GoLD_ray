// Non-template parts of Microfacet/MultiBounce.hpp: projected area, visible normals from either hemisphere, Zhang's segment term.
#include "pch.h"

#include "Microfacet/MultiBounce.hpp"

namespace rayt::multibounce {

    Real projectedArea(const GGXDistribution& dist, const Vector3& w) {
        const Real lam = dist.lambda(w);   // symmetric in w.z
        return w.z > Real(0) ? (Real(1) + lam) * w.z : lam * (-w.z);
    }

    Vector3 sampleVisibleNormal(Real alphaX, Real alphaY, const Vector3& w, Real u1, Real u2) {
        // stretch to alpha = 1, where the microsurface's normals are those of a hemisphere and the visible ones form a spherical cap
        const Vector3 ws = glm::normalize(Vector3(alphaX * w.x, alphaY * w.y, w.z));
        const Real phi = Real(2) * constants::PI * u1;
        const Real z = (Real(1) - u2) * (Real(1) + ws.z) - ws.z;
        const Real s = math::safe_sqrt(Real(1) - z * z);
        const Vector3 h = Vector3(s * std::cos(phi), s * std::sin(phi), z) + ws;
        // unstretch
        const Vector3 m(alphaX * h.x, alphaY * h.y, std::max(h.z, Real(0)));
        const Real len = glm::length(m);
        return len > Real(0) ? m / len : Vector3(0, 0, 1);
    }

    Real zhangSegment(const Real* lam, const bool* up, int n) {
        // The term of a sub-path e[a..b] depends only on a and b (the paper proves it reciprocal, so a sub-path and its reverse
        // agree): a table over (a, b). Type 3 needs the sub-paths without the first and without the last bounce, which always
        // start downward and end upward again and are one bounce shorter.
        constexpr int kMax = kMaxBounces + 2;
        if (n < 1 || n + 1 > kMax) return 0;
        // seg[a][b] for the sub-paths that start downward and end upward, shortest first
        Real seg[kMax][kMax] = {};
        for (int len = 1; len <= n; ++len) {
            for (int a = 0; a + len <= n; ++a) {
                const int b = a + len;
                if (up[a] || !up[b]) continue;
                Real r = 1;
                if (len == 1) {
                    r = Real(1) / (Real(1) + lam[a] + lam[b]);
                }
                else if (!up[b - 1]) {           // type 1: turns upward at the last bounce
                    for (int i = a; i < b; ++i) r /= Real(1) + lam[i] + lam[b];
                }
                else if (up[a + 1]) {            // type 2: turns upward at the first bounce
                    for (int i = a + 1; i <= b; ++i) r /= Real(1) + lam[i] + lam[a];
                }
                else {                           // type 3: in between (forward and inverse sub-paths, Eq. 25)
                    r = (seg[a + 1][b] + seg[a][b - 1]) / ((Real(1) + lam[a]) * (Real(1) + lam[b]));
                }
                seg[a][b] = r;
            }
        }
        return seg[0][n];
    }

} // namespace rayt::multibounce
