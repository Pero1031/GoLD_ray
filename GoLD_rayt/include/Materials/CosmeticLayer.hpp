#pragma once
// A layer of foundation make-up: a thin slab of a scattering medium with two kinds of particles, after
//   D. Lanza, J. Padron-Griffe, A. Pranovich, A. Munoz, J. R. Frisvad, A. Jarabo, "Practical Appearance Model for Foundation
//   Cosmetics", Computer Graphics Forum 43(4), 2024 (EGSR).
// Written from the paper (its code is not used). docs/COSMETICS.md says what is taken over and what differs.
//
//   - spherical diffusers (titanium dioxide, silica): a mixture of two Henyey-Greenstein lobes, the matte part of a finish;
//   - platelets (talc, kaolin, mica): mirror-like flakes whose normals follow an SGGX distribution, the glossy part.
//
// The slab is index-matched with the air (no interface, no refraction: the paper's fits gave an index close to 1), and light
// leaves where it entered (the layer is thin). Lengths are in mean free paths of the densest direction: the extinction is
//     sigma_t(w) = c_d + (1 - c_d) Chat_p(w),      Chat_p(w) = projected area of the platelets, 1 along their mean normal,
// so `thickness` is the optical thickness seen along the platelets' mean normal. Colour comes from absorption only: every
// collision multiplies the path by the albedo of the particle that was hit.
//
// Frame: z is the surface normal, the slab is 0 <= z <= thickness, the air is above (z > thickness), the skin below.
// Directions point the way the light travels.

#include <algorithm>
#include <cmath>

#include "Core/Constants.hpp"
#include "Core/Sampling.hpp"
#include "Core/Types.hpp"

namespace rayt {

class CosmeticLayer {
public:
    struct Params {
        Real thickness = Real(0.5);                 ///< optical thickness t (mean free paths along the platelets' normal)
        Real diffuserFraction = Real(0.9);          ///< c_d: share of the diffusers in the extinction, 0 (platelets only) to 1
        Vector3 diffuserAlbedo = Vector3(0.97, 0.94, 0.90);   ///< Lambda_d, single-scattering albedo of the diffusers (RGB)
        Vector3 plateletAlbedo = Vector3(0.85);     ///< Lambda_p, of the platelets
        Real g1 = Real(0.5);                        ///< mean cosine of the first Henyey-Greenstein lobe of the diffusers
        Real g2 = Real(-0.2);                       ///< of the second lobe
        Real lobeWeight = Real(0.7);                ///< w_g: weight of the first lobe
        Real plateletRoughness = Real(0.15);        ///< alpha_p: spread of the platelet normals (SGGX), small = aligned = sharp
        Real plateletTiltDeg = Real(0);             ///< theta_p: the mean platelet normal is tilted from the surface normal towards +y (degrees)
        int maxBounces = 512;                       ///< a path that scatters more often is dropped (the paper uses 512)
    };

    /// Null if the parameters are usable, else what is wrong with them.
    static const char* validate(const Params& p) {
        if (!(p.thickness >= 0)) return "thickness must be >= 0";
        if (!(p.diffuserFraction >= 0 && p.diffuserFraction <= 1)) return "diffuserFraction must be in [0, 1]";
        for (int c = 0; c < 3; ++c)
            if (!(p.diffuserAlbedo[c] >= 0 && p.diffuserAlbedo[c] <= 1 && p.plateletAlbedo[c] >= 0 && p.plateletAlbedo[c] <= 1)) return "the albedos must be in [0, 1]";
        if (!(std::abs(p.g1) < 1 && std::abs(p.g2) < 1)) return "g1 and g2 must be inside (-1, 1)";
        if (!(p.lobeWeight >= 0 && p.lobeWeight <= 1)) return "lobeWeight must be in [0, 1]";
        if (!(p.plateletRoughness >= Real(0.001) && p.plateletRoughness <= 1)) return "plateletRoughness must be in [0.001, 1]";
        if (!(std::abs(p.plateletTiltDeg) <= 89)) return "plateletTiltDeg must be within 89 degrees";
        if (p.maxBounces < 1) return "maxBounces must be >= 1";
        return nullptr;
    }

    explicit CosmeticLayer(const Params& p) : m_p(p) {
        // SGGX matrix S = R diag(a^2, a^2, 1) R^T of flakes whose mean normal m is the z axis turned about x by the tilt.
        const Real a2 = p.plateletRoughness * p.plateletRoughness;
        const Real th = p.plateletTiltDeg * constants::PI / Real(180);
        // (S = a^2 I + (1 - a^2) m m^T, so w^T S w and m^T S^-1 m have the closed forms used below.)
        // The paper does not say which way the normal tilts. Towards +y here: for light travelling towards +y the flakes lean
        // away from the light, and the highlight moves towards grazing. With the other sign the paper's "matte 2" (52 degrees)
        // would be brightest back towards the light, which its measurement does not show (docs/COSMETICS.md).
        m_mean = Vector3(0, std::sin(th), std::cos(th));
        m_a2 = a2;
    }

    const Params& params() const { return m_p; }

    // ---- the platelets (SGGX: Heitz, Dupuy, Crassin, Dachsbacher 2015) ------------------------------------------------

    /// Projected area of the platelets seen along w, normalised to 1 along their mean normal: sqrt(w^T S w).
    Real plateletArea(const Vector3& w) const {
        const Real c = glm::dot(w, m_mean);
        return std::sqrt(std::max(Real(0), m_a2 + (Real(1) - m_a2) * c * c));
    }

    /// Distribution of the platelet normals: D(m) = 1 / (pi sqrt|S| (m^T S^-1 m)^2).
    Real plateletD(const Vector3& n) const {
        const Real c = glm::dot(n, m_mean);
        const Real q = (Real(1) - c * c) / m_a2 + c * c;   // m^T S^-1 m, with S^-1 = I / a^2 + (1 - 1 / a^2) m m^T
        return Real(1) / (constants::PI * m_a2 * q * q);
    }

    /// Phase function of the platelets: mirror reflection off the flakes, without Fresnel (equation 10 of the paper).
    /// Light travelling along w leaves along wo; the flake normal is the half vector of -w and wo.
    Real plateletPhase(const Vector3& w, const Vector3& wo) const {
        const Vector3 h = wo - w;
        const Real len = glm::length(h);
        if (len < Real(1e-12)) return 0;   // straight on: no flake reflects this way
        return plateletD(h / len) / (Real(4) * plateletArea(w));
    }

    /// Samples the direction after a reflection off a platelet (the visible normals of SGGX): its density is plateletPhase().
    Vector3 samplePlatelet(const Vector3& w, Real u1, Real u2) const {
        const Vector3 wi = -w;   // towards where the light comes from
        // an orthonormal basis (wk, wj, wi)
        const Vector3 helper = std::abs(wi.z) < Real(0.9) ? Vector3(0, 0, 1) : Vector3(1, 0, 0);
        const Vector3 wk = glm::normalize(glm::cross(helper, wi));
        const Vector3 wj = glm::cross(wi, wk);
        auto quad = [&](const Vector3& a, const Vector3& b) {
            const Real am = glm::dot(a, m_mean), bm = glm::dot(b, m_mean);
            return m_a2 * glm::dot(a, b) + (Real(1) - m_a2) * am * bm;   // a^T S b
        };
        const Real Skk = quad(wk, wk), Sjj = quad(wj, wj), Sii = quad(wi, wi);
        const Real Skj = quad(wk, wj), Ski = quad(wk, wi), Sji = quad(wj, wi);
        const Real det = std::abs(Skk * Sjj * Sii - Skj * Skj * Sii - Ski * Ski * Sjj - Sji * Sji * Skk + Real(2) * Skj * Ski * Sji);
        const Real sqrtDet = std::sqrt(det);
        const Real invSqrtSii = Real(1) / std::sqrt(Sii);
        const Real tmp = std::sqrt(std::max(Real(1e-30), Sjj * Sii - Sji * Sji));
        const Vector3 Mk(sqrtDet / tmp, 0, 0);
        const Vector3 Mj(-invSqrtSii * (Ski * Sji - Skj * Sii) / tmp, invSqrtSii * tmp, 0);
        const Vector3 Mi(invSqrtSii * Ski, invSqrtSii * Sji, invSqrtSii * Sii);
        // a point of the unit half sphere, cosine distributed around wi
        const Real r = std::sqrt(u1), phi = constants::TWO_PI * u2;
        const Real u = r * std::cos(phi), v = r * std::sin(phi), ww = std::sqrt(std::max(Real(0), Real(1) - u * u - v * v));
        const Vector3 local = glm::normalize(u * Mk + v * Mj + ww * Mi);
        const Vector3 n = glm::normalize(local.x * wk + local.y * wj + local.z * wi);   // the flake normal
        return glm::normalize(w - Real(2) * glm::dot(w, n) * n);                         // mirror direction
    }

    // ---- the diffusers ---------------------------------------------------------------------------------------------

    static Real henyeyGreenstein(Real cosTheta, Real g) {
        const Real d = Real(1) + g * g - Real(2) * g * cosTheta;
        return (Real(1) - g * g) / (Real(4) * constants::PI * d * std::sqrt(d));
    }

    /// Phase function of the diffusers: two Henyey-Greenstein lobes (equation 9).
    Real diffuserPhase(const Vector3& w, const Vector3& wo) const {
        const Real c = glm::dot(w, wo);
        return m_p.lobeWeight * henyeyGreenstein(c, m_p.g1) + (Real(1) - m_p.lobeWeight) * henyeyGreenstein(c, m_p.g2);
    }

    Vector3 sampleDiffuser(const Vector3& w, Real uLobe, Real u1, Real u2) const {
        const Real g = uLobe < m_p.lobeWeight ? m_p.g1 : m_p.g2;
        Real c;
        if (std::abs(g) < Real(1e-4)) c = Real(1) - Real(2) * u1;
        else {
            const Real s = (Real(1) - g * g) / (Real(1) - g + Real(2) * g * u1);
            c = std::clamp((Real(1) + g * g - s * s) / (Real(2) * g), Real(-1), Real(1));
        }
        const Real sn = std::sqrt(std::max(Real(0), Real(1) - c * c)), phi = constants::TWO_PI * u2;
        const Vector3 helper = std::abs(w.z) < Real(0.9) ? Vector3(0, 0, 1) : Vector3(1, 0, 0);
        const Vector3 t = glm::normalize(glm::cross(helper, w));
        const Vector3 b = glm::cross(w, t);
        return glm::normalize(sn * std::cos(phi) * t + sn * std::sin(phi) * b + c * w);
    }

    // ---- the medium ------------------------------------------------------------------------------------------------

    /// Extinction for light travelling along w, in units of the base extinction (equation 5).
    Real sigmaT(const Vector3& w) const { return m_p.diffuserFraction + (Real(1) - m_p.diffuserFraction) * plateletArea(w); }

    /// What a collision of light travelling along w sends towards wo, per steradian: albedo times phase function of the
    /// mixture (equations 7 and 8 multiplied).
    Vector3 scatter(const Vector3& w, const Vector3& wo) const {
        const Real cd = m_p.diffuserFraction, cp = (Real(1) - cd) * plateletArea(w);
        if (cd + cp <= 0) return Vector3(0);
        Vector3 f(0);
        if (cd > 0) f += m_p.diffuserAlbedo * (cd * diffuserPhase(w, wo));
        if (cp > 0) f += m_p.plateletAlbedo * (cp * plateletPhase(w, wo));
        return f / (cd + cp);
    }

    /// One path through the layer.
    struct Walk {
        bool up = false;            ///< left through the top (reflected) or through the bottom (transmitted to the skin)
        bool lost = false;          ///< dropped after maxBounces collisions (its weight is not counted)
        int bounces = 0;            ///< collisions inside the layer (0: the light passed straight through)
        Vector3 direction{ 0 };     ///< direction of travel when it left
        Vector3 weight{ 0 };        ///< what the absorption left of it (RGB)
    };

    /// Follows light that enters the layer travelling along `w`: from the air (fromTop, w.z < 0) or from the skin (w.z > 0).
    Walk walk(const Vector3& w, bool fromTop = true) const {
        Walk out;
        Vector3 dir = w;
        Real z = fromTop ? m_p.thickness : Real(0);
        Vector3 weight(1.0);
        for (int bounce = 0; ; ++bounce) {
            const Real s = -std::log(std::max(Real(1e-300), Real(1) - sampling::Random())) / sigmaT(dir);
            z += s * dir.z;
            if (z >= m_p.thickness || z <= 0) {
                out.up = z >= m_p.thickness;
                out.bounces = bounce;
                out.direction = dir;
                out.weight = weight;
                return out;
            }
            if (bounce >= m_p.maxBounces) { out.lost = true; out.bounces = bounce; return out; }
            // which particle was hit: by its share of the extinction; its albedo is what survives
            const Real cd = m_p.diffuserFraction, cp = (Real(1) - cd) * plateletArea(dir);
            if (sampling::Random() * (cd + cp) < cd) {
                weight *= m_p.diffuserAlbedo;
                dir = sampleDiffuser(dir, sampling::Random(), sampling::Random(), sampling::Random());
            }
            else {
                weight *= m_p.plateletAlbedo;
                dir = samplePlatelet(dir, sampling::Random(), sampling::Random());
            }
            if (weight.x + weight.y + weight.z <= 0) { out.lost = true; out.bounces = bounce + 1; return out; }
        }
    }

    /// One sample of the reflection of the layer alone (nothing below it): light arrives travelling along `w` (w.z < 0) and
    /// is seen leaving along `wo` (wo.z > 0). The mean of many samples is f(w, wo) cos(theta_o), the BSDF times the cosine
    /// of the leaving direction. Every collision of the path is connected to wo (next event estimation, as in position-free
    /// Monte Carlo: Guo, Hasan, Zhao 2018).
    Vector3 reflectionSample(const Vector3& w, const Vector3& wo) const {
        Vector3 sum(0);
        Vector3 dir = w;
        Real z = m_p.thickness;
        Vector3 weight(1.0);
        const Real sigmaOut = sigmaT(wo);
        for (int bounce = 0; bounce < m_p.maxBounces; ++bounce) {
            const Real s = -std::log(std::max(Real(1e-300), Real(1) - sampling::Random())) / sigmaT(dir);
            z += s * dir.z;
            if (z >= m_p.thickness || z <= 0) break;
            sum += weight * scatter(dir, wo) * std::exp(-sigmaOut * (m_p.thickness - z) / wo.z);
            const Real cd = m_p.diffuserFraction, cp = (Real(1) - cd) * plateletArea(dir);
            if (sampling::Random() * (cd + cp) < cd) {
                weight *= m_p.diffuserAlbedo;
                dir = sampleDiffuser(dir, sampling::Random(), sampling::Random(), sampling::Random());
            }
            else {
                weight *= m_p.plateletAlbedo;
                dir = samplePlatelet(dir, sampling::Random(), sampling::Random());
            }
        }
        return sum;
    }

    /// The part of the light travelling along w that crosses the whole layer without a collision.
    Real unscattered(const Vector3& w) const { return std::exp(-sigmaT(w) * m_p.thickness / std::max(Real(1e-9), std::abs(w.z))); }

private:
    Params m_p;
    Vector3 m_mean{ 0, 0, 1 };
    Real m_a2 = Real(1);
};

} // namespace rayt
