/**
 * @file Materials/MixMaterial.hpp
 * @brief Two materials mixed over the surface by a mask: gilt leather (metal leaf on the ridges, leather in the grooves),
 *        worn paint, a metal inlay. Header only.
 *
 * At every hit the mask value m(uv) (the weight of material A) is used as a probability: the hit behaves entirely as A with
 * probability m, otherwise entirely as B, as pbrt-v4's MixMaterial does. The choice is a hash of the hit point and of the
 * direction toward the viewer, so every call for the same hit (eval for light samples, sample, pdf, the polarization) sees the
 * same material and multiple importance sampling stays consistent. On average the surface reflects m f_A + (1 - m) f_B, without
 * bias; where the mask is between 0 and 1 the choice adds some noise (a blend region a few texels wide costs little).
 *
 * The mask is the red channel of an image laid on the uv, read as the value written in the file (gray 128 = 0.5: the sRGB
 * decoding of the image loader is undone), or a constant weight. Emitting, subsurface and medium-boundary materials are refused
 * by the scene loader: those need more than eval/sample at the hit.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>

#include "Core/Texture.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    class MixMaterial : public Material {
    public:
        /// `mask` (may be null): weight of A at the uv, from the red channel; without a mask the weight is `weight` everywhere.
        MixMaterial(std::shared_ptr<const Material> a, std::shared_ptr<const Material> b, std::shared_ptr<const Texture2D> mask, Real weight = Real(0.5))
            : m_a(std::move(a)), m_b(std::move(b)), m_mask(std::move(mask)), m_weight(std::clamp(weight, Real(0), Real(1))) {}

        /// Weight of material A at this hit, in [0, 1].
        Real weightA(const SurfaceInteraction& rec) const {
            if (!m_mask) return m_weight;
            const Real lin = std::clamp(m_mask->rgb(rec.uv.x, rec.uv.y).x, Real(0), Real(1));
            // back to the value stored in the file (the loader decoded sRGB): a painted gray of 128 means one half
            const Real v = lin <= Real(0.0031308) ? Real(12.92) * lin : Real(1.055) * std::pow(lin, Real(1) / Real(2.4)) - Real(0.055);
            return std::clamp(v, Real(0), Real(1));
        }

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, TransportMode mode = TransportMode::Radiance) const override {
            return pick(rec, wo).eval(rec, wo, wi, lambda, mode);
        }

        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
            const SampledWavelengths& lambda, TransportMode mode = TransportMode::Radiance) const override {
            return pick(rec, wo).sample(rec, wo, u, lambda, mode);
        }

        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const override {
            return pick(rec, wo).pdf(rec, wo, wi, lambda);
        }

        void prepareWavelengths(SampledWavelengths& lambda) const override {
            // the hit's material is not known here: a dispersive material on either side drops the secondary wavelengths
            m_a->prepareWavelengths(lambda);
            m_b->prepareWavelengths(lambda);
        }

        void perturbShading(SurfaceInteraction& rec, const Vector3& wo) const override {
            pick(rec, wo).perturbShading(rec, wo);
        }

        bool scatterMueller(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, pol::ScatterEvent& out) const override {
            return pick(rec, wo).scatterMueller(rec, wo, wi, lambda, out);
        }

        /// Specular only when both sides are (the integrators skip light sampling on specular surfaces).
        bool isSpecular() const override { return m_a->isSpecular() && m_b->isSpecular(); }

        Real opacity(const SurfaceInteraction& rec) const override {
            const Real m = weightA(rec);
            return m * m_a->opacity(rec) + (Real(1) - m) * m_b->opacity(rec);
        }

        Vector3 aovAlbedo(const SurfaceInteraction& rec) const override {
            const Real m = weightA(rec);
            return m * m_a->aovAlbedo(rec) + (Real(1) - m) * m_b->aovAlbedo(rec);
        }

        const Material& materialA() const { return *m_a; }
        const Material& materialB() const { return *m_b; }

        /// The material this hit behaves as (exposed for the checks).
        const Material& pick(const SurfaceInteraction& rec, const Vector3& wo) const {
            const Real m = weightA(rec);
            if (m >= Real(1)) return *m_a;
            if (m <= Real(0)) return *m_b;
            return hashUniform(rec.p, wo) < m ? *m_a : *m_b;
        }

    private:
        static std::uint64_t bits(Real x) {
            const double d = double(x);
            std::uint64_t b;
            std::memcpy(&b, &d, sizeof b);
            return b;
        }
        static std::uint64_t mix64(std::uint64_t z) {   // splitmix64 finalizer
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        }
        /// A uniform number in [0, 1) fixed by the hit point and the direction toward the viewer.
        static Real hashUniform(const Point3& p, const Vector3& wo) {
            std::uint64_t h = 0x9e3779b97f4a7c15ULL;
            for (const Real v : { p.x, p.y, p.z, wo.x, wo.y, wo.z }) h = mix64(h ^ bits(v));
            return Real(double(h >> 11) * (1.0 / 9007199254740992.0));
        }

        std::shared_ptr<const Material> m_a, m_b;
        std::shared_ptr<const Texture2D> m_mask;
        Real m_weight = Real(0.5);
    };

} // namespace rayt
