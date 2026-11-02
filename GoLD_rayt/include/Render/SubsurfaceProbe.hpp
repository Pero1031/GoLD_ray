/**
 * @file Render/SubsurfaceProbe.hpp
 * @brief Samples the exit point of a subsurface path with probe rays (PBRT v3 style BSSRDF sampling).
 *
 * Given the entry point po the exit point pi is drawn in three steps:
 *  1. pick a radius r from the profile of one color channel (chosen uniformly), and an angle phi;
 *  2. cast a probe ray through the sphere of radius rMax around po, along the surface normal
 *     (probability 1/2) or one of the two tangents (1/4 each), at offset (r, phi): this finds
 *     surface points that are about r away even where the surface curves or folds;
 *  3. pick one of the probe hits that belong to the same SubsurfaceProfile.
 * The weight is Sp(|po - pi|) / pdf, where pdf is the density of choosing exactly that pi: the
 * sum over all three axes and all channels (so a point found "the other way" still gets the right
 * weight), divided by the number of candidates. The Fresnel factors are part of the entry event
 * (CoatedDiffuse) and of the exit material (SubsurfaceExit).
 */

#pragma once

#include <array>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Sampling.hpp"
#include "Core/Interaction.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Materials/Subsurface.hpp"
#include "Scene/Scene.hpp"

namespace rayt {

    /**
     * @brief Moves a path that entered a subsurface material at `po` to a sampled exit point.
     *
     * @param exitVertex receives the exit vertex: geometry of the hit, normal pointing out of the
     *                   surface, and the material set to the SubsurfaceExit of `po`'s material
     * @param weight     Sp / pdf per channel; multiply it onto the path throughput
     * @return false if no exit point was found (the path ends)
     */
    inline bool SampleSubsurfaceExit(const Scene& scene, const SurfaceInteraction& po,
        const SampledWavelengths& lambda, SurfaceInteraction& exitVertex, Spectrum& weight)
    {
        const Material* material = po.matPtr;
        const SubsurfaceProfile* profile = material ? material->subsurface() : nullptr;
        if (!profile) return false;

        SubsurfaceParams params;
        material->subsurfaceParams(po, lambda, params);

        const frame::Frame frame(po.n);
        const std::array<Vector3, 3> axes = { frame.s, frame.t, frame.n };
        constexpr std::array<Real, 3> axisProb = { Real(0.25), Real(0.25), Real(0.5) };

        // --- 1. radius and angle -------------------------------------------------------------
        const int channel = std::min(static_cast<int>(sampling::Random() * Real(kSssChannels)), kSssChannels - 1);
        const Real r = sss::sampleRadius(sampling::Random(), params.d[channel]);
        Real rMax = 0;
        for (int c = 0; c < kSssChannels; ++c) rMax = std::max(rMax, sss::sampleRadius(Real(0.999), params.d[c]));
        if (r >= rMax) return false;   // the 0.1 % tail that is cut off
        const Real phi = Real(2) * constants::PI * sampling::Random();

        // --- 2. probe ray along one of the three axes ----------------------------------------
        const Real uAxis = sampling::Random();
        const int a = uAxis < Real(0.5) ? 2 : (uAxis < Real(0.75) ? 0 : 1);
        const Vector3& vz = axes[static_cast<size_t>(a)];
        const Vector3& vx = axes[static_cast<size_t>((a + 1) % 3)];
        const Vector3& vy = axes[static_cast<size_t>((a + 2) % 3)];

        const Real l = Real(2) * std::sqrt(rMax * rMax - r * r);
        Ray ray(po.p + (vx * std::cos(phi) + vy * std::sin(phi)) * r - vz * (l * Real(0.5)), vz, Real(0));
        ray.tMax = l;
        ray.time = po.time;

        // --- 3. all hits of the same profile along the segment -------------------------------
        constexpr int kMaxCandidates = 8, kMaxHits = 32;
        std::array<SurfaceInteraction, kMaxCandidates> found;
        int nFound = 0;
        for (int hits = 0; hits < kMaxHits && ray.tMax > Real(0); ++hits) {
            SurfaceInteraction h;
            if (!scene.hit(ray, h)) break;
            if (h.matPtr && h.matPtr->subsurface() == profile && nFound < kMaxCandidates) found[static_cast<size_t>(nFound++)] = h;
            const Real advance = h.t + constants::RAY_EPSILON;   // continue behind this hit
            ray.o = ray.o + ray.d * advance;
            ray.tMax -= advance;
        }
        if (nFound == 0) return false;

        const int selected = std::min(static_cast<int>(sampling::Random() * Real(nFound)), nFound - 1);
        const SurfaceInteraction& pi = found[static_cast<size_t>(selected)];

        // --- density of this exit point over everything that could have produced it ----------
        const Vector3 dv = po.p - pi.p;
        const Vector3 dLocal(glm::dot(dv, frame.s), glm::dot(dv, frame.t), glm::dot(dv, frame.n));
        const Vector3 nLocal(glm::dot(pi.gn, frame.s), glm::dot(pi.gn, frame.t), glm::dot(pi.gn, frame.n));
        const std::array<Real, 3> rProj = {
            std::sqrt(dLocal.y * dLocal.y + dLocal.z * dLocal.z),   // probe along x: the offset plane is (y, z)
            std::sqrt(dLocal.z * dLocal.z + dLocal.x * dLocal.x),
            std::sqrt(dLocal.x * dLocal.x + dLocal.y * dLocal.y) };
        Real pdf = 0;
        for (int axis = 0; axis < 3; ++axis) {
            Real channels = 0;
            for (int c = 0; c < kSssChannels; ++c) channels += sss::profile(rProj[static_cast<size_t>(axis)], params.d[c]);
            pdf += channels / Real(kSssChannels) * std::abs(nLocal[axis]) * axisProb[static_cast<size_t>(axis)];
        }
        pdf /= Real(nFound);
        if (!(pdf > Real(0)) || !std::isfinite(pdf)) return false;

        const Real distance = glm::length(dv);
        weight = Spectrum(0.0);
        for (int c = 0; c < kSssChannels; ++c) weight[c] = params.albedo[c] * sss::profile(distance, params.d[c]) / pdf;

        // --- the exit vertex: outward normal, the exit lobe as its material --------------------
        exitVertex = pi;
        exitVertex.gn = pi.frontFace ? pi.gn : -pi.gn;   // scene.hit() flips the normals toward the probe ray
        exitVertex.n = pi.frontFace ? pi.n : -pi.n;
        exitVertex.frontFace = true;
        exitVertex.matPtr = material->subsurfaceExit();
        exitVertex.areaLight = nullptr;
        exitVertex.medium = po.medium;
        exitVertex.time = po.time;
        return exitVertex.matPtr != nullptr;
    }

} // namespace rayt
