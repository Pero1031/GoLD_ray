/**
 * @file IO/PmxPose.hpp
 * @brief Poses a PMX model: bone transforms (FK + grant / append bones) and linear blend skinning.
 *
 * Pipeline, in MMD's own (left-handed) model space; the caller converts to the renderer's space
 * afterwards, so no quaternion mirroring is needed here:
 *
 *  1. Local animation per bone from the pose: translation t_i and rotation q_i (identity by default).
 *  2. Grant ("append") bones, evaluated in deform order (layer, then index):
 *        q_eff = q_i * slerp(identity, q_eff[grantParent], weight)     (rotation grant)
 *        t_eff = t_i + t_eff[grantParent] * weight                     (translation grant)
 *     Many models skin their legs / twist joints to such bones (leg "D" bones, arm twist bones),
 *     so this step is essential, not optional.
 *  3. Global transforms from the root down:
 *        G_i = G_parent * Translate(rest_i - rest_parent + t_eff_i) * Rotate(q_eff_i)
 *  4. Skin matrix S_i = G_i * Translate(-rest_i), and  v' = sum_k w_k S_k v  (normals with the
 *     rotation part, renormalized).
 *
 * NOT done: IK solving (the FK bones of an MMD pose file already hold the IK-solved result), fixed
 * axis / local axis bones, SDEF (treated as BDEF2), morphs, physics.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "IO/PmxLoader.hpp"
#include "IO/VpdLoader.hpp"

namespace rayt::io {

    struct PoseReport {
        int applied = 0;                                   ///< VPD bones found in the model
        std::vector<std::string> missing;                  ///< VPD bones the model does not have
        std::vector<glm::dquat> effectiveRotation;         ///< per model bone (for tests)
        std::vector<glm::dmat4> skin;                      ///< per model bone
    };

    namespace detail {
        /// q^w (rotation by w times the angle about the same axis); w may be negative or fractional.
        inline glm::dquat quatPow(const glm::dquat& q, Real w) {
            const glm::dquat n = glm::normalize(q);
            const Real angle = Real(2) * std::acos(std::clamp(n.w, Real(-1), Real(1)));
            const Real s = std::sqrt(std::max(Real(0), Real(1) - n.w * n.w));
            if (s < Real(1e-9) || std::abs(angle) < Real(1e-12)) return glm::dquat(1, 0, 0, 0);
            const glm::dvec3 axis(n.x / s, n.y / s, n.z / s);
            return glm::angleAxis(angle * w, axis);
        }
    } // namespace detail

    /**
     * @brief Returns the model with its vertices (positions, normals) moved to the pose.
     * @param model  rest-pose model
     * @param pose   VPD pose; bones are matched by name, unmatched model bones stay at rest
     */
    inline PmxModel applyPose(const PmxModel& model, const VpdPose& pose, PoseReport* report = nullptr) {
        const size_t nb = model.bones.size();
        std::map<std::string, int> byName;
        for (size_t i = 0; i < nb; ++i) byName.emplace(model.bones[i].name, static_cast<int>(i));

        std::vector<Vector3> animT(nb, Vector3(0.0));
        std::vector<glm::dquat> animQ(nb, glm::dquat(1, 0, 0, 0));
        PoseReport local;
        PoseReport& rep = report ? *report : local;
        rep = PoseReport();
        for (const VpdBone& vb : pose.bones) {
            const auto it = byName.find(vb.name);
            if (it == byName.end()) { rep.missing.push_back(vb.name); continue; }
            animT[static_cast<size_t>(it->second)] = vb.translation;
            animQ[static_cast<size_t>(it->second)] = glm::normalize(glm::dquat(vb.qw, vb.qx, vb.qy, vb.qz));   // glm order: w, x, y, z
            ++rep.applied;
        }

        // Deform order: by layer, then by index.
        std::vector<int> order(nb);
        for (size_t i = 0; i < nb; ++i) order[i] = static_cast<int>(i);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return model.bones[static_cast<size_t>(a)].layer < model.bones[static_cast<size_t>(b)].layer; });

        std::vector<glm::dquat> qEff(nb, glm::dquat(1, 0, 0, 0));
        std::vector<Vector3> tEff(nb, Vector3(0.0));
        for (int i : order) {
            const PmxBone& b = model.bones[static_cast<size_t>(i)];
            glm::dquat q = animQ[static_cast<size_t>(i)];
            Vector3 t = animT[static_cast<size_t>(i)];
            if (b.grantParent >= 0 && b.grantParent < static_cast<int>(nb)) {
                if (b.inheritsRotation()) q = q * detail::quatPow(qEff[static_cast<size_t>(b.grantParent)], b.grantWeight);
                if (b.inheritsTranslation()) t += tEff[static_cast<size_t>(b.grantParent)] * b.grantWeight;
            }
            qEff[static_cast<size_t>(i)] = glm::normalize(q);
            tEff[static_cast<size_t>(i)] = t;
        }

        // Global transforms (recursive, memoized: a parent may have a larger index than its child).
        std::vector<glm::dmat4> G(nb, glm::dmat4(1.0));
        std::vector<char> done(nb, 0);
        std::vector<int> stack;
        for (size_t start = 0; start < nb; ++start) {
            stack.clear();
            int cur = static_cast<int>(start);
            while (cur >= 0 && !done[static_cast<size_t>(cur)]) { stack.push_back(cur); cur = model.bones[static_cast<size_t>(cur)].parent; if (stack.size() > nb) break; }
            while (!stack.empty()) {
                const int i = stack.back(); stack.pop_back();
                const PmxBone& b = model.bones[static_cast<size_t>(i)];
                const bool hasParent = b.parent >= 0 && b.parent < static_cast<int>(nb);
                const Vector3 rel = hasParent ? Vector3(b.position) - Vector3(model.bones[static_cast<size_t>(b.parent)].position) : Vector3(b.position);
                const glm::dmat4 local = glm::translate(glm::dmat4(1.0), glm::dvec3(rel + tEff[static_cast<size_t>(i)])) * glm::mat4_cast(qEff[static_cast<size_t>(i)]);
                G[static_cast<size_t>(i)] = hasParent ? G[static_cast<size_t>(b.parent)] * local : local;
                done[static_cast<size_t>(i)] = 1;
            }
        }

        std::vector<glm::dmat4> S(nb);
        for (size_t i = 0; i < nb; ++i) S[i] = G[i] * glm::translate(glm::dmat4(1.0), -glm::dvec3(model.bones[i].position));
        rep.effectiveRotation = qEff;
        rep.skin = S;

        PmxModel out = model;
        for (PmxVertex& v : out.vertices) {
            glm::dvec4 p(0.0);
            glm::dvec3 n(0.0);
            Real wsum = 0;
            for (int k = 0; k < 4; ++k) {
                if (v.bones[k] < 0 || v.weights[k] == 0) continue;
                const glm::dmat4& m = S[static_cast<size_t>(v.bones[k])];
                p += v.weights[k] * (m * glm::dvec4(glm::dvec3(v.position), 1.0));
                n += v.weights[k] * (glm::dmat3(m) * glm::dvec3(v.normal));
                wsum += v.weights[k];
            }
            if (wsum <= Real(0)) continue;                 // not skinned: stays where it is
            v.position = Point3(p.x, p.y, p.z);
            const Real len = glm::length(n);
            if (len > Real(1e-12)) v.normal = Vector3(n / len);
        }
        return out;
    }

} // namespace rayt::io
