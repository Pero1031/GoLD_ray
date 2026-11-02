/**
 * @file IO/PmxMorph.hpp
 * @brief Applies PMX morphs (facial expressions and material switches) to a model.
 *
 * A morph is a named set of offsets scaled by a weight (usually 0..1):
 *  - vertex   : adds weight * offset to vertex positions (eye blink, mouth shapes, brows ...),
 *  - uv       : adds to texture coordinates (type 3 = the main UV set),
 *  - material : changes material parameters, with add or multiply semantics (a blush layer is switched
 *               on by raising its diffuse alpha from 0),
 *  - bone     : a translation / rotation of a bone (e.g. eye direction), returned as pose bones,
 *  - group    : a weighted list of other morphs.
 *
 * Vertex morphs must be applied to the REST pose, before skinning, so call applyMorphs() first and
 * applyPose() (IO/PmxPose.hpp) on its result.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <type_traits>
#include <string>
#include <utility>
#include <vector>

#include <glm/gtc/quaternion.hpp>

#include "IO/PmxLoader.hpp"
#include "IO/VpdLoader.hpp"

namespace rayt::io {

    namespace detail {
        inline void accumulateMorph(const PmxModel& model, int index, Real weight, int depth, std::vector<Real>& eff) {
            if (index < 0 || index >= static_cast<int>(model.morphs.size()) || depth > 8) return;
            eff[static_cast<size_t>(index)] += weight;
            const PmxMorph& m = model.morphs[static_cast<size_t>(index)];
            if (m.type == 0) for (const PmxGroupOffset& g : m.groupOffsets) accumulateMorph(model, g.morph, weight * g.rate, depth + 1, eff);
        }
    } // namespace detail

    /**
     * @param weights  morph index -> weight
     * @param boneOut  receives the bone morphs as pose bones (translation, rotation) if not null
     */
    inline PmxModel applyMorphs(const PmxModel& model, const std::map<int, Real>& weights, std::vector<VpdBone>* boneOut = nullptr) {
        PmxModel out = model;
        std::vector<Real> eff(model.morphs.size(), Real(0));
        for (const auto& [index, w] : weights) detail::accumulateMorph(model, index, w, 0, eff);

        for (size_t mi = 0; mi < model.morphs.size(); ++mi) {
            const Real w = eff[mi];
            if (w == Real(0)) continue;
            const PmxMorph& m = model.morphs[mi];
            switch (m.type) {
            case 1:
                for (const auto& o : m.vertexOffsets)
                    if (o.vertex >= 0 && o.vertex < static_cast<int>(out.vertices.size()))
                        out.vertices[static_cast<size_t>(o.vertex)].position += Point3(o.offset * w);
                break;
            case 3:
                for (const auto& o : m.uvOffsets)
                    if (o.vertex >= 0 && o.vertex < static_cast<int>(out.vertices.size()))
                        out.vertices[static_cast<size_t>(o.vertex)].uv += Point2(o.offset.x * w, o.offset.y * w);
                break;
            case 8:
                for (const auto& o : m.materialOffsets) {
                    for (size_t k = 0; k < out.materials.size(); ++k) {
                        if (o.material >= 0 && o.material != static_cast<int>(k)) continue;
                        PmxMaterial& mat = out.materials[k];
                        const Real before = mat.diffuse.w;
                        auto blend = [&](auto& value, const auto& delta) {
                            if (o.op == 1) value += delta * w;                         // add
                            else value = value * (Real(1) + (delta - std::decay_t<decltype(value)>(Real(1))) * w);   // multiply: lerp(1, delta, w)
                        };
                        blend(mat.diffuse, o.diffuse);
                        blend(mat.specular, o.specular);
                        blend(mat.ambient, o.ambient);
                        blend(mat.edgeColor, o.edgeColor);
                        if (o.op == 1) { mat.shininess += o.shininess * w; mat.edgeSize += o.edgeSize * w; }
                        else { mat.shininess *= Real(1) + (o.shininess - Real(1)) * w; mat.edgeSize *= Real(1) + (o.edgeSize - Real(1)) * w; }
                        if (std::abs(mat.diffuse.w - before) > Real(1e-6)) mat.alphaFromMorph = true;
                    }
                }
                break;
            case 2:
                if (boneOut) {
                    for (const auto& o : m.boneOffsets) {
                        if (o.bone < 0 || o.bone >= static_cast<int>(model.bones.size())) continue;
                        VpdBone vb;
                        vb.name = model.bones[static_cast<size_t>(o.bone)].name;
                        vb.translation = o.translation * w;
                        const glm::dquat q(o.rotation.w, o.rotation.x, o.rotation.y, o.rotation.z);
                        const glm::dquat part = glm::slerp(glm::dquat(1, 0, 0, 0), glm::normalize(q), w);
                        vb.qx = part.x; vb.qy = part.y; vb.qz = part.z; vb.qw = part.w;
                        boneOut->push_back(vb);
                    }
                }
                break;
            default: break;   // group handled above; flip / impulse / extra UV sets are ignored
            }
        }
        return out;
    }

    /// Adds `extra` bone transforms (from bone morphs) to a pose: same bone = translations add, rotations multiply.
    inline void mergePoseBones(VpdPose& pose, const std::vector<VpdBone>& extra) {
        for (const VpdBone& e : extra) {
            auto it = std::find_if(pose.bones.begin(), pose.bones.end(), [&](const VpdBone& b) { return b.name == e.name; });
            if (it == pose.bones.end()) { pose.bones.push_back(e); continue; }
            it->translation += e.translation;
            const glm::dquat q = glm::dquat(e.qw, e.qx, e.qy, e.qz) * glm::dquat(it->qw, it->qx, it->qy, it->qz);
            it->qx = q.x; it->qy = q.y; it->qz = q.z; it->qw = q.w;
        }
    }

    /**
     * @brief Reads an expression file: UTF-8 text, one "morph name = weight" per line ('#' starts a comment).
     *        A name that is a plain integer is a morph index.
     */
    inline bool loadExpression(const std::filesystem::path& file, std::vector<std::pair<std::string, Real>>& out, std::string* outError = nullptr) {
        std::ifstream f(file);
        if (!f) { if (outError) *outError = "cannot open the expression file"; return false; }
        out.clear();
        std::string line;
        while (std::getline(f, line)) {
            const size_t hash = line.find('#');
            if (hash != std::string::npos) line.resize(hash);
            const size_t eq = line.find('=');
            if (eq == std::string::npos) continue;
            const std::string name = detail::trim(line.substr(0, eq));
            const std::string val = detail::trim(line.substr(eq + 1));
            if (name.empty() || val.empty()) continue;
            out.emplace_back(name, Real(std::strtod(val.c_str(), nullptr)));
        }
        return !out.empty();
    }

} // namespace rayt::io
