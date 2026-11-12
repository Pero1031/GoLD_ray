// PMX skeleton and pose checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "IO/PmxPose.hpp"
#include "Core/Constants.hpp"

namespace rayt::cli {

bool RunPmxChecks(const io::PmxModel& model) {
    bool allPassed = true;
    auto check = [&](const std::string& what, bool ok, const std::string& detail) {
        std::cout << "[Pose]   " << what << ": " << detail << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };

    const size_t nb = model.bones.size();
    int ik = 0, grantRot = 0, grantTrans = 0, afterPhysics = 0, roots = 0;
    for (const auto& b : model.bones) {
        ik += b.isIk(); grantRot += b.inheritsRotation(); grantTrans += b.inheritsTranslation(); afterPhysics += b.afterPhysics();
        roots += (b.parent < 0);
    }
    int unskinned = 0;
    for (const auto& v : model.vertices) if (v.bones[0] < 0) ++unskinned;
    std::cout << "[Pose] skeleton: " << nb << " bones (" << roots << " roots, " << ik << " IK, " << grantRot << " rotation-grant, " << grantTrans
        << " translation-grant, " << afterPhysics << " deformed after physics), " << unskinned << " vertices without bones\n";
    {
        std::map<int, int> usage;   // bones that actually move vertices
        for (const auto& v : model.vertices) for (int k = 0; k < 4; ++k) if (v.bones[k] >= 0 && v.weights[k] > 0) ++usage[v.bones[k]];
        int usedGrant = 0, usedIk = 0;
        for (const auto& [bone, n] : usage) { usedGrant += model.bones[static_cast<size_t>(bone)].inheritsRotation(); usedIk += model.bones[static_cast<size_t>(bone)].isIk(); }
        std::cout << "[Pose] skinning uses " << usage.size() << " bones, " << usedGrant << " of them grant-rotation bones, " << usedIk << " IK bones\n";
    }

    // 1. identity pose must reproduce the rest pose
    {
        io::VpdPose none;
        io::VpdBone dummy; dummy.name = "__none__";
        none.bones.push_back(dummy);
        const io::PmxModel same = io::applyPose(model, none);
        Real worst = 0;
        for (size_t i = 0; i < model.vertices.size(); ++i) {
            worst = std::max(worst, glm::length(Vector3(same.vertices[i].position) - Vector3(model.vertices[i].position)));
            worst = std::max(worst, glm::length(Vector3(same.vertices[i].normal) - Vector3(model.vertices[i].normal)));
        }
        check("rest pose reproduced (max vertex / normal difference; normals are renormalized)", worst < 1e-5, std::to_string(worst));
    }

    // 2. rotate the root-most "center" bone by 90 degrees about +Y: everything under it must turn rigidly about its rest position
    {
        int center = -1;
        for (size_t i = 0; i < nb; ++i) if (model.bones[i].name == "センター") { center = static_cast<int>(i); break; }
        if (center < 0) { check("center bone test", false, "no bone named センター"); }
        else {
            io::VpdPose p;
            io::VpdBone b; b.name = "センター"; b.qy = std::sin(constants::PI / 4); b.qw = std::cos(constants::PI / 4);
            p.bones.push_back(b);
            const io::PmxModel moved = io::applyPose(model, p);
            const Vector3 c(model.bones[static_cast<size_t>(center)].position);
            const glm::dmat3 R = glm::dmat3(glm::mat4_cast(glm::angleAxis(constants::PI / 2, glm::dvec3(0, 1, 0))));
            // descendants of the center bone
            std::vector<char> under(nb, 0);
            for (size_t i = 0; i < nb; ++i) { int cur = static_cast<int>(i); for (size_t guard = 0; guard <= nb && cur >= 0; ++guard) { if (cur == center) { under[i] = 1; break; } cur = model.bones[static_cast<size_t>(cur)].parent; } }
            int total = 0, matching = 0;
            for (size_t i = 0; i < model.vertices.size(); ++i) {
                const auto& v = model.vertices[i];
                bool all = v.bones[0] >= 0;
                for (int k = 0; k < 4 && all; ++k) if (v.bones[k] >= 0 && v.weights[k] > 0 && !under[static_cast<size_t>(v.bones[k])]) all = false;
                if (!all) continue;
                ++total;
                const Vector3 expected = c + R * (Vector3(v.position) - c);
                if (glm::length(Vector3(moved.vertices[i].position) - expected) < 1e-6) ++matching;
            }
            check("rotating センター by 90 deg about +Y turns " + std::to_string(total) + " vertices rigidly about its rest position", total > 0 && matching == total,
                std::to_string(matching) + " of " + std::to_string(total) + " match");
        }
    }

    // 3. grant bones: with an identity own rotation, a bone inheriting rotation with weight w must end up with q^w of its grant parent
    {
        int tested = 0, passed = 0;
        for (size_t i = 0; i < nb && tested < 6; ++i) {
            const auto& b = model.bones[i];
            if (!b.inheritsRotation() || b.grantParent < 0) continue;
            const auto& gp = model.bones[static_cast<size_t>(b.grantParent)];
            if (gp.inheritsRotation()) continue;   // keep the test simple: a plain FK source
            io::VpdPose p;
            io::VpdBone vb; vb.name = gp.name; vb.qx = 0.3; vb.qy = 0.2; vb.qz = 0.1; vb.qw = std::sqrt(1 - 0.09 - 0.04 - 0.01);
            p.bones.push_back(vb);
            io::PoseReport rep;
            io::applyPose(model, p, &rep);
            const glm::dquat src = rep.effectiveRotation[static_cast<size_t>(b.grantParent)];
            const glm::dquat expected = io::detail::quatPow(src, b.grantWeight);
            const glm::dquat got = rep.effectiveRotation[i];
            ++tested;
            if (std::abs(std::abs(glm::dot(expected, got)) - 1.0) < 1e-9) ++passed;
            std::cout << "[Pose]     grant " << model.bones[i].name << " <- " << gp.name << " x" << b.grantWeight << ": "
                << (std::abs(std::abs(glm::dot(expected, got)) - 1.0) < 1e-9 ? "ok" : "MISMATCH") << "\n";
        }
        check("grant (append) rotation of " + std::to_string(tested) + " sample bones", tested == 0 || passed == tested, std::to_string(passed) + " of " + std::to_string(tested));
    }
    std::cout << "[Pose] checks: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

} // namespace rayt::cli
