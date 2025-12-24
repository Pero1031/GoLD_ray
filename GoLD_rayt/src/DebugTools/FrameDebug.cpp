#include "pch.h"                  
#include "Core/Types.hpp"
#include "Core/Core.hpp"
#include "Core/Constants.hpp"
#include "Geometry/Frame.hpp"           // 置いてる場所に合わせて include
#include "DebugTools/FrameDebug.hpp"
#include <iostream>
#include <limits>
#include <algorithm>

namespace rayt::debug {

    void TestFrameRoundTrip() {
        using rayt::Vector3;
        using rayt::Real;

        std::cout << "\n[Debug] Frame round-trip test: worldToLocal(localToWorld(a))\n";

        Vector3 n_world = glm::normalize(Vector3(0.3, 0.9, 0.1));
        Vector3 tangent_world = glm::normalize(Vector3(1.0, 0.2, 0.0));

        rayt::frame::Frame F;
        F.buildFromNormalAndTangent(n_world, tangent_world);

        auto test_one = [&](const Vector3& a_local, const char* name) {
            Vector3 w = F.localToWorld(a_local);
            Vector3 b = F.worldToLocal(w);
            Vector3 diff = b - a_local;
            Real err = glm::length(diff);

            std::cout << "  " << name
                << "  a=" << a_local.x << "," << a_local.y << "," << a_local.z
                << "  b=" << b.x << "," << b.y << "," << b.z
                << "  |b-a|=" << err << "\n";
            return err;
            };

        Real max_err = 0.0;
        max_err = std::max(max_err, test_one(glm::normalize(Vector3(1, 0, 0)), "ex"));
        max_err = std::max(max_err, test_one(glm::normalize(Vector3(0, 1, 0)), "ey"));
        max_err = std::max(max_err, test_one(glm::normalize(Vector3(0, 0, 1)), "ez"));
        max_err = std::max(max_err, test_one(glm::normalize(Vector3(1, 1, 1)), "e111"));
        max_err = std::max(max_err, test_one(glm::normalize(Vector3(0.3, 0.4, 0.8660254)), "misc"));

        Real ss = glm::dot(F.s, F.s);
        Real tt = glm::dot(F.t, F.t);
        Real nn = glm::dot(F.n, F.n);
        Real st = glm::dot(F.s, F.t);
        Real sn = glm::dot(F.s, F.n);
        Real tn = glm::dot(F.t, F.n);

        std::cout << "  [ONB] |s|^2=" << ss << " |t|^2=" << tt << " |n|^2=" << nn << "\n";
        std::cout << "  [ONB] s·t=" << st << " s·n=" << sn << " t·n=" << tn << "\n";
        std::cout << "  [ONB] max round-trip error = " << max_err << "\n";

        if (max_err > Real(1e-5) ||
            std::abs(ss - Real(1.0)) > Real(1e-4) ||
            std::abs(tt - Real(1.0)) > Real(1e-4) ||
            std::abs(nn - Real(1.0)) > Real(1e-4)) {
            std::cout << "  [WARN] Frame may not be orthonormal enough.\n";
        }
        else {
            std::cout << "  [OK] Frame looks orthonormal.\n";
        }
    }

} // namespace rayt::debug
