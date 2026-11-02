/**
 * @file IO/VpdLoader.hpp
 * @brief Reader for MikuMikuDance pose files (.vpd, "Vocaloid Pose Data file").
 *
 * A VPD is a Shift-JIS text file with one block per bone:
 *
 *     Bone0{センター
 *       0.000000,0.000000,0.000000;      // trans x,y,z
 *       0.000000,0.000000,0.000000,1.000000;   // Quaternion x,y,z,w
 *     }
 *
 * - the bone name follows "{" on the same line (so it is matched against the model's bone names),
 * - translation is relative to the bone's rest position, in MMD units, left-handed,
 * - the quaternion is the bone's LOCAL rotation relative to its rest orientation (x, y, z, w).
 *
 * MMD saves the final values of every bone, including the IK-solved result in the FK bones of the
 * legs / arms, so a pose can be reproduced without running the IK solver (as long as the model has
 * roughly the proportions of the model the pose was made with).
 */

#pragma once

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "Core/Types.hpp"

#if defined(_WIN32)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace rayt::io {

    struct VpdBone {
        std::string name;                                  ///< UTF-8
        Vector3 translation = Vector3(0.0);
        Real qx = 0, qy = 0, qz = 0, qw = 1;               ///< local rotation quaternion
    };

    struct VpdPose {
        std::string parentModel;                           ///< "宮舞モカ.osm" etc. (informational)
        std::vector<VpdBone> bones;
    };

    /// Shift-JIS (code page 932) -> UTF-8. Text that is already valid UTF-8 is returned unchanged.
    inline std::string sjisToUtf8(const std::string& s) {
#if defined(_WIN32)
        if (s.empty()) return s;
        // Valid UTF-8 (strict) is kept as is.
        if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, s.data(), static_cast<int>(s.size()), nullptr, 0) > 0) return s;
        const int wlen = MultiByteToWideChar(932, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
        if (wlen <= 0) return s;
        std::wstring w(static_cast<size_t>(wlen), L'\0');
        MultiByteToWideChar(932, 0, s.data(), static_cast<int>(s.size()), w.data(), wlen);
        const int ulen = WideCharToMultiByte(CP_UTF8, 0, w.data(), wlen, nullptr, 0, nullptr, nullptr);
        std::string out(static_cast<size_t>(ulen), '\0');
        WideCharToMultiByte(CP_UTF8, 0, w.data(), wlen, out.data(), ulen, nullptr, nullptr);
        return out;
#else
        return s;
#endif
    }

    namespace detail {
        inline std::string trim(const std::string& s) {
            size_t a = 0, b = s.size();
            while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) ++a;
            while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) --b;
            return s.substr(a, b - a);
        }

        /// Numbers before the first ';' of a line ("x,y,z;  // comment"). Returns how many were read.
        inline int numbers(const std::string& line, double* out, int maxCount) {
            std::string body = line.substr(0, line.find(';'));
            int n = 0;
            const char* p = body.c_str();
            while (*p && n < maxCount) {
                char* end = nullptr;
                const double v = std::strtod(p, &end);
                if (end == p) { ++p; continue; }
                out[n++] = v;
                p = end;
                while (*p == ',' || *p == ' ' || *p == '\t') ++p;
            }
            return n;
        }
    } // namespace detail

    inline bool loadVpd(const std::filesystem::path& file, VpdPose& out, std::string* outError = nullptr) {
        auto fail = [&](const std::string& msg) { if (outError) *outError = msg; return false; };

        std::ifstream f(file, std::ios::binary);
        if (!f) return fail("cannot open the VPD file");
        const std::string raw((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        const std::string text = sjisToUtf8(raw);
        if (text.compare(0, 6, "Vocalo") != 0) return fail("not a VPD file (missing the 'Vocaloid Pose Data file' header)");

        out = VpdPose();
        std::vector<std::string> lines;
        {
            size_t start = 0;
            while (start <= text.size()) {
                size_t end = text.find('\n', start);
                if (end == std::string::npos) end = text.size();
                lines.push_back(text.substr(start, end - start));
                start = end + 1;
            }
        }

        for (size_t i = 0; i < lines.size(); ++i) {
            const std::string line = detail::trim(lines[i]);
            if (line.compare(0, 4, "Bone") == 0 && line.find('{') != std::string::npos) {
                VpdBone b;
                b.name = detail::trim(line.substr(line.find('{') + 1));
                double t[3] = {}, q[4] = { 0, 0, 0, 1 };
                // the next two non-empty lines: translation, quaternion
                size_t j = i + 1;
                auto nextLine = [&]() { while (j < lines.size() && detail::trim(lines[j]).empty()) ++j; return j < lines.size() ? detail::trim(lines[j++]) : std::string(); };
                if (detail::numbers(nextLine(), t, 3) != 3) return fail("bad translation line for bone " + b.name);
                if (detail::numbers(nextLine(), q, 4) != 4) return fail("bad quaternion line for bone " + b.name);
                b.translation = Vector3(t[0], t[1], t[2]);
                b.qx = q[0]; b.qy = q[1]; b.qz = q[2]; b.qw = q[3];
                out.bones.push_back(std::move(b));
                i = j - 1;
            }
            else if (line.find("osm;") != std::string::npos || line.find("pmx;") != std::string::npos || line.find("pmd;") != std::string::npos) {
                out.parentModel = detail::trim(line.substr(0, line.find(';')));
            }
        }
        if (out.bones.empty()) return fail("no Bone blocks found");
        return true;
    }

} // namespace rayt::io
