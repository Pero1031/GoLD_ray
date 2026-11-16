// hpbrdf_pca.hpp - header-only decoder for PCA-compressed BRDF tables (.hpbpca, format v2). No dependencies, C++17.
// Home: brdfkit/cpp/hpbrdf_pca.hpp (the exporter is brdfkit/export/hpbpca.py; tests/ check both agree).
//
// File layout (little endian), written by  python -m brdfkit export OUT.npz OUT.hpbpca
//   char     magic[8] = "HPBPCA02"
//   uint32   n_phi, n_thd, n_thh, n_ch, k1, k2, coef_type (0 = float16, 1 = float32), flags
//            flags bit 0: polarized (Mueller part present);  bit 1: channels are RGB (else wavelengths in nm)
//   float    eps
//   float    phi_d[n_phi], theta_d[n_thd], theta_h[n_thh]        // radians, ascending
//   float    channels[n_ch]                                      // nm, or 0 1 2 for RGB
//   float    mean_logm00[n_ch],  basis_logm00[n_ch][k1]
//   float    mean_ratio[n_ch*15], basis_ratio[n_ch*15][k2]       // only when polarized; ratio index = ch*15 + (e-1), e = row*4+col
//   uint8    valid[n_phi][n_thd][n_thh]                          // 1 where the source table had data
//   coef_t   coef[n_phi][n_thd][n_thh][k1+k2]                    // latent vector per angular bin
//
// Decode:  m00(ch)  = exp( mean_logm00(ch) + basis_logm00(ch) . c[0:k1] )
//          m_e(ch)  = m00(ch) . ( mean_ratio(ch,e) + basis_ratio(ch,e) . c[k1:k1+k2] )
//
// Conventions (identical to Mitsuba 3 `measured_polarized`, see brdfkit/geometry.py):
//   * local shading frame, normal = +z, unit vectors
//   * to_light : surface -> light   (light propagates along -to_light)
//   * to_view  : surface -> viewer  (light leaves along +to_view)
//   * returned Mueller matrix / m00 is f_r without the cosine factor
//   * M is expressed in Stokes frames whose x reference vectors (x_in for the incident beam along -to_light,
//     x_out for the outgoing beam along +to_view) lie in the plane of reflection. Use rotate_mueller_basis()
//     to move it into your renderer's frames. Sign conventions of rotator()/unit_angle follow Mitsuba 3.
//   * the channel coordinate passed to m00()/mueller() is a wavelength in nm, or 0 / 1 / 2 for RGB tables.
#pragma once
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace hpbrdf {

struct Vec3 {
    float x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
};
inline Vec3  operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3  operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3  operator-(Vec3 a) { return {-a.x, -a.y, -a.z}; }
inline Vec3  operator*(Vec3 a, float s) { return {a.x * s, a.y * s, a.z * s}; }
inline float dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3  cross(Vec3 a, Vec3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline float length(Vec3 a) { return std::sqrt(dot(a, a)); }
inline Vec3  normalize(Vec3 a) { float l = length(a); return l > 0 ? a * (1.f / l) : a; }
inline Vec3  rotate_z(Vec3 v, float ang) { float c = std::cos(ang), s = std::sin(ang); return {v.x * c - v.y * s, v.x * s + v.y * c, v.z}; }
inline float safe_acos(float x) { return std::acos(x < -1.f ? -1.f : (x > 1.f ? 1.f : x)); }

struct Mueller { float m[4][4]; };   // m[row][col]

// ---------------------------------------------------------------------------- geometry
struct Geometry {
    float phi_d, theta_d, theta_h;   // Rusinkiewicz parameters as stored in the table
    Vec3  x_in, x_out;               // Stokes reference vectors of the table's Mueller matrix
};

// Mitsuba directions_to_rusinkiewicz(i, o): i = to_light (std frame), o = to_view (std frame)
inline void rusinkiewicz(Vec3 i, Vec3 o, float &phi_d, float &theta_h, float &theta_d) {
    Vec3 h = normalize(i + o), n(0, 0, 1);
    Vec3 b = cross(n, h);
    b = length(b) < 1e-7f ? Vec3(0, 1, 0) : normalize(b);
    Vec3 t = normalize(cross(b, h));
    theta_d = safe_acos(dot(h, i));
    theta_h = safe_acos(h.z);
    Vec3 i_prj = normalize(i - h * dot(i, h));
    float cp = dot(t, i_prj), sp = dot(b, i_prj);
    cp = cp < -1 ? -1 : (cp > 1 ? 1 : cp); sp = sp < -1 ? -1 : (sp > 1 ? 1 : sp);
    phi_d = std::atan2(sp, cp);
}

inline Geometry geometry(Vec3 to_light, Vec3 to_view) {
    float phi_std = std::atan2(to_view.y, to_view.x);
    if (phi_std < 0) phi_std += 2.f * 3.14159265358979f;
    Vec3 wo = rotate_z(to_light, -phi_std), wi = rotate_z(to_view, -phi_std);
    Geometry g;
    rusinkiewicz(wo, wi, g.phi_d, g.theta_h, g.theta_d);
    Vec3 d = wo - wi;
    if (length(d) < 1e-7f) d = Vec3(0, 1, 0);           // retro-reflection: plane of reflection undefined
    Vec3 zo = -wo, to = normalize(cross(d, zo)), yo = normalize(cross(to, zo)), xo = cross(yo, zo);
    Vec3 zi = wi,  ti = normalize(cross(-d, zi)), yi = normalize(cross(ti, zi)), xi = cross(yi, zi);
    g.x_in  = rotate_z(xo, phi_std);
    g.x_out = rotate_z(xi, phi_std);
    return g;
}

// ---------------------------------------------------------------------------- Stokes basis helpers (Mitsuba 3)
inline Vec3 default_basis(Vec3 w) {               // mitsuba::mueller::stokes_basis(w) = coordinate_system(w).first
    float sign = w.z >= 0 ? 1.f : -1.f, a = -1.f / (sign + w.z), b = w.x * w.y * a;
    return {w.x * w.x * a * sign + 1.f, b * sign, -w.x * sign};
}
inline Mueller rotator(float theta) {
    float s = std::sin(2 * theta), c = std::cos(2 * theta);
    Mueller R{}; R.m[0][0] = 1; R.m[3][3] = 1; R.m[1][1] = c; R.m[1][2] = s; R.m[2][1] = -s; R.m[2][2] = c;
    return R;
}
inline Mueller rotate_stokes_basis(Vec3 forward, Vec3 cur, Vec3 tgt) {
    cur = normalize(cur); tgt = normalize(tgt);
    float h = 0.5f * length(tgt - cur); h = h > 1 ? 1 : h;
    float theta = 2.f * std::asin(h);
    if (dot(forward, cross(cur, tgt)) < 0) theta = -theta;
    return rotator(theta);
}
inline Mueller mul(const Mueller &A, const Mueller &B) {
    Mueller C{};
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) { float s = 0; for (int k = 0; k < 4; ++k) s += A.m[i][k] * B.m[k][j]; C.m[i][j] = s; }
    return C;
}
inline Mueller transpose(const Mueller &A) { Mueller T; for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) T.m[i][j] = A.m[j][i]; return T; }
// M expressed in (in_cur, out_cur) reference frames  ->  M in (in_tgt, out_tgt)
inline Mueller rotate_mueller_basis(const Mueller &M, Vec3 in_forward, Vec3 in_cur, Vec3 in_tgt,
                                    Vec3 out_forward, Vec3 out_cur, Vec3 out_tgt) {
    Mueller R_in = rotate_stokes_basis(in_forward, in_cur, in_tgt), R_out = rotate_stokes_basis(out_forward, out_cur, out_tgt);
    return mul(mul(R_out, M), transpose(R_in));
}

// ---------------------------------------------------------------------------- float16
inline float half_to_float(uint16_t h) {
    uint32_t s = (h >> 15) & 1u, e = (h >> 10) & 0x1fu, m = h & 0x3ffu, f;
    if (e == 0) { float v = (float) m * (1.f / 1024.f) * 6.103515625e-05f; return s ? -v : v; }
    if (e == 31) f = (s << 31) | 0x7f800000u | (m << 13);
    else         f = (s << 31) | ((e + 112u) << 23) | (m << 13);
    float out; std::memcpy(&out, &f, 4); return out;
}

// ---------------------------------------------------------------------------- the table
class PCATable {
public:
    uint32_t n_phi = 0, n_thd = 0, n_thh = 0, n_ch = 0, k1 = 0, k2 = 0;
    bool polarized = false, rgb = false;
    float eps = 1e-4f;
    std::vector<float> phi_d, theta_d, theta_h, channels;
    std::vector<float> mean_logm00, basis_logm00;   // [ch], [ch][k1]
    std::vector<float> mean_ratio, basis_ratio;     // [ch*15+e], [ch*15+e][k2]   (empty when !polarized)
    std::vector<uint8_t> valid;                     // [phi][thd][thh]
    std::vector<float> coef;                        // [phi][thd][thh][k]  (converted to float32 on load)

    uint32_t k() const { return k1 + k2; }

    bool load(const std::string &path, std::string *err = nullptr) {
#ifdef _MSC_VER
        FILE *f = nullptr;
        if (fopen_s(&f, path.c_str(), "rb") != 0) f = nullptr;
#else
        FILE *f = std::fopen(path.c_str(), "rb");
#endif
        if (!f) { if (err) *err = "cannot open " + path; return false; }
        char magic[8]; uint32_t hdr[8];
        if (std::fread(magic, 1, 8, f) != 8 || std::memcmp(magic, "HPBPCA02", 8) != 0) {
            if (err) *err = std::memcmp(magic, "HPBPCA", 6) == 0 ? "old .hpbpca format: re-export with brdfkit" : "bad magic";
            std::fclose(f); return false;
        }
        bool ok = std::fread(hdr, 4, 8, f) == 8 && std::fread(&eps, 4, 1, f) == 1;
        n_phi = hdr[0]; n_thd = hdr[1]; n_thh = hdr[2]; n_ch = hdr[3]; k1 = hdr[4]; k2 = hdr[5];
        uint32_t coef_type = hdr[6], flags = hdr[7];
        polarized = (flags & 1u) != 0; rgb = (flags & 2u) != 0;
        auto rd = [&](std::vector<float> &v, size_t n) { v.resize(n); return n == 0 || std::fread(v.data(), 4, n, f) == n; };
        ok = ok && rd(phi_d, n_phi) && rd(theta_d, n_thd) && rd(theta_h, n_thh) && rd(channels, n_ch) &&
             rd(mean_logm00, n_ch) && rd(basis_logm00, (size_t) n_ch * k1);
        if (polarized) ok = ok && rd(mean_ratio, (size_t) n_ch * 15) && rd(basis_ratio, (size_t) n_ch * 15 * k2);
        size_t nb = (size_t) n_phi * n_thd * n_thh;
        valid.resize(nb); ok = ok && std::fread(valid.data(), 1, nb, f) == nb;
        coef.resize(nb * k());
        if (coef_type == 1) ok = ok && std::fread(coef.data(), 4, coef.size(), f) == coef.size();
        else {
            std::vector<uint16_t> tmp(coef.size());
            ok = ok && std::fread(tmp.data(), 2, tmp.size(), f) == tmp.size();
            for (size_t i = 0; i < tmp.size(); ++i) coef[i] = half_to_float(tmp[i]);
        }
        std::fclose(f);
        if (!ok && err) *err = "truncated file";
        return ok;
    }

    // --- axis lookup: index j and fraction t such that x ~ axis[j]*(1-t) + axis[j+1]*t
    static void locate(const std::vector<float> &ax, float x, uint32_t &j, float &t) {
        uint32_t n = (uint32_t) ax.size();
        if (n < 2) { j = 0; t = 0; return; }
        if (x <= ax[0]) { j = 0; t = 0; return; }
        if (x >= ax[n - 1]) { j = n - 2; t = 1; return; }
        uint32_t lo = 0, hi = n - 1;
        while (hi - lo > 1) { uint32_t mid = (lo + hi) >> 1; if (ax[mid] <= x) lo = mid; else hi = mid; }
        j = lo; float w = ax[j + 1] - ax[j]; t = w > 0 ? (x - ax[j]) / w : 0.f;
    }

    const float *bin(uint32_t ip, uint32_t id, uint32_t ih) const { return &coef[(((size_t) ip * n_thd + id) * n_thh + ih) * k()]; }

    // --- trilinear interpolation of the latent vector; `c` must hold k() floats
    void coeffs(float phi_d_, float theta_d_, float theta_h_, float *c) const {
        uint32_t ip, id, ih; float tp, td, th;
        locate(phi_d, phi_d_, ip, tp); locate(theta_d, theta_d_, id, td); locate(theta_h, theta_h_, ih, th);
        for (uint32_t q = 0; q < k(); ++q) c[q] = 0;
        for (int dp = 0; dp < 2; ++dp) for (int dd = 0; dd < 2; ++dd) for (int dh = 0; dh < 2; ++dh) {
            float w = (dp ? tp : 1 - tp) * (dd ? td : 1 - td) * (dh ? th : 1 - th);
            if (w == 0) continue;
            const float *b = bin(ip + dp, id + dd, ih + dh);
            for (uint32_t q = 0; q < k(); ++q) c[q] += w * b[q];
        }
    }

    // --- channel lookup (wavelength in nm, or 0 / 1 / 2 for RGB)
    void locate_channel(float ch, uint32_t &j, float &t) const {
        if (n_ch == 1) { j = 0; t = 0; return; }
        locate(channels, ch, j, t);
    }
    const float *mean1_at(uint32_t j) const { return &mean_logm00[j]; }

    // m00 (= unpolarized f_r) at channel coordinate ch
    float m00(const float *c, float ch) const {
        uint32_t j; float t; locate_channel(ch, j, t);
        uint32_t j1 = (n_ch == 1) ? j : j + 1;
        float v = (1 - t) * mean_logm00[j] + t * mean_logm00[j1];
        const float *b0 = &basis_logm00[(size_t) j * k1], *b1 = &basis_logm00[(size_t) j1 * k1];
        for (uint32_t q = 0; q < k1; ++q) v += ((1 - t) * b0[q] + t * b1[q]) * c[q];
        return std::exp(v);
    }
    // m00 at table channel index w (no interpolation)
    float m00_band(const float *c, uint32_t w) const {
        float v = mean_logm00[w]; const float *b = &basis_logm00[(size_t) w * k1];
        for (uint32_t q = 0; q < k1; ++q) v += b[q] * c[q];
        return std::exp(v);
    }
    // full Mueller matrix (table frame) at channel coordinate ch; for an unpolarized table only m[0][0] is set
    Mueller mueller(const float *c, float ch) const {
        Mueller M{}; M.m[0][0] = m00(c, ch);
        if (!polarized) return M;
        uint32_t j; float t; locate_channel(ch, j, t);
        uint32_t j1 = (n_ch == 1) ? j : j + 1;
        const float *c2 = c + k1; float i00 = M.m[0][0];
        for (uint32_t e = 1; e < 16; ++e) {
            size_t r0 = (size_t) j * 15 + (e - 1), r1 = (size_t) j1 * 15 + (e - 1);
            float v = (1 - t) * mean_ratio[r0] + t * mean_ratio[r1];
            const float *b0 = &basis_ratio[r0 * k2], *b1 = &basis_ratio[r1 * k2];
            for (uint32_t q = 0; q < k2; ++q) v += ((1 - t) * b0[q] + t * b1[q]) * c2[q];
            M.m[e >> 2][e & 3] = v * i00;
        }
        return M;
    }

    // --- convenience: full evaluation from directions
    struct Result { Geometry g; float m00; Mueller M; };
    Result eval(Vec3 to_light, Vec3 to_view, float ch, bool want_polarized = true, float *scratch = nullptr) const {
        Result r; r.g = geometry(to_light, to_view);
        std::vector<float> tmp; float *c = scratch;
        if (!c) { tmp.resize(k()); c = tmp.data(); }
        coeffs(r.g.phi_d, r.g.theta_d, r.g.theta_h, c);
        if (want_polarized && polarized) { r.M = mueller(c, ch); r.m00 = r.M.m[0][0]; }
        else { r.m00 = m00(c, ch); r.M = Mueller{}; r.M.m[0][0] = r.m00; }
        return r;
    }
    // Mueller matrix rotated into Mitsuba's implicit Stokes frames (default_basis of the propagation directions).
    Mueller eval_default_frames(Vec3 to_light, Vec3 to_view, float ch) const {
        Result r = eval(to_light, to_view, ch, true);
        Vec3 fin = -to_light, fout = to_view;
        return rotate_mueller_basis(r.M, fin, r.g.x_in, default_basis(fin), fout, r.g.x_out, default_basis(fout));
    }
};

} // namespace hpbrdf
