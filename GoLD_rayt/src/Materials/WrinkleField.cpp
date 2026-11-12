// The construction of the wrinkle tile: the pore graph, its simulation and the drawing of the height map (Materials/WrinkleField.hpp).
#include "pch.h"
#include "Materials/WrinkleField.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>

namespace rayt {

    std::vector<WrinkleField::Level> WrinkleField::skinDefault() {
        std::vector<Level> v(3);
        v[0].spacingMm = Real(1.0);  v[0].angleDeg = Real(0);   v[0].base = Real(0.35); v[0].distance = Real(0.7);
        v[0].continuation = Real(0.5); v[0].similarity = Real(0.4); v[0].crossing = Real(0.6); v[0].deposit = Real(0.35);
        v[0].width = Real(0.13); v[0].poreWidth = Real(0.18); v[0].poreMin = Real(0.5); v[0].depthUm = Real(40); v[0].wobble = Real(0.2);
        v[1].spacingMm = Real(0.55); v[1].angleDeg = Real(60);  v[1].base = Real(0.5);  v[1].distance = Real(0.7);
        v[1].continuation = Real(0.4); v[1].similarity = Real(0.4); v[1].crossing = Real(0.5); v[1].deposit = Real(0.3);
        v[1].width = Real(0.16); v[1].poreWidth = Real(0.2); v[1].poreMin = Real(0.5); v[1].depthUm = Real(28); v[1].wobble = Real(0.2);
        v[2].spacingMm = Real(0.3);  v[2].angleDeg = Real(120); v[2].base = Real(1.0);  v[2].distance = Real(0.7);
        v[2].continuation = Real(0.3); v[2].similarity = Real(0.3); v[2].crossing = Real(0.5); v[2].deposit = Real(0.3);
        v[2].width = Real(0.3); v[2].poreWidth = Real(0.3); v[2].poreMin = Real(0.3); v[2].depthUm = Real(24); v[2].wobble = Real(0.25);
        v[2].bulge = false; v[2].darkPores = false;
        return v;
    }

    WrinkleField::WrinkleField(const Params& p) : m_p(p) {
        if (m_p.levels.empty()) m_p.levels = skinDefault();
        m_p.resolution = std::max(m_p.resolution, 64);
        build();
    }

    // ------------------------------------------------------------------------------------------------------------------------------
    // torus geometry (the tile is periodic)
    // ------------------------------------------------------------------------------------------------------------------------------
    Real WrinkleField::wrapDelta(Real d) const { const Real L = m_p.tileMm; d = std::fmod(d, L); if (d > L / 2) d -= L; if (d < -L / 2) d += L; return d; }
    WrinkleField::V2 WrinkleField::delta(const V2& a, const V2& b) const { return V2{ wrapDelta(b.x - a.x), wrapDelta(b.y - a.y) }; }   // from a to b, the shortest way
    Real WrinkleField::angDiff2Pi(Real a, Real b) {   // |a - b| wrapped to [0, pi]
        Real d = std::fmod(std::abs(a - b), Real(2) * constants::PI);
        return d > constants::PI ? Real(2) * constants::PI - d : d;
    }
    Real WrinkleField::angDiffPi(Real a, Real b) {    // |a - b| as undirected lines, [0, pi / 2]
        const Real d = angDiff2Pi(a, b);
        return std::min(d, constants::PI - d);
    }

    /// A slow, periodic noise in 0..1 (a few plane waves with whole numbers of periods per tile): where pores are dense.
    Real WrinkleField::density(const V2& q) const {
        Real s = 0;
        for (const auto& w : m_waves) s += std::cos(Real(2) * constants::PI * (w[0] * q.x + w[1] * q.y) / m_p.tileMm + w[2]);
        return std::clamp(Real(0.5) + s / Real(2) / Real(std::max<size_t>(m_waves.size(), 1)) * Real(1.6), Real(0), Real(1));
    }

    /// Poisson-disk sampling on the torus with minimum distance r; the nodes already in `nodes` are kept.
    void WrinkleField::addNodes(std::vector<V2>& nodes, Real r, std::mt19937& rng) const {
        std::uniform_real_distribution<Real> U(0, 1);
        const Real L = m_p.tileMm;
        const int G = std::max(1, static_cast<int>(std::floor(L / (r / Real(std::sqrt(2.0))))));
        const Real cell = L / Real(G);
        std::vector<std::vector<int>> grid(static_cast<size_t>(G) * G);
        auto cellOf = [&](const V2& p) { return std::pair<int, int>{ std::min(G - 1, int(p.x / cell)), std::min(G - 1, int(p.y / cell)) }; };
        for (size_t i = 0; i < nodes.size(); ++i) { const auto c = cellOf(nodes[i]); grid[static_cast<size_t>(c.second) * G + c.first].push_back(int(i)); }
        const long tries = static_cast<long>((L / r) * (L / r) * 14);
        for (long t = 0; t < tries; ++t) {
            const V2 q{ U(rng) * L, U(rng) * L };
            if (U(rng) > Real(1) - m_p.densityVariation + m_p.densityVariation * density(q)) continue;   // thinner in some regions
            const auto c = cellOf(q);
            bool ok = true;
            for (int dy = -2; dy <= 2 && ok; ++dy)
                for (int dx = -2; dx <= 2 && ok; ++dx) {
                    const int cx = (c.first + dx + G * 4) % G, cy = (c.second + dy + G * 4) % G;
                    for (int idx : grid[static_cast<size_t>(cy) * G + cx]) {
                        const V2 d = delta(nodes[static_cast<size_t>(idx)], q);
                        if (d.x * d.x + d.y * d.y < r * r) { ok = false; break; }
                    }
                }
            if (!ok) continue;
            nodes.push_back(q);
            grid[static_cast<size_t>(c.second) * G + c.first].push_back(int(nodes.size()) - 1);
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------------
    // the graph and its simulation
    // ------------------------------------------------------------------------------------------------------------------------------
    WrinkleField::Graph WrinkleField::makeGraph(const std::vector<V2>& nodes, int k, Real spacing) const {
        Graph g;
        g.nodes = nodes;
        const int N = int(nodes.size());
        g.adj.assign(static_cast<size_t>(N), {});
        std::vector<std::pair<Real, int>> cand;   // (not named "near": a macro of the Windows headers)
        std::vector<std::pair<int, int>> keys;
        for (int i = 0; i < N; ++i) {
            cand.clear();
            for (int j = 0; j < N; ++j) {
                if (j == i) continue;
                const V2 d = delta(nodes[static_cast<size_t>(i)], nodes[static_cast<size_t>(j)]);
                const Real d2 = d.x * d.x + d.y * d.y;
                if (d2 < spacing * spacing * Real(9)) cand.emplace_back(d2, j);
            }
            const int kk = std::min<int>(k, int(cand.size()));
            std::partial_sort(cand.begin(), cand.begin() + kk, cand.end());
            for (int a = 0; a < kk; ++a) keys.emplace_back(std::min(i, cand[static_cast<size_t>(a)].second), std::max(i, cand[static_cast<size_t>(a)].second));
        }
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
        for (const auto& kv : keys) {
            const int e = int(g.edges.size());
            g.edges.push_back({ kv.first, kv.second });
            const V2 d = delta(nodes[static_cast<size_t>(kv.first)], nodes[static_cast<size_t>(kv.second)]);
            g.edgeVec.push_back(d);
            const Real th = std::atan2(d.y, d.x), len = std::sqrt(d.x * d.x + d.y * d.y);
            g.adj[static_cast<size_t>(kv.first)].push_back({ e, kv.second, th, len });
            g.adj[static_cast<size_t>(kv.second)].push_back({ e, kv.first, th + constants::PI, len });
        }
        // crossings: bin the edge midpoints, test the neighbours
        const int E = int(g.edges.size());
        g.crossing.assign(static_cast<size_t>(E), {});
        Real maxLen = spacing;
        for (const V2& d : g.edgeVec) maxLen = std::max(maxLen, std::sqrt(d.x * d.x + d.y * d.y));
        const Real L = m_p.tileMm;
        const int G = std::max(1, int(L / (maxLen * Real(1.0))));
        const Real cell = L / Real(G);
        std::vector<std::vector<int>> grid(static_cast<size_t>(G) * G);
        std::vector<V2> mid(static_cast<size_t>(E));
        for (int e = 0; e < E; ++e) {
            const V2& a = g.nodes[static_cast<size_t>(g.edges[static_cast<size_t>(e)][0])];
            const V2& d = g.edgeVec[static_cast<size_t>(e)];
            V2 m{ std::fmod(a.x + d.x / 2 + L * 2, L), std::fmod(a.y + d.y / 2 + L * 2, L) };
            mid[static_cast<size_t>(e)] = m;
            grid[static_cast<size_t>(std::min(G - 1, int(m.y / cell))) * G + std::min(G - 1, int(m.x / cell))].push_back(e);
        }
        auto cross2 = [](const V2& a, const V2& b) { return a.x * b.y - a.y * b.x; };
        for (int e = 0; e < E; ++e) {
            const int gx = std::min(G - 1, int(mid[static_cast<size_t>(e)].x / cell)), gy = std::min(G - 1, int(mid[static_cast<size_t>(e)].y / cell));
            const V2 p0 = g.nodes[static_cast<size_t>(g.edges[static_cast<size_t>(e)][0])];
            const V2 r = g.edgeVec[static_cast<size_t>(e)];
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    const int cx = (gx + dx + G) % G, cy = (gy + dy + G) % G;
                    for (int f : grid[static_cast<size_t>(cy) * G + cx]) {
                        if (f <= e) continue;
                        const auto& ee = g.edges[static_cast<size_t>(e)];
                        const auto& ff = g.edges[static_cast<size_t>(f)];
                        if (ee[0] == ff[0] || ee[0] == ff[1] || ee[1] == ff[0] || ee[1] == ff[1]) continue;
                        const V2 q0 = [&] { const V2 d = delta(p0, g.nodes[static_cast<size_t>(ff[0])]); return d; }();   // relative to p0
                        const V2 s = g.edgeVec[static_cast<size_t>(f)];
                        const Real den = cross2(r, s);
                        if (std::abs(den) < Real(1e-12)) continue;
                        const Real t = cross2(q0, s) / den, u = cross2(q0, r) / den;
                        if (t > Real(1e-6) && t < Real(1) - Real(1e-6) && u > Real(1e-6) && u < Real(1) - Real(1e-6)) {
                            g.crossing[static_cast<size_t>(e)].push_back(f);
                            g.crossing[static_cast<size_t>(f)].push_back(e);
                        }
                    }
                }
        }
        return g;
    }

    std::vector<Real> WrinkleField::simulate(const Graph& g, const Level& L, std::mt19937& rng) const {
        const int E = int(g.edges.size()), N = int(g.nodes.size());
        std::vector<Real> W(static_cast<size_t>(E), Real(0));
        if (E == 0) return W;
        std::uniform_real_distribution<Real> U(0, 1);
        const Real aTheta = L.angleDeg * constants::PI / Real(180);
        const Real sigma = Real(0.3), dt = Real(0.03), decay = Real(0.5);
        auto gauss = [&](Real x) { return std::exp(-x * x / (Real(2) * sigma * sigma)); };
        std::vector<Real> score;
        for (int step = 0; step < m_p.steps; ++step) {
            for (int ev = 0; ev < E; ++ev) {
                const int v = std::min(N - 1, int(U(rng) * Real(N)));
                const auto& hs = g.adj[static_cast<size_t>(v)];
                if (hs.empty()) continue;
                Real dmin = hs[0].dist;
                for (const auto& h : hs) dmin = std::min(dmin, h.dist);
                score.assign(hs.size(), 0);
                Real total = 0;
                for (size_t a = 0; a < hs.size(); ++a) {
                    const HalfEdge& h = hs[a];
                    Real pcont = 0;
                    for (size_t b = 0; b < hs.size(); ++b) {
                        if (b == a) continue;
                        const Real w = W[static_cast<size_t>(hs[b].edge)];
                        pcont += w * (L.continuation * gauss(angDiff2Pi(h.theta, hs[b].theta + constants::PI)) - L.similarity * gauss(angDiff2Pi(h.theta, hs[b].theta)));
                    }
                    Real pcross = 0;
                    for (int f : g.crossing[static_cast<size_t>(h.edge)]) pcross += W[static_cast<size_t>(f)];
                    const Real ppref = (L.base + gauss(angDiffPi(h.theta, aTheta))) / (Real(1) + L.base);
                    const Real s = Real(1) - L.distance * (h.dist / dmin - Real(1)) + pcont + m_p.preference * ppref - L.crossing * pcross;
                    score[a] = std::pow(std::max(s, Real(1e-4)), m_p.selectivity);
                    total += score[a];
                }
                Real pick = U(rng) * total;
                size_t chosen = hs.size() - 1;
                for (size_t a = 0; a < hs.size(); ++a) { pick -= score[a]; if (pick <= 0) { chosen = a; break; } }
                Real& w = W[static_cast<size_t>(hs[chosen].edge)];
                w += dt * L.deposit * (Real(1) - w);
            }
            for (Real& w : W) w *= Real(1) - dt * decay;   // N / |E| = 1
        }
        return W;
    }

    // ------------------------------------------------------------------------------------------------------------------------------
    // drawing
    // ------------------------------------------------------------------------------------------------------------------------------
    /// The profile of a groove: valley (positive depth) and bulge (positive height) at the distance x in units of the width.
    void WrinkleField::profile(Real x, bool bulge, Real& valley, Real& ridge) {
        valley = 0; ridge = 0;
        if (bulge) {
            if (x < 1) valley = std::max(Real(0), -(x - 3) * (x - 3) * (x - 1) / Real(9));
            else if (x < 3) ridge = std::max(Real(0), (x - 3) * (x - 3) * (x - 1) / Real(9));
        }
        else if (x < Real(1.5)) { const Real t = Real(1) - Real(2) * x / Real(3); valley = t * t * t; }
    }

    Real WrinkleField::hash1(std::uint32_t a, std::uint32_t b) {
        std::uint32_t h = a * 0x9e3779b1u ^ b * 0x85ebca6bu;
        h ^= h >> 15; h *= 0x2c1b3c6du; h ^= h >> 12; h *= 0x297a2d39u; h ^= h >> 15;
        return Real(h) / Real(4294967296.0);
    }
    /// 1-D value noise in [-1, 1].
    Real WrinkleField::noise1(Real x, std::uint32_t seed) {
        const Real f = std::floor(x), t = x - f;
        const Real s = t * t * (Real(3) - Real(2) * t);
        const auto i = static_cast<std::uint32_t>(static_cast<std::int64_t>(f));
        return (hash1(i, seed) * (1 - s) + hash1(i + 1, seed) * s) * 2 - 1;
    }

    /// Adds the fragments of one primitive: `local` holds its valley / ridge (max over its pieces) on a pixel box.
    void WrinkleField::flush(Canvas& cv, int x0, int y0, int w, int h, const std::vector<float>& lv, const std::vector<float>& lb, Real scale) const {
        const int R = m_p.resolution;
        const Real beta = Real(20);
        for (int j = 0; j < h; ++j)
            for (int i = 0; i < w; ++i) {
                const float v = lv[static_cast<size_t>(j) * w + i], b = lb[static_cast<size_t>(j) * w + i];
                if (v <= 1e-4f && b <= 1e-4f) continue;
                const int px = ((x0 + i) % R + R) % R, py = ((y0 + j) % R + R) % R;
                const size_t id = static_cast<size_t>(py) * R + px;
                cv.accV[id] += float(std::exp(beta * scale * v));
                cv.accB[id] += float(std::exp(beta * scale * b));
                if (cv.count[id] < 255) ++cv.count[id];
            }
    }

    void WrinkleField::drawEdge(Canvas& cv, const V2& a, const V2& d, Real strength, const Level& L, Real depthNorm, std::uint32_t seed) const {
        const int R = m_p.resolution;
        const Real pxMm = m_p.tileMm / Real(R);
        const Real len = std::sqrt(d.x * d.x + d.y * d.y);
        if (len <= 0) return;
        const Real half = std::max(L.width * L.spacingMm * strength, pxMm * Real(0.9));   // width of the valley (mm) = x of 1 (not below a texel)
        const Real reach = (L.bulge ? Real(3) : Real(1.5)) * half + pxMm * Real(1.5);
        const int pieces = std::max(2, int(len / (pxMm * Real(3))));
        // the wobbling center line
        std::vector<V2> pts(static_cast<size_t>(pieces) + 1);
        const V2 perp{ -d.y / len, d.x / len };
        const Real amp = L.wobble * L.spacingMm;
        for (int i = 0; i <= pieces; ++i) {
            const Real t = Real(i) / Real(pieces);
            const Real o = amp * std::sin(constants::PI * t) * noise1(t * len / L.spacingMm * L.wobbleFreq, seed);
            pts[static_cast<size_t>(i)] = V2{ a.x + d.x * t + perp.x * o, a.y + d.y * t + perp.y * o };
        }
        Real minx = pts[0].x, maxx = pts[0].x, miny = pts[0].y, maxy = pts[0].y;
        for (const V2& p : pts) { minx = std::min(minx, p.x); maxx = std::max(maxx, p.x); miny = std::min(miny, p.y); maxy = std::max(maxy, p.y); }
        const int x0 = int(std::floor((minx - reach) / pxMm)), x1 = int(std::ceil((maxx + reach) / pxMm));
        const int y0 = int(std::floor((miny - reach) / pxMm)), y1 = int(std::ceil((maxy + reach) / pxMm));
        const int w = x1 - x0 + 1, h = y1 - y0 + 1;
        if (w <= 0 || h <= 0 || w > R || h > R) return;
        std::vector<float> lv(static_cast<size_t>(w) * h, 0.f), lb(static_cast<size_t>(w) * h, 0.f);
        for (int s = 0; s < pieces; ++s) {
            const V2& p = pts[static_cast<size_t>(s)];
            const V2 q{ pts[static_cast<size_t>(s) + 1].x - p.x, pts[static_cast<size_t>(s) + 1].y - p.y };
            const Real q2 = q.x * q.x + q.y * q.y;
            const int sx0 = int(std::floor((std::min(p.x, p.x + q.x) - reach) / pxMm)) - x0, sx1 = int(std::ceil((std::max(p.x, p.x + q.x) + reach) / pxMm)) - x0;
            const int sy0 = int(std::floor((std::min(p.y, p.y + q.y) - reach) / pxMm)) - y0, sy1 = int(std::ceil((std::max(p.y, p.y + q.y) + reach) / pxMm)) - y0;
            for (int j = std::max(0, sy0); j <= std::min(h - 1, sy1); ++j)
                for (int i = std::max(0, sx0); i <= std::min(w - 1, sx1); ++i) {
                    const Real X = (Real(x0 + i) + Real(0.5)) * pxMm - p.x, Y = (Real(y0 + j) + Real(0.5)) * pxMm - p.y;
                    const Real t = q2 > 0 ? std::clamp((X * q.x + Y * q.y) / q2, Real(0), Real(1)) : Real(0);
                    const Real dx = X - q.x * t, dy = Y - q.y * t;
                    const Real dist = std::sqrt(dx * dx + dy * dy);
                    if (dist > reach) continue;
                    Real v, r;
                    profile(dist / half, L.bulge, v, r);
                    float& lvv = lv[static_cast<size_t>(j) * w + i];
                    float& lbb = lb[static_cast<size_t>(j) * w + i];
                    lvv = std::max(lvv, float(v * strength));
                    lbb = std::max(lbb, float(r * strength));
                }
        }
        flush(cv, x0, y0, w, h, lv, lb, depthNorm);
    }

    void WrinkleField::drawPore(Canvas& cv, const V2& c, Real strength, const Level& L, Real depthNorm) const {
        const int R = m_p.resolution;
        const Real pxMm = m_p.tileMm / Real(R);
        const Real rad = std::max(L.poreWidth * L.spacingMm * strength, pxMm * Real(0.8));
        const Real skew = Real(1) + Real(2) * L.poreSkew;
        const Real ca = std::cos(L.angleDeg * constants::PI / Real(180)), sa = std::sin(L.angleDeg * constants::PI / Real(180));
        const Real reach = (L.bulge ? Real(3) : Real(1.5)) * rad * skew + pxMm;
        const int x0 = int(std::floor((c.x - reach) / pxMm)), x1 = int(std::ceil((c.x + reach) / pxMm));
        const int y0 = int(std::floor((c.y - reach) / pxMm)), y1 = int(std::ceil((c.y + reach) / pxMm));
        const int w = x1 - x0 + 1, h = y1 - y0 + 1;
        if (w <= 0 || h <= 0 || w > R || h > R) return;
        std::vector<float> lv(static_cast<size_t>(w) * h, 0.f), lb(static_cast<size_t>(w) * h, 0.f);
        for (int j = 0; j < h; ++j)
            for (int i = 0; i < w; ++i) {
                const Real X = (Real(x0 + i) + Real(0.5)) * pxMm - c.x, Y = (Real(y0 + j) + Real(0.5)) * pxMm - c.y;
                const Real u = (X * ca + Y * sa) / skew, v = -X * sa + Y * ca;   // elongated along the primary orientation
                const Real dist = std::sqrt(u * u + v * v);
                Real val, rid;
                profile(dist / rad, L.bulge, val, rid);
                lv[static_cast<size_t>(j) * w + i] = float(val * strength);
                lb[static_cast<size_t>(j) * w + i] = float(rid * strength);
                if (val > 0 && L.darkPores) {   // the pore mask
                    const int px = ((x0 + i) % R + R) % R, py = ((y0 + j) % R + R) % R;
                    float& pm = cv.pore[static_cast<size_t>(py) * R + px];
                    pm = std::max(pm, float(val * std::min(Real(1), strength * Real(1.5))));
                }
            }
        flush(cv, x0, y0, w, h, lv, lb, depthNorm);
    }

    // ------------------------------------------------------------------------------------------------------------------------------
    void WrinkleField::build() {
        const auto t0 = std::chrono::steady_clock::now();
        const int R = m_p.resolution;
        std::mt19937 rng(m_p.seed);
        {
            std::uniform_real_distribution<Real> W01(0, 1);
            for (int i = 0; i < 6; ++i) {
                const Real ang = W01(rng) * Real(2) * constants::PI, k = Real(1) + std::floor(W01(rng) * Real(3));   // 1..3 periods per tile
                m_waves.push_back({ std::round(k * std::cos(ang)), std::round(k * std::sin(ang)), W01(rng) * Real(2) * constants::PI });
            }
        }
        Canvas cv;
        cv.accV.assign(static_cast<size_t>(R) * R, 0.f);
        cv.accB.assign(static_cast<size_t>(R) * R, 0.f);
        cv.pore.assign(static_cast<size_t>(R) * R, 0.f);
        cv.count.assign(static_cast<size_t>(R) * R, 0);
        Real depthMax = 0;
        for (const Level& L : m_p.levels) depthMax = std::max(depthMax, L.depthUm);
        depthMax *= m_p.depthScale;
        std::vector<V2> nodes;
        int li = 0;
        for (const Level& L : m_p.levels) {
            addNodes(nodes, L.spacingMm * Real(0.8), rng);
            const Graph g = makeGraph(nodes, 5, L.spacingMm);
            std::vector<Real> W = simulate(g, L, rng);
            Real wmax = Real(1e-6);
            for (Real w : W) wmax = std::max(wmax, w);
            const Real depthNorm = (L.depthUm * m_p.depthScale) / std::max(depthMax, Real(1e-6));
            // the edges
            int drawn = 0;
            for (size_t e = 0; e < g.edges.size(); ++e) {
                const Real s = W[e] / wmax;
                if (s < m_p.carveThreshold) continue;
                ++drawn;
                drawEdge(cv, g.nodes[static_cast<size_t>(g.edges[e][0])], g.edgeVec[e], s, L, depthNorm, std::uint32_t(m_p.seed * 7919u + li * 104729u + e));
            }
            // the pores: every node, deeper where wrinkles meet
            for (size_t v = 0; v < g.nodes.size(); ++v) {
                Real wm = 0, ws = 0;
                for (const auto& h : g.adj[v]) { const Real s = W[static_cast<size_t>(h.edge)] / wmax; wm = std::max(wm, s); ws += s; }
                const Real wv = std::min(Real(1), std::max(wm + L.poreMin, wm + L.poreBlend * (ws - wm)));
                drawPore(cv, g.nodes[v], wv, L, depthNorm);
            }
            if (m_p.levels.size() > 0 && !m_quiet)
                std::cout << "[Wrinkles] level " << li << ": spacing " << L.spacingMm << " mm, " << g.nodes.size() << " pores, " << g.edges.size() << " edges, " << drawn << " carved\n";
            ++li;
        }
        // MellowMax of the fragments (padded with zero-depth ones up to Nf = 16), valley and ridge apart
        m_height.assign(static_cast<size_t>(R) * R, 0.f);
        const Real Nf = Real(16), beta = Real(20);
        for (size_t i = 0; i < m_height.size(); ++i) {
            const Real n = std::min<Real>(Real(cv.count[i]), Nf);
            const Real mv = std::log((Nf - n + cv.accV[i]) / Nf) / beta, mb = std::log((Nf - n + cv.accB[i]) / Nf) / beta;
            m_height[i] = float((mb - mv) * depthMax * Real(1e-6));   // metres
        }
        m_pore = std::move(cv.pore);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (!m_quiet) std::cout << "[Wrinkles] " << m_p.levels.size() << " levels, tile " << m_p.tileMm << " mm, " << R << "^2 texels, built in " << sec << " s\n";
    }

} // namespace rayt
