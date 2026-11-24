// The PROSPECT leaf model (Materials/Prospect.hpp). The steps follow the authors' R package (jbferet/prospect, prospect.R and
// calctav.R), written again from the equations.
#include "pch.h"

#include "Materials/Prospect.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <mutex>
#include <sstream>

namespace rayt::leaf {

    Real expintE1(Real x) {
        if (!(x > 0)) return std::numeric_limits<Real>::infinity();
        constexpr Real euler = Real(0.57721566490153286061);
        if (x <= 1) {   // E1(x) = -gamma - ln x - sum_k (-x)^k / (k k!)
            Real sum = 0, term = 1;
            for (int k = 1; k < 60; ++k) {
                term *= -x / Real(k);
                const Real add = term / Real(k);
                sum += add;
                if (std::abs(add) < Real(1e-17) * std::abs(sum)) break;
            }
            return -euler - std::log(x) - sum;
        }
        // Continued fraction (modified Lentz): E1(x) = e^-x / (x + 1 / (1 + 1 / (x + 2 / (1 + 2 / (x + ...)))))
        constexpr Real tiny = Real(1e-300);
        Real b = x + 1, c = Real(1) / tiny, d = Real(1) / b, h = d;
        for (int i = 1; i < 200; ++i) {
            const Real a = -Real(i) * Real(i);
            b += 2;
            d = Real(1) / (a * d + b);
            c = b + a / c;
            const Real del = c * d;
            h *= del;
            if (std::abs(del - 1) < Real(1e-16)) break;
        }
        return h * std::exp(-x);
    }

    Real tav(Real alphaDeg, Real nr) {
        const Real rd = Real(3.14159265358979323846) / Real(180);
        const Real n2 = nr * nr, np = n2 + 1, nm = n2 - 1;
        const Real a = (nr + 1) * (nr + 1) / 2;
        const Real k = -(n2 - 1) * (n2 - 1) / 4;
        const Real sa = std::sin(alphaDeg * rd);
        const Real b2 = sa * sa - np / 2;
        const Real b1 = alphaDeg == Real(90) ? Real(0) : std::sqrt(b2 * b2 + k);
        const Real b = b1 - b2;
        const Real b3 = b * b * b, a3 = a * a * a;
        const Real ts = (k * k / (6 * b3) + k / b - b / 2) - (k * k / (6 * a3) + k / a - a / 2);
        const Real tp1 = -2 * n2 * (b - a) / (np * np);
        const Real tp2 = -2 * n2 * np * std::log(b / a) / (nm * nm);
        const Real tp3 = n2 * (Real(1) / b - Real(1) / a) / 2;
        const Real tp4 = 16 * n2 * n2 * (n2 * n2 + 1) * std::log((2 * np * b - nm * nm) / (2 * np * a - nm * nm)) / (np * np * np * nm * nm);
        const Real tp5 = 16 * n2 * n2 * n2 * (Real(1) / (2 * np * b - nm * nm) - Real(1) / (2 * np * a - nm * nm)) / (np * np * np);
        return (ts + tp1 + tp2 + tp3 + tp4 + tp5) / (2 * sa * sa);
    }

    std::shared_ptr<const ProspectTable> ProspectTable::load(const std::string& path, std::string* err) {
        std::ifstream in(path);
        if (!in) {
            if (err) *err = "cannot open " + path + " (run: python tools/data/leaf_optics.py fetch)";
            return nullptr;
        }
        std::string line;
        std::getline(in, line);
        std::map<std::string, int> col;
        {
            std::istringstream hs(line);
            std::string name;
            for (int i = 0; hs >> name; ++i) col[name] = i;
        }
        const char* names[] = { "lambda", "nrefrac", "sac_chl", "sac_car", "sac_ant", "sac_brown", "sac_ewt", "sac_lma" };
        for (const char* nm : names)
            if (!col.count(nm)) {
                if (err) *err = path + ": no column \"" + nm + "\" in the header";
                return nullptr;
            }
        auto t = std::make_shared<ProspectTable>();
        Real expect = -1;
        while (std::getline(in, line)) {
            std::istringstream ls(line);
            std::vector<Real> v;
            Real x;
            while (ls >> x) v.push_back(x);
            if (v.empty()) continue;
            if (int(v.size()) < int(col.size())) {
                if (err) *err = path + ": a short row";
                return nullptr;
            }
            const Real lam = v[col["lambda"]];
            if (expect < 0) t->m_lambda0 = lam;
            else if (std::abs(lam - expect) > Real(1e-6)) {
                if (err) *err = path + ": the rows must be 1 nm apart";
                return nullptr;
            }
            expect = lam + 1;
            Row r;
            r.n = v[col["nrefrac"]];
            r.chl = v[col["sac_chl"]]; r.car = v[col["sac_car"]]; r.ant = v[col["sac_ant"]];
            r.brown = v[col["sac_brown"]]; r.ewt = v[col["sac_ewt"]]; r.lma = v[col["sac_lma"]];
            r.tav90 = tav(Real(90), r.n);
            r.tav40 = tav(Real(40), r.n);
            t->m_rows.push_back(r);
        }
        if (t->m_rows.size() < 2) {
            if (err) *err = path + ": no data";
            return nullptr;
        }
        return t;
    }

    std::shared_ptr<const ProspectTable> ProspectTable::standard(std::string* err) {
        static std::mutex lock;
        static std::shared_ptr<const ProspectTable> cached;
        static std::string cachedErr;
        static bool tried = false;
        std::lock_guard<std::mutex> g(lock);
        if (!tried) {
            tried = true;
            cached = load(standardPath(), &cachedErr);
        }
        if (!cached && err) *err = cachedErr;
        return cached;
    }

    ProspectTable::Row ProspectTable::at(Real lambdaNm) const {
        const Real x = std::clamp(lambdaNm - m_lambda0, Real(0), Real(m_rows.size() - 1));
        const size_t i = std::min(size_t(x), m_rows.size() - 2);
        const Real f = x - Real(i);
        const Row& a = m_rows[i];
        const Row& b = m_rows[i + 1];
        auto mix = [f](Real p, Real q) { return p + (q - p) * f; };
        Row r;
        r.n = mix(a.n, b.n);
        r.chl = mix(a.chl, b.chl); r.car = mix(a.car, b.car); r.ant = mix(a.ant, b.ant);
        r.brown = mix(a.brown, b.brown); r.ewt = mix(a.ewt, b.ewt); r.lma = mix(a.lma, b.lma);
        r.tav90 = mix(a.tav90, b.tav90);
        r.tav40 = mix(a.tav40, b.tav40);
        return r;
    }

    ReflTrans prospect(const ProspectTable& table, const Pigments& p, Real lambdaNm, Real alphaDeg) {
        const ProspectTable::Row c = table.at(lambdaNm);
        const Real N = std::max(p.N, Real(1));
        const Real k = (p.chl * c.chl + p.car * c.car + p.ant * c.ant + p.brown * c.brown + p.ewt * c.ewt + p.lma * c.lma) / N;

        // transmission of one plate for isotropic light (Allen et al. 1969)
        const Real tau = k > 0 ? (Real(1) - k) * std::exp(-k) + k * k * expintE1(k) : Real(1);

        // the faces
        const Real talf = alphaDeg == Real(90) ? c.tav90 : (alphaDeg == Real(40) ? c.tav40 : tav(alphaDeg, c.n));
        const Real ralf = Real(1) - talf;
        const Real t12 = c.tav90, r12 = Real(1) - t12;
        const Real t21 = t12 / (c.n * c.n), r21 = Real(1) - t21;

        // the first plate seen from the light (top) and from inside (bottom)
        Real denom = Real(1) - r21 * r21 * tau * tau;
        const Real Ta = talf * tau * t21 / denom;
        const Real Ra = ralf + r21 * tau * Ta;
        const Real t = t12 * tau * t21 / denom;
        const Real r = r12 + r21 * tau * t;

        // the other N - 1 plates (Stokes 1862, N real)
        Real Rsub, Tsub;
        if (r + t >= Real(1) - Real(1e-12)) {   // no absorption
            Tsub = t / (t + (Real(1) - t) * (N - 1));
            Rsub = Real(1) - Tsub;
        }
        else {
            const Real D = std::sqrt(std::max((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t), Real(0)));
            const Real rq = r * r, tq = t * t;
            const Real a = (1 + rq - tq + D) / (2 * r);
            const Real b = (1 - rq + tq + D) / (2 * t);
            const Real bNm1 = std::pow(b, N - 1);
            const Real bN2 = bNm1 * bNm1;
            const Real a2 = a * a;
            const Real den = a2 * bN2 - 1;
            if (!(t > Real(1e-300)) || !std::isfinite(den) || std::abs(den) < Real(1e-300)) {   // an opaque plate: nothing gets through
                Rsub = Real(1) / a;
                Tsub = 0;
            }
            else {
                Rsub = a * (bN2 - 1) / den;
                Tsub = bNm1 * (a2 - 1) / den;
            }
        }

        denom = Real(1) - Rsub * r;
        ReflTrans out;
        out.T = Ta * Tsub / denom;
        out.R = Ra + Ta * Rsub * t / denom;
        out.rFace = ralf;
        out.t12 = t12;
        out.n = c.n;
        return out;
    }

} // namespace rayt::leaf
