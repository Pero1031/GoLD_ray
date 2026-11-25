// GPU checks (--scene gpu, docs/GPU.md): the kernel code (include/Kernel) on the CPU against the engine's functions,
// and on the GPU against the CPU. Without the GPU module or a CUDA device the GPU part prints SKIP.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Frame.hpp"
#include "Gpu/GpuDevice.hpp"
#include "Kernel/SelfTest.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace rayt::cli {

namespace {

    bool Report(const std::string& name, bool ok, const std::string& detail) {
        std::cout << "[GPU]    " << name << "\n             " << detail << (ok ? "   OK" : "   FAIL") << "\n";
        return ok;
    }

    std::string Fmt(double v, int prec = 3) {
        std::ostringstream s;
        s << std::setprecision(prec) << v;
        return s.str();
    }

    /// Cases for evalMathCase: random ones plus edges (u = 0, poles of the basis, -0, total internal reflection, cos < 0).
    std::vector<kernel::MathCase> MakeCases(int n) {
        std::vector<kernel::MathCase> cases;
        kernel::Pcg32 rng;
        rng.seed(20261126, 7);
        auto add = [&](float u1, float u2, float nx, float ny, float nz, float cosT, float etaI, float etaT) {
            kernel::MathCase c{};
            c.seed = 0x5eedull + uint64_t(cases.size()) * 977u;
            c.pixel = uint32_t(cases.size() * 31u);
            c.sample = uint32_t(cases.size() % 1024u);
            c.u1 = u1; c.u2 = u2; c.nx = nx; c.ny = ny; c.nz = nz; c.cosTheta = cosT; c.etaI = etaI; c.etaT = etaT;
            cases.push_back(c);
        };
        add(0, 0, 0, 0, 1, 1, 1, 1.5f);
        add(0.5f, 0.999999f, 0, 0, -1, 0, 1.5f, 1);
        add(0.25f, 0.5f, 0.6f, 0.8f, -0.0f, -0.3f, 1, 1.33f);   // n.z = -0: copysign picks the lower branch
        add(0.75f, 0.25f, 1, 0, 0, 0.2f, 1.33f, 1);
        const float etas[4][2] = { { 1, 1.5f }, { 1.5f, 1 }, { 1, 1.33f }, { 1.33f, 1 } };
        while (int(cases.size()) < n) {
            const float nz = 1 - 2 * rng.nextFloat(), phi = 6.2831853f * rng.nextFloat(), r = std::sqrt(std::max(0.0f, 1 - nz * nz));
            const int e = int(rng.nextUint() % 4u);
            add(rng.nextFloat(), rng.nextFloat(), r * std::cos(phi), r * std::sin(phi), nz, -0.2f + 1.2f * rng.nextFloat(), etas[e][0], etas[e][1]);
        }
        return cases;
    }

    const char* ValueName(int k) {
        static const char* names[kernel::kMathValueCount] = { "rng float 0", "rng float 1", "hemisphere x", "hemisphere y", "hemisphere z",
            "sphere x", "sphere y", "sphere z", "basis s.x", "basis s.y", "basis s.z", "basis t.x", "basis t.y", "basis t.z", "Fresnel" };
        return names[k];
    }

    /// The engine's functions for the same case, in double: what the kernel's double values must equal.
    void EngineValues(const kernel::MathCase& c, double* v) {
        const Point2 u(Real(c.u1), Real(c.u2));
        const Vector3 h = sampling::CosineSampleHemisphere(u);
        const Vector3 s = sampling::UniformSampleSphere(u);
        Vector3 n(Real(c.nx), Real(c.ny), Real(c.nz));
        n = n / std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);   // as the kernel normalizes (glm::normalize rounds differently)
        Vector3 bs, bt;
        frame::makeOrthonormalBasis(n, bs, bt);
        const double vals[] = { h.x, h.y, h.z, s.x, s.y, s.z, bs.x, bs.y, bs.z, bt.x, bt.y, bt.z,
            fresnel::fresnelDielectric(Real(c.cosTheta), Real(c.etaI), Real(c.etaT)) };
        for (int k = kernel::kHemiX; k < kernel::kMathValueCount; ++k) v[k] = vals[k - kernel::kHemiX];
    }

} // namespace

int RunGpuTests(const CliOptions& /*opt*/) {
    bool ok = true;
    const int n = 4096;
    const std::vector<kernel::MathCase> cases = MakeCases(n);
    std::vector<kernel::MathOut> cpu(cases.size());
    for (size_t i = 0; i < cases.size(); ++i) kernel::evalMathCase(cases[i], cpu[i]);

    // 1. The kernel's functions on the CPU in double are the engine's functions (same formulas, same order).
    {
        double maxDiff = 0;
        size_t differing = 0, nans = 0;
        for (size_t i = 0; i < cases.size(); ++i) {
            double e[kernel::kMathValueCount];
            EngineValues(cases[i], e);
            bool same = true;
            for (int k = kernel::kHemiX; k < kernel::kMathValueCount; ++k) {
                const double d = std::abs(cpu[i].d[k] - e[k]);
                if (std::isnan(d)) ++nans; else maxDiff = std::max(maxDiff, d);
                if (!(cpu[i].d[k] == e[k])) same = false;
            }
            if (!same) ++differing;
        }
        ok &= Report("kernel (CPU, double) = engine functions: cosine hemisphere, uniform sphere, orthonormal basis, dielectric Fresnel", differing == 0,
            std::to_string(cases.size()) + " cases, " + std::to_string(differing) + " not bit-identical (max |diff| " + Fmt(maxDiff) + ", " +
            std::to_string(nans) + " NaN)");
    }

    // 2. The GPU module.
    std::string why;
    const gpu::Module* mod = gpu::LoadModule(why);
    if (!mod) {
        std::cout << "[GPU]    SKIP the GPU part: " << why << "\n";
        return ok ? 0 : 2;
    }
    RaytGpuInfo info{};
    if (mod->query(&info) != 0) {
        std::cout << "[GPU]    SKIP the GPU part: " << info.error << " (module " << mod->path << ")\n";
        return ok ? 0 : 2;
    }
    std::cout << "[GPU]    device 0: " << info.name << ", compute capability " << info.computeMajor << "." << info.computeMinor << ", "
              << info.multiprocessors << " multiprocessors, " << Fmt(double(info.totalMemory) / (1024.0 * 1024.0 * 1024.0), 3) << " GB; driver supports CUDA "
              << info.driverVersion / 1000 << "." << (info.driverVersion % 1000) / 10 << ", module built with CUDA " << info.runtimeVersion / 1000 << "."
              << (info.runtimeVersion % 1000) / 10 << "\n";

    std::vector<kernel::MathOut> gpu(cases.size());
    char err[512] = {};
    if (mod->mathTest(cases.data(), gpu.data(), int(cases.size()), err, int(sizeof(err))) != 0)
        return Report("evalMathCase on the GPU", false, err) ? 0 : 2;

    // 3. Random numbers: integer arithmetic, so the same bits.
    {
        size_t differing = 0;
        for (size_t i = 0; i < cases.size(); ++i) {
            bool same = true;
            for (int k = 0; k < 4; ++k) same &= cpu[i].rng[k] == gpu[i].rng[k];
            for (int k = kernel::kRngFloat0; k <= kernel::kRngFloat1; ++k) same &= cpu[i].f[k] == gpu[i].f[k] && cpu[i].d[k] == gpu[i].d[k];
            if (!same) ++differing;
        }
        ok &= Report("random numbers (PCG32 streams keyed by seed, pixel, sample): GPU = CPU bit for bit", differing == 0,
            std::to_string(cases.size()) + " streams, " + std::to_string(differing) + " differ");
    }

    // 4. Floating point: the same formulas, but CUDA's sinf / cosf and its fused multiply-adds round differently from MSVC.
    auto compare = [&](bool useFloat, double tol) {
        double maxDiff = 0;
        size_t identical = 0, values = 0, nans = 0, worstCase = 0;
        int worstValue = kernel::kHemiX;
        for (size_t i = 0; i < cases.size(); ++i)
            for (int k = kernel::kHemiX; k < kernel::kMathValueCount; ++k) {
                const double a = useFloat ? double(cpu[i].f[k]) : cpu[i].d[k];
                const double b = useFloat ? double(gpu[i].f[k]) : gpu[i].d[k];
                if (std::isnan(a) || std::isnan(b)) ++nans;
                else if (std::abs(a - b) > maxDiff) { maxDiff = std::abs(a - b); worstCase = i; worstValue = k; }
                identical += a == b;
                ++values;
            }
        const kernel::MathCase& w = cases[worstCase];
        std::ostringstream where;
        where << ValueName(worstValue) << " of case " << worstCase << " (u " << w.u1 << ", " << w.u2 << "; n " << w.nx << ", " << w.ny << ", " << w.nz
              << "; cos " << w.cosTheta << ", eta " << w.etaI << " -> " << w.etaT << ")";
        return std::make_tuple(maxDiff <= tol && nans == 0, maxDiff, double(identical) / double(values), where.str());
    };
    {
        // Some values amplify rounding (the Fresnel term just past the critical angle: cos theta_t = sqrt(small)), so a fixed tolerance on
        // |GPU - CPU| is not the right test. Instead: the GPU's float values are as close to the double values as the CPU's float values.
        double errCpu = 0, errGpu = 0;
        for (size_t i = 0; i < cases.size(); ++i)
            for (int k = kernel::kHemiX; k < kernel::kMathValueCount; ++k) {
                errCpu = std::max(errCpu, std::abs(double(cpu[i].f[k]) - cpu[i].d[k]));
                errGpu = std::max(errGpu, std::abs(double(gpu[i].f[k]) - cpu[i].d[k]));
            }
        const auto [unused, maxDiff, same, where] = compare(true, 1.0);
        (void)unused;
        ok &= Report("float: the GPU's float values are as accurate as the CPU's (sampling, basis, Fresnel)", errGpu <= 2 * errCpu + 1e-7,
            "max |float - double|: GPU " + Fmt(errGpu) + ", CPU " + Fmt(errCpu) + " (GPU must be <= 2 x CPU); max |GPU - CPU| " + Fmt(maxDiff) + " at " + where +
            "; " + Fmt(100 * same, 4) + " % of the values bit-identical");
    }
    {
        const auto [pass, maxDiff, same, where] = compare(false, 1e-12);
        ok &= Report("double: GPU = CPU to rounding", pass,
            "max |diff| " + Fmt(maxDiff) + " (tolerance 1e-12) at " + where + ", " + Fmt(100 * same, 4) + " % of the values bit-identical");
    }
    return ok ? 0 : 2;
}

} // namespace rayt::cli
