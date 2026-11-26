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
    // The CPU computes the float values (as the GPU does) and the double reference (which only the CPU computes).
    std::vector<kernel::MathOut> cpu(cases.size());
    std::vector<std::array<double, kernel::kMathValueCount>> ref(cases.size());
    for (size_t i = 0; i < cases.size(); ++i) {
        kernel::evalMathCase(cases[i], cpu[i]);
        kernel::evalMathCaseDouble(cases[i], ref[i].data());
    }

    // 1. The kernel's functions on the CPU in double are the engine's functions (same formulas, same order).
    {
        double maxDiff = 0;
        size_t differing = 0, nans = 0;
        for (size_t i = 0; i < cases.size(); ++i) {
            double e[kernel::kMathValueCount];
            EngineValues(cases[i], e);
            bool same = true;
            for (int k = kernel::kHemiX; k < kernel::kMathValueCount; ++k) {
                const double d = std::abs(ref[i][k] - e[k]);
                if (std::isnan(d)) ++nans; else maxDiff = std::max(maxDiff, d);
                if (!(ref[i][k] == e[k])) same = false;
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
        std::cout << "[GPU]    SKIP the GPU part: " << why << std::endl;
        return ok ? 0 : 2;
    }
    RaytGpuInfo info{};
    if (mod->query(&info) != 0) {
        std::cout << "[GPU]    SKIP the GPU part: " << info.error << " (module " << mod->path << ")" << std::endl;
        return ok ? 0 : 2;
    }
    std::cout << "[GPU]    device 0: " << info.name << ", compute capability " << info.computeMajor << "." << info.computeMinor << ", "
              << info.multiprocessors << " multiprocessors, " << Fmt(double(info.totalMemory) / (1024.0 * 1024.0 * 1024.0), 3) << " GB; driver supports CUDA "
              << info.driverVersion / 1000 << "." << (info.driverVersion % 1000) / 10 << "; module built with CUDA " << info.runtimeVersion / 1000 << "."
              << (info.runtimeVersion % 1000) / 10 << ", float only, "
              << (info.fastMath ? "fast math (--use_fast_math)" : "IEEE division / square root, no implicit FMA") << std::endl;

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
            for (int k = kernel::kRngFloat0; k <= kernel::kRngFloat1; ++k) same &= cpu[i].f[k] == gpu[i].f[k];
            if (!same) ++differing;
        }
        ok &= Report("random numbers (PCG32 streams keyed by seed, pixel, sample): GPU = CPU bit for bit", differing == 0,
            std::to_string(cases.size()) + " streams, " + std::to_string(differing) + " differ");
    }

    // 4. Float values. The default module computes the same operations in the same order as the CPU (IEEE division and square root, no
    //    implicit FMA, our own sin / cos), so its floats must be the CPU's bit for bit. A fast-math module (approximate division and
    //    square root, implicit FMAs) cannot be; it must be as close to the double reference as the CPU is. (A fixed tolerance on
    //    |GPU - CPU| would be the wrong test: the Fresnel term just past the critical angle, cos theta_t = sqrt(small), amplifies a
    //    rounding difference in its input to 1.4e-5.)
    {
        struct Group { const char* name; int first, last; };
        const Group groups[] = { { "cosine hemisphere", kernel::kHemiX, kernel::kHemiZ }, { "uniform sphere", kernel::kSphereX, kernel::kSphereZ },
            { "orthonormal basis", kernel::kBasisSX, kernel::kBasisTZ }, { "dielectric Fresnel", kernel::kFresnel, kernel::kFresnel } };
        for (const Group& g : groups) {
            double errCpu = 0, errGpu = 0, maxDiff = 0;
            size_t identical = 0, values = 0, nans = 0, worst = 0;
            for (size_t i = 0; i < cases.size(); ++i)
                for (int k = g.first; k <= g.last; ++k) {
                    const double c = cpu[i].f[k], d = gpu[i].f[k];
                    if (std::isnan(c) || std::isnan(d)) { ++nans; continue; }
                    errCpu = std::max(errCpu, std::abs(c - ref[i][k]));
                    errGpu = std::max(errGpu, std::abs(d - ref[i][k]));
                    if (std::abs(c - d) > maxDiff) { maxDiff = std::abs(c - d); worst = i; }
                    identical += c == d;
                    ++values;
                }
            std::ostringstream detail;
            detail << Fmt(100.0 * double(identical) / double(values), 4) << " % of the values bit-identical, max |GPU - CPU| " << Fmt(maxDiff) << " (case " << worst
                   << "); max |float - double|: GPU " << Fmt(errGpu) << ", CPU " << Fmt(errCpu);
            if (nans) detail << ", " << nans << " NaN";
            const bool pass = nans == 0 && (info.fastMath ? errGpu <= 2 * errCpu + 1e-7 : identical == values);
            ok &= Report(std::string("float, ") + g.name + (info.fastMath ? ": GPU (fast math) as accurate as the CPU (<= 2 x CPU + 1e-7)" : ": GPU = CPU bit for bit"),
                pass, detail.str());
        }
    }
    return ok ? 0 : 2;
}

} // namespace rayt::cli
