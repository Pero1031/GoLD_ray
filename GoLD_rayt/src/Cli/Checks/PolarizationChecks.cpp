// Polarization checks (--scene polarization): Fresnel amplitudes against closed forms, the chain of Mueller matrices and basis rotations
// against a direct propagation of the electric field in 3D, and the Stokes vector of a path tracer ray reflected by smooth glass.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/SceneParts.hpp"
#include "Color/SpectralConversion.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Polarization.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Dielectric.hpp"
#include "Render/Integrator.hpp"

namespace rayt::cli {

namespace {

using pol::Complex;
constexpr Real kPi = Real(3.14159265358979323846);

struct CVec3 { Complex v[3]; };   // an electric field: three complex components in world space

Complex dotC(const CVec3& e, const Vector3& a) { return e.v[0] * a.x + e.v[1] * a.y + e.v[2] * a.z; }
CVec3 scaled(const Vector3& a, Complex c) { return { { c * a.x, c * a.y, c * a.z } }; }
CVec3 sum(const CVec3& a, const CVec3& b) { return { { a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2] } }; }

/// The Stokes vector of a field in the basis (x, y = k x x): the definitions of Core/Polarization.hpp.
void stokesOf(const CVec3& e, const Vector3& x, const Vector3& k, Real out[4]) {
    const Complex ex = dotC(e, x), ey = dotC(e, glm::cross(k, x));
    const Complex w = ex * std::conj(ey);
    out[0] = std::norm(ex) + std::norm(ey);
    out[1] = std::norm(ex) - std::norm(ey);
    out[2] = Real(2) * std::real(w);
    out[3] = -Real(2) * std::imag(w);
}

Vector3 randomUnit() {
    for (;;) {
        const Vector3 v(Real(2) * sampling::Random() - 1, Real(2) * sampling::Random() - 1, Real(2) * sampling::Random() - 1);
        const Real l = glm::length(v);
        if (l > Real(0.1) && l <= Real(1)) return v / l;
    }
}

Vector3 perpendicularTo(const Vector3& k) {
    const Vector3 r = randomUnit();
    return glm::normalize(r - k * glm::dot(r, k));
}

Real mean3(const Vector3& v) { return (v.x + v.y + v.z) / Real(3); }

} // namespace

int RunPolarizationTests() {
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok, const std::string& detail = "") {
        std::cout << "[Polarization]   " << name << (detail.empty() ? "" : "  " + detail) << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    std::cout << "[Polarization] Stokes vectors, Mueller matrices and the Fresnel amplitudes\n";

    // (1) The amplitudes against the closed forms of a dielectric: normal incidence, Brewster's angle, the unpolarized reflectance,
    //     energy conservation, total internal reflection.
    {
        const Real n = 1.5;
        const pol::Amplitudes a0 = pol::fresnelAmplitudes(1.0, Complex(n));
        report("normal incidence: rs = (1 - n) / (1 + n) = -rp",
            std::abs(a0.rs - Complex((1 - n) / (1 + n))) < 1e-12 && std::abs(a0.rp + a0.rs) < 1e-12);

        const Real thetaB = std::atan(n);
        const pol::Amplitudes ab = pol::fresnelAmplitudes(std::cos(thetaB), Complex(n));
        report("Brewster's angle atan(n): rp = 0", std::abs(ab.rp) < 1e-12, "|rp| = " + std::to_string(std::abs(ab.rp)));

        Real worstF = 0, worstE = 0;
        for (Real deg = 0; deg < 89.5; deg += 3.5) {
            const Real ci = std::cos(deg * kPi / 180);
            const pol::Amplitudes a = pol::fresnelAmplitudes(ci, Complex(n));
            worstF = std::max(worstF, std::abs(Real(0.5) * (std::norm(a.rs) + std::norm(a.rp)) - fresnel::fresnelDielectric(ci, 1.0, n)));
            const Complex ct = std::sqrt(Complex(1) - Complex(1 - ci * ci) / Complex(n * n));
            const Real g = std::real(Complex(n) * ct) / ci;   // the energy factor of transmission
            worstE = std::max(worstE, std::abs(std::norm(a.rs) + g * std::norm(a.ts) - 1));
            worstE = std::max(worstE, std::abs(std::norm(a.rp) + g * std::norm(a.tp) - 1));
        }
        report("unpolarized reflectance = fresnelDielectric, 26 angles", worstF < 1e-12, "max difference " + std::to_string(worstF));
        report("R + T = 1 for s and p, 26 angles", worstE < 1e-12, "max difference " + std::to_string(worstE));

        // Inside the glass beyond the critical angle (41.8 deg for n = 1.5): all the light is reflected, with a phase difference.
        const pol::Amplitudes t = pol::fresnelAmplitudes(std::cos(55 * kPi / 180), Complex(1 / n));
        const Complex c = t.rs * std::conj(t.rp);
        report("total internal reflection (55 deg inside): |rs| = |rp| = 1, phase difference", std::abs(std::abs(t.rs) - 1) < 1e-12 && std::abs(std::abs(t.rp) - 1) < 1e-12 && std::abs(std::imag(c)) > 0.1,
            "arg(rs conj(rp)) = " + std::to_string(std::arg(c) * 180 / kPi) + " deg");

        // A conductor: the unpolarized reflectance of the amplitudes is the engine's conductor Fresnel.
        Real worstM = 0;
        for (Real deg = 0; deg < 89.5; deg += 7) {
            const Real ci = std::cos(deg * kPi / 180);
            const pol::Amplitudes a = pol::fresnelAmplitudes(ci, Complex(0.27, 2.95));
            worstM = std::max(worstM, std::abs(Real(0.5) * (std::norm(a.rs) + std::norm(a.rp)) - fresnel::fresnelConductor(ci, Real(0.27), Real(2.95))));
        }
        report("metal (0.27 + 2.95 i): unpolarized reflectance = fresnelConductor, 13 angles", worstM < 1e-9, "max difference " + std::to_string(worstM));
    }

    // (2) A chain of two reflections, the second on a metal, with random directions: the Stokes vector at the camera from the matrix chain
    //     (the integrator's two lines) against the electric field propagated with the Fresnel amplitudes, vector by vector.
    {
        Real worst = 0;
        int trials = 0;
        for (int trial = 0; trial < 400; ++trial) {
            Vector3 n1 = randomUnit(), n2 = randomUnit();
            const Vector3 d0 = randomUnit();
            const Vector3 wo1 = -d0;
            if (glm::dot(wo1, n1) < 0) n1 = -n1;
            const Vector3 wi1 = Real(2) * glm::dot(wo1, n1) * n1 - wo1;
            const Vector3 wo2 = -wi1;
            if (glm::dot(wo2, n2) < 0) n2 = -n2;
            const Vector3 wi2 = Real(2) * glm::dot(wo2, n2) * n2 - wo2;
            const Real c1 = glm::dot(wo1, n1), c2 = glm::dot(wo2, n2);
            if (c1 < 0.05 || c2 < 0.05) continue;   // a grazing bounce: the s axis is ill-conditioned
            ++trials;

            const pol::Amplitudes a1 = pol::fresnelAmplitudes(c1, Complex(1.5));
            const pol::Amplitudes a2 = pol::fresnelAmplitudes(c2, Complex(0.3, 2.8));
            const Vector3 s1 = pol::planeNormalS(wo1, n1), s2 = pol::planeNormalS(wo2, n2);

            // The light that arrives at the second surface, along -wi2, with a random (elliptical) polarization.
            const Vector3 kIn2 = -wi2;
            const Vector3 ea = perpendicularTo(kIn2), eb = glm::cross(kIn2, ea);
            const Complex ca(sampling::Random() - Real(0.5), sampling::Random() - Real(0.5));
            const Complex cb(sampling::Random() - Real(0.5), sampling::Random() - Real(0.5));
            const CVec3 e0 = sum(scaled(ea, ca), scaled(eb, cb));

            // The field through the two reflections: s and p components with their own signs (p = k x s on each side).
            const Vector3 pIn2 = glm::cross(kIn2, s2), pOut2 = glm::cross(wo2, s2);
            const CVec3 e1 = sum(scaled(s2, a2.rs * dotC(e0, s2)), scaled(pOut2, a2.rp * dotC(e0, pIn2)));
            const Vector3 kIn1 = -wi1;
            const Vector3 pIn1 = glm::cross(kIn1, s1), pOut1 = glm::cross(wo1, s1);
            const CVec3 e2 = sum(scaled(s1, a1.rs * dotC(e1, s1)), scaled(pOut1, a1.rp * dotC(e1, pIn1)));

            const Vector3 k0 = wo1;
            const Vector3 xCam = perpendicularTo(k0);
            Real sCam[4], sIn[4];
            stokesOf(e2, xCam, k0, sCam);
            stokesOf(e0, s2, kIn2, sIn);   // in the basis of the incoming light at the last event

            // The chain, as updatePolarization() builds it.
            pol::Mueller M1 = pol::interfaceMueller(a1.rs, a1.rp);
            pol::Mueller T = pol::Mueller::identity();
            Vector3 xPrev = xCam;
            T = T * pol::basisChange(s1, xPrev, wo1) * M1;
            xPrev = s1;
            // The second event is the metal: its matrix uses one amplitude pair for every component here.
            pol::Mueller M2 = pol::interfaceMueller(a2.rs, a2.rp);
            T = T * pol::basisChange(s2, xPrev, wo2) * M2;
            const Real gain = Real(0.25) * (std::norm(a1.rs) + std::norm(a1.rp)) * (std::norm(a2.rs) + std::norm(a2.rp));   // M00 of both
            Real got[4];
            for (int i = 0; i < 4; ++i) {
                Real v = 0;
                for (int j = 0; j < 4; ++j) v += T.e[i][j][0] * sIn[j];
                got[i] = v * gain;
            }
            for (int i = 0; i < 4; ++i) worst = std::max(worst, std::abs(got[i] - sCam[i]) / std::max(Real(1e-9), sIn[0]));
        }
        report("two reflections (glass, metal) with random directions: Mueller chain = field propagation, " + std::to_string(trials) + " cases",
            worst < 1e-9, "max difference " + std::to_string(worst) + " of the input intensity");
    }

    // (3) The path tracer: a camera ray reflected by smooth glass under a uniform sky. Only the reflection counts (depth 2: the refracted
    //     ray is not followed), so the Stokes vector is F (1, B, 0, 0) in the basis (s, p) of the reflection: horizontal polarization
    //     for a horizontal surface seen from the side, for the basis x = along s, and the same turned for another reference direction.
    {
        auto world = std::make_shared<Aggregate>();
        const Real n = 1.5;
        world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0, std::make_shared<Dielectric>(n, 0.0)));
        auto scene = std::make_shared<Scene>(world);
        auto sky = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
        scene->finalize("spatial");

        for (const Real deg : { Real(56.3099), Real(40) }) {
            const Real ti = deg * kPi / 180;
            const Real tt = std::asin(std::sin(ti) / n);
            const Real rs = std::pow(std::sin(ti - tt) / std::sin(ti + tt), Real(2));
            const Real rp = std::pow(std::tan(ti - tt) / std::tan(ti + tt), Real(2));
            const Vector3 d(std::sin(ti), -std::cos(ti), 0);   // plane of incidence: x-y, so s is the z axis
            const Ray ray(Point3(0, 0, 0) - d * Real(10), d);
            const Vector3 k = -d;
            for (const Real turnDeg : { Real(0), Real(30), Real(90) }) {
                const Real turn = turnDeg * kPi / 180;
                const Vector3 sAxis(0, 0, 1), u = glm::cross(k, sAxis);
                const Vector3 xRef = std::cos(turn) * sAxis + std::sin(turn) * u;   // the reference direction of the image, turned by `turn` from s
                const Vector3 yRef = glm::cross(k, xRef);
                const Real alpha = std::atan2(glm::dot(sAxis, yRef), glm::dot(sAxis, xRef));   // angle of s in the image basis
                const Real e0 = Real(0.5) * (rs + rp), e1 = Real(0.5) * (rs - rp) * std::cos(2 * alpha), e2 = Real(0.5) * (rs - rp) * std::sin(2 * alpha);

                PathIntegrator integrator(nullptr, 2, 1);
                const int N = 80000;
                Real sum[4] = { 0, 0, 0, 0 }, sum2[4] = { 0, 0, 0, 0 };
                for (int i = 0; i < N; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
                    SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
                    SampledWavelengths lambda;
#endif
                    PathStokes st;
                    st.xRef = xRef;
                    integrator.Li(ray, *scene, lambda, nullptr, &st);
                    for (int c = 0; c < 4; ++c) {
                        const Real v = mean3(SpectrumToRGB(st.S[c], lambda));
                        sum[c] += v;
                        sum2[c] += v * v;
                    }
                }
                const Real expected[4] = { e0, e1, e2, 0 };
                bool ok = true;
                std::string detail;
                for (int c = 0; c < 4; ++c) {
                    const Real m = sum[c] / N, se = std::sqrt(std::max(Real(0), sum2[c] / N - m * m) / N);
                    ok = ok && std::abs(m - expected[c]) <= 4 * se + Real(1e-3) * e0;
                    detail += (c ? ", " : "") + std::to_string(m) + " (" + std::to_string(expected[c]) + ")";
                }
                report("glass at " + std::to_string(int(std::lround(deg))) + " deg, image axis turned " + std::to_string(int(turnDeg)) + " deg from s: S0..S3 (expected)", ok, detail);
            }
        }
    }

    std::cout << "[Polarization] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
