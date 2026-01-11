#pragma once
#include "Core/Precision.hpp"
#include "Core/Types.hpp"
#include "Color/ComplexIorTable.hpp"

namespace rayt::color {

    /**
     * @brief Extracts RGB-like eta/k by sampling the IOR table at 3 representative wavelengths.
     *
     * This is a pragmatic approximation for RGB mode. It will be replaced by proper
     * spectral-to-RGB conversion (CMF integration) later.
     */
    inline void approxIorToRgb_3Wavelengths(
        const ComplexIorTable& table,
        Vector3& outEtaRgb,
        Vector3& outKRgb,
        Real lambdaB_nm = Real(450),
        Real lambdaG_nm = Real(550),
        Real lambdaR_nm = Real(650))
    {
        auto [etaB, kB] = table.evaluateNK(lambdaB_nm);
        auto [etaG, kG] = table.evaluateNK(lambdaG_nm);
        auto [etaR, kR] = table.evaluateNK(lambdaR_nm);

        // Note: We pack as (R,G,B).
        outEtaRgb = Vector3(etaR, etaG, etaB);
        outKRgb = Vector3(kR, kG, kB);
    }

} // namespace rayt::color