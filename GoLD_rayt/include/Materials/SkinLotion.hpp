/**
 * @file Materials/SkinLotion.hpp
 * @brief The time course of a lotion on skin: how thick the film is, what it is made of, and how hydrated the stratum corneum is (docs/SKIN_MOISTURE.md, 2.3).
 *
 * A water-based lotion is mostly water with a little humectant (glycerin). Of the liquid that was applied (a film of mean thickness `dose`)
 *  - the water part, dose (1 - glycerin), evaporates: water(t) = dose (1 - glycerin) exp(-t / waterTau);
 *  - the humectant stays: residue = dose glycerin.
 * The film is both: thickness = water + residue, the refractive index and the roughness are their volume-weighted means (water 1.33 and smooth,
 * the humectant 1.47 and a little rougher). The stratum corneum takes up water, then loses it: the hydration rises and falls,
 *     H(t) = peak (exp(-t / decay) - exp(-t / rise)) / max,   H(0) = 0, one maximum `peak` at t = ln(decay / rise) / (1 / rise - 1 / decay).
 * None of the time constants is measured (the film can be weighed as it dries, the hydration measured with a corneometer); they are parameters to vary.
 * The state at one time is what the film (Materials/SkinFilm.hpp) and the scatter (Materials/SkinScatter.hpp) of the skin material are built with.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"

namespace rayt {

    class SkinLotion {
    public:
        struct Params {
            Real dose = Real(10e-6);              ///< mean thickness of the liquid as applied (m)
            Real glycerin = Real(0.05);           ///< volume fraction of the liquid that does not evaporate (humectant)
            Real waterIor = Real(1.33), humectantIor = Real(1.47);
            Real waterTau = Real(120);            ///< time constant of the evaporation of the water film (s)
            Real waterRoughness = Real(0.02), residueRoughness = Real(0.08);   ///< perceptual roughness of the film surface
            Real hydrationPeak = Real(0.8);       ///< largest hydration of the stratum corneum, 0..1
            Real hydrationRise = Real(60);        ///< time constant of the uptake (s)
            Real hydrationDecay = Real(1800);     ///< time constant of the loss (s); longer than the rise
        };

        struct State {
            Real thickness = 0;     ///< mean thickness of the film (m); 0 = nothing left
            Real water = 0;         ///< ... of which water
            Real residue = 0;       ///< ... of which humectant
            Real ior = Real(1);     ///< refractive index of the film (volume mean)
            Real roughness = 0;     ///< roughness of its surface
            Real hydration = 0;     ///< 0..1
        };

        /// Valid when everything is positive where it must be, the humectant fraction is below 1 and the hydration decays slower than it rises. Empty string = valid.
        static const char* validate(const Params& p) {
            if (!(p.dose >= 0)) return "doseUm must not be negative";
            if (!(p.glycerin >= 0 && p.glycerin < 1)) return "glycerin must be in [0, 1)";
            if (!(p.waterIor >= 1 && p.humectantIor >= 1)) return "the refractive indices must be at least 1";
            if (!(p.waterTau > 0 && p.hydrationRise > 0 && p.hydrationDecay > 0)) return "the time constants must be above 0";
            if (!(p.hydrationDecay > p.hydrationRise)) return "hydrationDecayS must be longer than hydrationRiseS";
            if (!(p.hydrationPeak >= 0 && p.hydrationPeak <= 1)) return "hydrationPeak must be in [0, 1]";
            return "";
        }

        /// Time at which the hydration is largest (s).
        static Real hydrationPeakTime(const Params& p) {
            return std::log(p.hydrationDecay / p.hydrationRise) / (Real(1) / p.hydrationRise - Real(1) / p.hydrationDecay);
        }

        static Real hydration(const Params& p, Real t) {
            if (t <= 0) return Real(0);
            const Real tp = hydrationPeakTime(p);
            const Real norm = std::exp(-tp / p.hydrationDecay) - std::exp(-tp / p.hydrationRise);
            return p.hydrationPeak * (std::exp(-t / p.hydrationDecay) - std::exp(-t / p.hydrationRise)) / norm;
        }

        /// The state `t` seconds after the lotion was applied.
        static State at(const Params& p, Real t) {
            State s;
            t = std::max(t, Real(0));
            s.water = p.dose * (Real(1) - p.glycerin) * std::exp(-t / p.waterTau);
            s.residue = p.dose * p.glycerin;
            s.thickness = s.water + s.residue;
            if (s.thickness > Real(1e-12)) {
                s.ior = (p.waterIor * s.water + p.humectantIor * s.residue) / s.thickness;
                s.roughness = (p.waterRoughness * s.water + p.residueRoughness * s.residue) / s.thickness;
            }
            else s.thickness = 0;
            s.hydration = hydration(p, t);
            return s;
        }
    };

} // namespace rayt
