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
 *
 * Several products, one after the other (a lotion, then an emulsion: docs/SKIN_MOISTURE.md 5.6). Each is a `Layer` with the time it was applied.
 *  - The film is one mixture (a lotion and an oil-in-water emulsion both have water as the continuous phase): water, humectant and oil add up.
 *  - H(t) above is the solution of a store of water on the surface that empties into the stratum corneum, which then loses it:
 *        dS/dt = -S / rise,   dH/dt = k S - H / decay.
 *    Every layer adds to the store when it is applied (so that alone it would reach its `hydrationPeak`), and H never exceeds 1.
 *  - Occlusion: the oil that has merged into a film slows the loss, dH/dt = k S - H (1 - occ) / decay, occ = occlusion (1 - exp(-merged oil / occlusionOil)).
 *    Oil that is still droplets in water does not count (an oil-in-water emulsion becomes occlusive after its water has gone), and the evaporation of the
 *    water of the film is not slowed. `occlusion` and `occlusionOil` are not measured.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

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
            // An emulsion (milky lotion) also carries oil, as droplets in the water. The oil stays (as droplets, then as a film), the droplets scatter white light
            // while they are in the water, and they merge into a clear film with time (docs/SKIN_MOISTURE.md 5.5).
            Real oil = Real(0);                   ///< volume fraction of the applied liquid that is oil (0 = a plain lotion)
            Real oilIor = Real(1.45);
            Real dropletDiameter = Real(1e-6);    ///< diameter of the oil droplets (m)
            Real coalesceTau = Real(600);         ///< time constant of the merging of the droplets (s)
            Real dropletG = Real(0.9);            ///< asymmetry of the scattering by a droplet (large droplets scatter forwards)
            // The merged oil film slows the loss of water from the stratum corneum (docs/SKIN_MOISTURE.md 5.6).
            Real occlusion = Real(0.3);           ///< fraction of the loss that a thick oil film stops, 0..1 (0 = oil does not hold water back)
            Real occlusionOil = Real(1e-6);       ///< thickness of merged oil at which 63 % of that is reached (m)
        };

        /// One product put on the skin at the time `at`. The other parameters (evaporation, uptake and loss, the indices of water and humectant) belong to the skin
        /// and the room and are shared: they come from the `Params` given with the layers.
        struct Layer {
            Real at = 0;                          ///< when it was applied (s)
            Real dose = Real(10e-6), glycerin = Real(0.05), hydrationPeak = Real(0.8);
            Real oil = Real(0), oilIor = Real(1.45), dropletDiameter = Real(1e-6), coalesceTau = Real(600), dropletG = Real(0.9);
        };

        struct State {
            Real thickness = 0;     ///< mean thickness of the film (m); 0 = nothing left
            Real water = 0;         ///< ... of which water
            Real residue = 0;       ///< ... of which humectant
            Real oil = 0;           ///< ... of which oil
            Real ior = Real(1);     ///< refractive index of the film (volume mean)
            Real roughness = 0;     ///< roughness of its surface
            Real hydration = 0;     ///< 0..1
            Real emulsionTau = 0;   ///< transport optical depth (1 - g) tau of the droplets across the film of mean thickness; 0 = none
            Real occlusion = 0;     ///< fraction of the water loss of the stratum corneum that the merged oil stops now, 0..1
        };

        /// Efficiency of a droplet of relative index m and size parameter x for the extinction (anomalous diffraction, van de Hulst: good for m near 1 and x above a few).
        static Real dropletEfficiency(Real m, Real x) {
            const Real rho = Real(2) * x * (m - Real(1));
            if (std::abs(rho) < Real(0.05)) return rho * rho / Real(2) - rho * rho * rho * rho / Real(36);
            return Real(2) - (Real(4) / rho) * std::sin(rho) + (Real(4) / (rho * rho)) * (Real(1) - std::cos(rho));
        }

        /// Valid when everything is positive where it must be, the humectant fraction is below 1 and the hydration decays slower than it rises. Empty string = valid.
        static const char* validate(const Params& p) {
            if (!(p.dose >= 0)) return "doseUm must not be negative";
            if (!(p.glycerin >= 0 && p.oil >= 0 && p.glycerin + p.oil < 1)) return "glycerin and oil must not be negative and must add up to less than 1";
            if (!(p.oilIor >= 1 && p.dropletDiameter > 0 && p.coalesceTau > 0 && p.dropletG >= 0 && p.dropletG < 1)) return "oilIor must be at least 1, dropletUm and coalesceTauS above 0, dropletG in [0, 1)";
            if (!(p.waterIor >= 1 && p.humectantIor >= 1)) return "the refractive indices must be at least 1";
            if (!(p.waterTau > 0 && p.hydrationRise > 0 && p.hydrationDecay > 0)) return "the time constants must be above 0";
            if (!(p.hydrationDecay > p.hydrationRise)) return "hydrationDecayS must be longer than hydrationRiseS";
            if (!(p.hydrationPeak >= 0 && p.hydrationPeak <= 1)) return "hydrationPeak must be in [0, 1]";
            if (!(p.occlusion >= 0 && p.occlusion <= 1 && p.occlusionOil > 0)) return "occlusion must be in [0, 1] and occlusionOilUm above 0";
            return "";
        }

        /// The product of `p` as a layer applied at `at`.
        static Layer layerOf(const Params& p, Real at = Real(0)) {
            Layer l;
            l.at = at; l.dose = p.dose; l.glycerin = p.glycerin; l.hydrationPeak = p.hydrationPeak;
            l.oil = p.oil; l.oilIor = p.oilIor; l.dropletDiameter = p.dropletDiameter; l.coalesceTau = p.coalesceTau; l.dropletG = p.dropletG;
            return l;
        }

        /// The shared parameters with the product of one layer (to validate a layer with `validate`).
        static Params withLayer(Params shared, const Layer& l) {
            shared.dose = l.dose; shared.glycerin = l.glycerin; shared.hydrationPeak = l.hydrationPeak;
            shared.oil = l.oil; shared.oilIor = l.oilIor; shared.dropletDiameter = l.dropletDiameter; shared.coalesceTau = l.coalesceTau; shared.dropletG = l.dropletG;
            return shared;
        }

        /// Time after an application at which the hydration it gives is largest (s).
        static Real hydrationPeakTime(const Params& p) {
            return std::log(p.hydrationDecay / p.hydrationRise) / (Real(1) / p.hydrationRise - Real(1) / p.hydrationDecay);
        }

        /// The hydration of one application with a peak of 1, `t` seconds after it (no occlusion).
        static Real hydrationShape(const Params& p, Real t) {
            if (t <= 0) return Real(0);
            const Real tp = hydrationPeakTime(p);
            const Real norm = std::exp(-tp / p.hydrationDecay) - std::exp(-tp / p.hydrationRise);
            return (std::exp(-t / p.hydrationDecay) - std::exp(-t / p.hydrationRise)) / norm;
        }

        static Real hydration(const Params& p, Real t) { return p.hydrationPeak * hydrationShape(p, t); }

        /// Thickness of the oil that has merged into a film at the time `t` (m): the droplets of every layer merge with its own time constant.
        static Real mergedOil(const std::vector<Layer>& layers, Real t) {
            Real film = 0;
            for (const Layer& l : layers)
                if (t > l.at && l.oil > 0) film += l.dose * l.oil * (Real(1) - std::exp(-(t - l.at) / l.coalesceTau));
            return film;
        }

        /// Fraction of the water loss of the stratum corneum that the merged oil stops at the time `t`.
        static Real occlusionAt(const Params& shared, const std::vector<Layer>& layers, Real t) {
            if (!(shared.occlusion > 0)) return Real(0);
            return shared.occlusion * (Real(1) - std::exp(-mergedOil(layers, t) / shared.occlusionOil));
        }

        /// The hydration at the time `t`, by integrating dH/dt = (uptake from the store of every layer) - H (1 - occlusion) / decay in steps of at most 0.5 s. H stays at or below 1.
        /// Within a step the water a layer delivers is exact (also in the step it is applied in) and is put at the middle of the step; the loss rate is the one at the middle.
        static Real hydrationAt(const Params& shared, const std::vector<Layer>& layers, Real t) {
            if (t <= 0 || layers.empty()) return Real(0);
            const Real tp = hydrationPeakTime(shared);
            const Real norm = std::exp(-tp / shared.hydrationDecay) - std::exp(-tp / shared.hydrationRise);
            const Real gain = (Real(1) / shared.hydrationRise - Real(1) / shared.hydrationDecay) / norm;   // k S0 of a layer with a peak of 1
            const int steps = std::clamp(int(std::ceil(double(t) / 0.5)), 1, 400000);
            const Real dt = t / Real(steps);
            Real h = 0;
            for (int i = 0; i < steps; ++i) {
                const Real t0 = Real(i) * dt, t1 = t0 + dt, tm = t0 + Real(0.5) * dt;
                Real delivered = 0;
                for (const Layer& l : layers) {
                    if (t1 <= l.at) continue;
                    const Real from = std::max(t0, l.at);
                    delivered += l.hydrationPeak * gain * shared.hydrationRise * (std::exp(-(from - l.at) / shared.hydrationRise) - std::exp(-(t1 - l.at) / shared.hydrationRise));
                }
                const Real lambda = (Real(1) - occlusionAt(shared, layers, tm)) / shared.hydrationDecay;
                h = h * std::exp(-lambda * dt) + delivered * std::exp(-lambda * Real(0.5) * dt);
                h = std::min(h, Real(1));
            }
            return h;
        }

        /// The state at the time `t` when the layers were applied at their times `at` (a layer with at > t is not on the skin yet).
        static State at(const Params& shared, const std::vector<Layer>& layers, Real t) {
            State s;
            t = std::max(t, Real(0));
            Real oilIndexVolume = 0, peaks = 0;
            for (const Layer& l : layers) {
                if (t < l.at) continue;
                s.water += l.dose * (Real(1) - l.glycerin - l.oil) * std::exp(-(t - l.at) / shared.waterTau);
                s.residue += l.dose * l.glycerin;
                s.oil += l.dose * l.oil;
                oilIndexVolume += l.oilIor * l.dose * l.oil;
                peaks += l.hydrationPeak;
            }
            s.thickness = s.water + s.residue + s.oil;
            if (s.thickness > Real(1e-12)) {
                s.ior = (shared.waterIor * s.water + shared.humectantIor * s.residue + oilIndexVolume) / s.thickness;
                s.roughness = (shared.waterRoughness * s.water + shared.residueRoughness * (s.residue + s.oil)) / s.thickness;
            }
            else s.thickness = 0;

            // Without occlusion the model is linear: the layers add up (and stay below 1 when their peaks do). Otherwise integrate.
            if ((!(shared.occlusion > 0) || s.oil <= Real(1e-12)) && peaks <= Real(1)) {
                for (const Layer& l : layers)
                    if (t > l.at) s.hydration += l.hydrationPeak * hydrationShape(shared, t - l.at);
            }
            else s.hydration = hydrationAt(shared, layers, t);
            s.occlusion = occlusionAt(shared, layers, t);

            if (s.oil > Real(1e-12)) {   // the droplets scatter in proportion to the oil per area, with an efficiency that falls as the water around them dries (their index gets closer to the medium's)
                const Real mediumVolume = s.water + s.residue;
                const Real nMedium = mediumVolume > Real(1e-12) ? (shared.waterIor * s.water + shared.humectantIor * s.residue) / mediumVolume : shared.humectantIor;
                for (const Layer& l : layers) {
                    if (t < l.at || !(l.oil > 0)) continue;
                    const Real x = Real(3.14159265358979) * l.dropletDiameter * nMedium / Real(550e-9);
                    const Real q = dropletEfficiency(l.oilIor / nMedium, x);
                    s.emulsionTau += Real(1.5) * q * l.dose * l.oil / l.dropletDiameter * (Real(1) - l.dropletG) * std::exp(-(t - l.at) / l.coalesceTau);
                }
            }
            return s;
        }

        /// The state `t` seconds after the one product of `p` was applied.
        static State at(const Params& p, Real t) { return at(p, std::vector<Layer>{ layerOf(p) }, t); }
    };

} // namespace rayt
