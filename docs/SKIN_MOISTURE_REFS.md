# Optics of moisturizing (lotion and emulsion): record of literature values (2026-11-09)

**Status:** A record of a literature survey, based on **search summaries (and in part WebFetch)**. **No paper was read in full** (text extraction from PDFs was also tried, but only fragments came out). Check every number against the original before using it. Confidence marks: **○** = confirmed with a source in a search summary. **△** = only one source, or the source is vague. **×** = not found, or my own estimate. The model that uses these values is in `docs/SKIN_MOISTURE.md`.

Purpose: to back up the coefficients used in the simulation of lotion → emulsion → make-up (foundation) (the coefficients of `docs/SKIN.md` are typical values from memory or values set by eye).

## 1. Refractive index

| Material | n | Confidence | Source / note |
|---|---|---|---|
| Air / water | 1.0000 / 1.3300 | ○ | J. Soc. Cosmet. Chem. 9(1) (an old table; library.scconline.org/v009n01/26) |
| Glycerin | 1.47 | ○ | Same. Mixed into the water phase, it narrows the difference from the oil phase |
| White mineral oil | 1.475 | △ | Same. Varies with the grade |
| Squalane | 1.4555 (20 °C) | △ | A comparison table of cosmetic oils (table in OCL 2018). The original table was not checked |
| Dimethicone (linear) | 1.403 (20 °C). 1.387 to 1.403 | △ | A table in a review, patents |
| Petrolatum, sebum | — | × | Not found. For sebum, the only option is to use the lipid layer of tears as a substitute (n ≈ 1.33 + a small wavelength term; not a value for skin) |
| Dry stratum corneum | about 1.55 (Na D line). 1.54 | ○ | Cited by several papers (Saiko 2023 and others) |
| Stratum corneum (in vivo, OCT, one subject) | 1.52. Living epidermis 1.34, dermis 1.41 | △ | The same source notes that in vitro values differ |
| Stratum corneum (another OCT estimate, 1300 nm) | 1.51 ± 0.02 | △ | |
| Hydrated stratum corneum | Only "decreases". **No value was found** | × | A direct measurement of the change in refractive index was not found |

## 2. Dimensions of the stratum corneum and the surface

| Quantity | Value | Confidence | Source / note |
|---|---|---|---|
| Thickness of the stratum corneum | 10 to 20 µm (excluding flat, thick places such as the sole). Some papers give 10 to 40 µm | ○ | Several |
| Same ("apparent thickness" from the gradient of water concentration) | cheek about 17 µm, forearm about 23 µm | △ | Egawa et al. 2007 (confocal Raman). Search summary |
| Thickness when hydrated | dry about 9 µm → 40 minutes under occlusion with water about 19 µm → dried about 8 µm | △ | Visible-light OCT (back; one study). **"Extreme hydration by soaking in water"**, not a lotion |
| Change of thickness with the water content of the stratum corneum | One description gives 52 → 16 µm at low water content | △ | The source is fragmentary |
| Corneocytes | diameter about 25 to 50 µm, thickness 0.2 to 0.5 µm, 15 to 25 layers (some papers say 15 to 20) | ○ | Several (Holbrook & Odland's electron microscopy, and others) |
| Lipid film on the skin surface | about 0.5 µm or almost none on the limbs; on the sebum-rich face, locally over 4 µm. Amount 1 µg/cm² (leg) to 189 ± 43 µg/cm² (face) | ○ | An electron-microscopy study (search summary; author and year not checked) |
| Skin furrows (forearm, one 30-year-old woman) | maximum depth Rtm 63 to 81 µm, mean spacing of furrows Sm 385 to 443 µm | △ | **Values for one person.** Ra and Rz of the face and cheek were not found. There is a report that applying a moisturizing cream for 14 days (twice a day) reduces Rz and others (light-transmission method, elderly women) |

**Relation to my own `SkinRelief`:** the typical values in `SkinRelief.hpp` (cell 0.5 to 0.8 mm, furrow depth 30 to 50 µm, pore depth 60 to 100 µm) are of roughly the same order as the forearm values above (spacing 0.4 mm, maximum depth 60 to 80 µm).

## 3. Amount applied and formulation

| Quantity | Value | Confidence | Source / note |
|---|---|---|---|
| Usual amount of moisturizing cream | about 2 mg/cm² (Simion et al. 2006, as cited in a clinical-trial protocol) | △ | A description in a protocol. 2 mg/cm² is **a film of about 20 µm** at a density of 1 g/cm³ (**my own calculation**; depends on the density of the product) |
| Sunscreen in actual use | 0.5 to 1.3 mg/cm² (tests use 2 mg/cm²) | ○ | Not the amount of a lotion |
| Amount of lotion used, film thickness | — | × | **Not found** |
| Formulation of a lotion | water about 97%, glycerin 1% (an alcohol-free toner, an example from a raw-material maker). There is also an example with glycerin 5% | △ | **Sites of raw-material makers and enthusiasts.** Not peer-reviewed papers |
| When a lotion of only water and glycerin dries | the water evaporates and the glycerin remains (qualitative description) | △ | A forum description. **The drying time was not found** |
| Water loss through the skin (TEWL) | healthy adult forearm about 6 to 12 g/m²/h, the cheek and forehead the same order. Increases with age | ○ | A meta-analysis such as Kottner et al. 2013. **This is evaporation through the skin, not evaporation of an applied water film** |
| Speed of evaporation of an applied water film | × | × | **Not found.** Since 1 g/m² = 1 µm of water, if the evaporation were that of a free water surface (roughly several hundred g/m²/h; **my own estimate, unverified**), a 20 µm film would vanish within minutes. Absorption into the stratum corneum also matters |

### 3.1 Layering (an emulsion after a lotion): occlusion, evaporation, persistence of hydration (additional survey, 2026-11-09)

**Only search summaries and page summaries were obtained.** The full texts were not read. Where author and year could not be obtained, the source is given by volume, issue and page.

| Quantity | Value | Confidence | Source / note |
|---|---|---|---|
| Time-course comparison of lotion → emulsion layering with a lotion alone | — | × | **Not found** (in Japanese or English). There is a report that applying a moisturizing product after bathing gives a high water content of the stratum corneum from 1 to 60 minutes (a single product) |
| Occlusion (fraction by which evaporation is reduced): petrolatum | 98 to 99% | △ | Widely cited, but only secondary sources were found (trade journals, sales sites). The original was not checked |
| Same: mineral oil, lanolin, silicone | 20 to 30% | △ | The same secondary sources. For mineral oil other sources disagree with about 11% and about 45% |
| Same: relation to the amount of petrolatum | examined at 0.25 to 2 mg/cm², "very large even for a small amount". At 3 mg/cm² it decreases for 2 hours but is still high after 3 hours | △ | J. Soc. Cosmet. Chem. 30(6) (<https://library.scconline.org/v030n06/17>). The numbers in the figure could not be read. 2 subjects |
| Same: o/w emulsion | depends on the formulation. One with no occlusion for 90 minutes, and one (with a high oil-phase ratio) from 7.7 ± 2.1 of bare skin to below 5.6 g/m²h after 30 minutes (about 30% less). Occlusion appears "after most of the water has evaporated". Petrolatum occludes more strongly than mineral oil | △ | J. Soc. Cosmet. Chem. 39(2) (<https://library.scconline.org/v039n02/63>). Amount applied unknown |
| Same: one example of a w/o emulsion | 6.5 ± 2.5% | △ | Search summary only |
| Evaporation of water from o/w and w/o | faster from o/w | △ | J. Soc. Cosmet. Chem. 23(6) (search summary) |
| Evaporation of applied water (forearm, distilled water 4 to 6 mg/cm²) | the measured evaporation is high right after application and returns to the original value **within 30 minutes** | △ | J. Soc. Cosmet. Chem. 43(3) (<https://library.scconline.org/v043n03/31>). 4 to 6 mg/cm² is 40 to 60 µm (**my own conversion**). The order agrees with the default "10 µm, time constant 120 seconds" |
| Same: an o/w lotion with a humectant (propylene glycol) | still high after 40 minutes, and 5 to 6 times that of water after 2 hours | △ | Same. **The effect of the humectant holding water** (not in the present model) |
| Evaporation of water from an emulsion (on a gelatin film) | all in 10 to 15 minutes | △ | Search summary. Not skin |
| Persistence of hydration (single application) | high after 2 hours for many products. Still high after 24 hours for those rich in glycerin and those with a lamellar structure | △ | Electrical measurement (it also measures the humectant remaining on the stratum corneum). The first measurement is at 1 to 2 hours, so **there is no fine time course for the first 30 minutes** |
| Water loss after removing occlusion | decreases within minutes. One report fits it with a sum of two exponentials. The only time constant found is 3.1 minutes (another condition: a patent) | △ | A model with a concentration-dependent diffusion coefficient (Sci. Rep. 2022; abstract only) |
| Remaining film and evaporation, occlusivity and water retention | the effect of reducing evaporation depends on the thickness and nature of the remaining film. Moisturizing effect = occlusivity + water retention | △ | J. Soc. Cosmet. Technologists Japan (日本化粧品技術者会誌) 13(2) 1979, 21(2) 1987 (abstracts only) |

**How it is used in the model:** the default occlusion of 0.3 is the order of mineral oil and silicone, and of the o/w emulsions that do work. To agree with "acts after the water has evaporated", only oil that has merged into a film is counted for occlusion. The decay of hydration (default 30 minutes) lies between "minutes" and "hours", with nothing to decide it.

## 4. Optical constants of skin (inputs of the layered model)

| Quantity | Value | Confidence | Source / note |
|---|---|---|---|
| Reduced scattering coefficient μs′ of skin (whole tissue) | μs′(λ) = a (λ / 500 nm)^(−b). **a = 46.0 ± 13.7 (n = 8), b = 1.421 ± 0.517** (check the unit in the original paper; it should be cm⁻¹). In another form, a′ = 48.0, f_Rayleigh = 0.409, b_Mie = 0.702 | ○ | Table 2 of Jacques 2013 (*Phys. Med. Biol.* 58, R37). Obtained from the OMLC page (WebFetch). **An average over the whole skin, not a value for the stratum corneum alone** |
| Anisotropy g of each layer | 0.7 to 0.9 (strong forward scattering) | △ | A review by van Gemert, Jacques et al. 1989 (search summary) |
| Scattering coefficient of the stratum corneum | × | × | The numbers in the table of a review could not be read in the search summary (check the table of Bashkatov, Genina, Tuchin 2011 in the original) |
| Effect of hydration (in vivo) | the spectral reflectance and the scattering coefficient **increased with dehydration and decreased with hydration**. The skin became darker and redder, and the perceived translucency rose | △ | J. Soc. Cosmet. Chem. 44(5) (library.scconline.org/v044n05/20, search summary). **Authors and conditions not checked** |
| Why dry skin looks white | the scales diffusely reflect white light. The brightness is proportional to the thickness of the scales | △ | Image analysis of dry skin in J. Soc. Cosmet. Chem. (search summary) |
| Light confinement in the stratum corneum | because the stratum corneum has a higher refractive index than the tissue below, the light intensity is 11 to 17% higher (analytic model). The effect is notable when the stratum corneum is 150 µm or thicker (sole, callus) | △ | Saiko 2023 (arXiv:2309.13146, **abstract only**). It implies the effect is small for the thin stratum corneum of the face (10 to 20 µm) |
| Amount of sebum and reflection | the more sebum, the lower the directional reflectance at 400 to 600 nm, and the higher at 590 nm and above | △ | A 2023 in vivo study (search summary) |

## 5. Make-up (reference; the next stage)

- Huang et al. (*Comput. Animat. Virtual Worlds* 24(3-4), 2013): Kubelka-Munk + screen-space SSS. Foundation, blush and lipstick were measured. Pearl and glitter excluded.
- Doi et al. (2005, 2009): a two-layer model that estimates the reflectance after make-up from the reflectance of bare skin and the absorption and scattering coefficients of the foundation.
- Lanza et al. (EGSR 2024, *Comput. Graph. Forum* 43(4)): a medium with two kinds of scatterers (spheres and plates) per layer. Multilayer BSDF. Validated by measurement. The code is **non-commercial only**.
- Okamoto et al. (*J. Biomed. Opt.* 18(6), 2013): Monte Carlo of single scattering in a powder layer.

## 6. What is known and what is not

**What is known**
1. **The stratum corneum is thin** (10 to 20 µm; 15 to 25 layers of cells). It is 1 to 2 orders of magnitude thinner than the scattering length of skin (hundreds of µm to mm with the coefficients of `docs/SKIN.md`). So if "hydration changes the internal scattering of the stratum corneum" matters for the overall reflection, it is more natural to express it as **a "thin scattering layer just below the surface" and the refractive-index difference at the interfaces**, not as the scattering distance of SSS that includes the dermis.
2. Values that were obtained: refractive indices of water, glycerin, mineral oil, squalane and dimethicone, the refractive index of the stratum corneum (dry), the thickness of the stratum corneum, the dimensions of corneocytes, the dimensions of skin furrows (one person), the usual amount used (cream).
3. Order of magnitude of the water-film thickness: the 20 µm range (converted from the amount of cream). The skin furrows (maximum depth 60 to 80 µm on the forearm) **are not completely filled** (a water film alone should collect at the bottom of the furrows by capillary force; inference).

**What is not known (not found; needs measurement or further survey)**
- The change in the refractive index and the scattering coefficient of hydrated stratum corneum.
- The amount of lotion used (mg/cm²) and the time until it dries.
- The refractive index of sebum and petrolatum. The size of oil droplets in an emulsion, and the refractive index of an emulsion.
- The scattering and absorption coefficients of the stratum corneum alone. The dimensions of the skin furrows of the face (cheek).
- The optics of scales (scattering strength versus thickness).

## 7. References (found by search; bibliographic data not cross-checked)

- Jacques SL. *Optical properties of biological tissues: a review.* Phys Med Biol 58, R37, 2013. (<https://omlc.org/news/dec14/Jacques_PMB2013/index.html>)
- Bashkatov AN, Genina EA, Tuchin VV. *Optical properties of skin, subcutaneous, and muscle tissues: a review.* J Innov Opt Health Sci 4(1), 9-38, 2011. (<https://www.worldscientific.com/doi/10.1142/S1793545811001319>)
- Saiko G. *Light confinement in stratum corneum.* arXiv:2309.13146, 2023. (abstract only)
- Egawa M et al. *In vivo estimation of stratum corneum thickness from water concentration profiles obtained with Raman spectroscopy.* Acta Derm Venereol, 2007.
- Kottner J et al. *Transepidermal water loss in young and aged healthy humans: a systematic review and meta-analysis.* Arch Dermatol Res, 2013. (<https://pubmed.ncbi.nlm.nih.gov/23341028/>)
- J. Soc. Cosmet. Chem. 44(5) (appearance of moisturizing; <https://library.scconline.org/v044n05/20>), 9(1) (table of refractive indices; <https://library.scconline.org/v009n01/26>)
- Emulsified film on the skin surface and gloss (TU Delft; <https://repository.tudelft.nl/file/File_2d303454-6fec-4305-b05e-f229fc760d38>), hydration of the stratum corneum by visible-light OCT (White Rose; <https://eprints.whiterose.ac.uk/id/eprint/205782/1/Visible-light%20optical%20coherence%20tomography%20platform%20for%20the%20characterization%20of%20the%20skin%20barrier.pdf>).
- Huang et al. 2013, Doi et al., Lanza et al. 2024, Okamoto et al. 2013: section 5.
