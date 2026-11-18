# Polarization (minimal version)

The polarization state of light is followed as a Stokes vector and the action of an interface as a Mueller matrix, and the **polarization of light reflected and refracted by glass, water, metal and thin films** is turned into an image. Images with a polarizing filter in front of the camera can also be produced.

- Code: `include/Core/Polarization.hpp` (Stokes, Mueller, Fresnel amplitudes), `updatePolarization` in `Integrator.cpp`, `Dielectric::scatterMueller`
- Self-check: `GoLD_rayt.exe --scene polarization` (closed-form Fresnel, the Mueller chain checked against a direct 3D electric-field calculation, the Stokes vector of a ray reflected by glass, and the same under a polarized sky)
- Usage: `--stokes` (`<out>_S0..S3.pfm`, degree of polarization `_dolp.png`, angle of polarization `_aolp.png`), `--polarizer 0,45,90` (images through a polarizing filter, `<out>_pol0.png` ...), `--sky-pol P[,EL,AZ]` (a polarized clear sky, below)

## Approach

The path is traced from the camera, and a Mueller matrix is multiplied in at each event of a material that knows about polarization (reflection and refraction at smooth glass and water, reflection from metal, reflection and transmission of thin films). When the path hits a diffuse surface, a rough dielectric, or scattering in a medium, the light beyond that point is treated as unpolarized. Light that samples a light source directly (NEE) is also counted through the Mueller matrix of that vertex. Unpolarized light has Stokes vector (L, 0, 0, 0), so the polarization seen at the camera is "the first column of the chain × the scalar radiance". The existing scalar path tracing is unchanged; polarization is just accumulated alongside it. Images without `--stokes` do not change.

The plain image keeps the scalar radiance, so it leaves out what polarization does to the intensity: a second polarizing event that meets light polarized by the first, or polarized light from the sky. The intensity that includes it is `<out>_S0.pfm`; use S0 (or the polarizer images) when that matters.

## A polarized sky (`--sky-pol`)

Light from a clear sky is partially polarized by Rayleigh scattering, most strongly at 90° from the sun. A reflection then depends on how the sky's polarization lines up with the plane of incidence, and that changes its intensity, which the naked eye sees: glass at Brewster's angle reflects light polarized along s and none along p. `--sky-pol P` gives the environment light (`EnvLight::setSkyPolarization`) the polarization of single Rayleigh scattering: the light from a direction at angle γ from the sun has degree P sin²γ / (1 + cos²γ) and its electric field perpendicular to the plane of the sun and that direction. P is the degree at 90° from the sun (a clear sky has about 0.6 to 0.8; haze and clouds lower it). The sun is the brightest spot of the map, or `--sky-pol P,EL,AZ` gives its elevation and azimuth (degrees, azimuth from +z toward +x). The radiance of the map is not changed. Only directions above the horizon (+y) are polarized, and the degree is the same at every wavelength. The option implies `--stokes`.

A light can now say that it is polarized (`Light::polarization`: a normalized Stokes vector in a basis of its own). The integrator then sends that Stokes vector through the whole chain, turned into the basis of the last specular event, instead of only the chain's first column; it does so for the light that a path finds (`addLightFrom`) and for the light it samples (NEE), so the MIS weights of the two strategies stay consistent. Lights that do not say so, and every light without `--sky-pol`, are counted as before: with the same `--seed`, the outputs of a render without `--sky-pol` are byte for byte those of the renderer before this change (checked on a measured-material scene, `--scene sea` with polarizers and `--scene gold`).

Self-check (part 4 of `--scene polarization`): the sky seen straight has its own Stokes vector; glass at 56° and 40° and a narrow-lobe metal at 70° and 45° under a sky polarized along s, along p and at 45° give the Stokes vectors that the electric fields give when they are reflected with the Fresnel amplitudes (the polarized part and the unpolarized rest as incoherent fields), without the integrator's Mueller matrices or basis rotations. Glass at Brewster's angle reflects 1.7 times as much under a sky polarized along s and 0.3 times along p (P = 0.7); the metal, which both sampling strategies reach, agrees within 0.1 %.

Example: three spheres on an outdoor map (`bryanston_park_sunrise_2k.hdr`, sun 5° high), the camera at 90° from the sun, `--sky-pol 0.75`, same seed with and without. The upper half of a black glass sphere, which reflects the strongly polarized sky overhead, falls to 0.34 of its brightness at its darkest (1st percentile; 27 % of the sphere changes by more than 10 %); iron (Fresnel from n, k) falls to 0.65; aluminium hardly changes (within 6 %, except a thin rim at grazing incidence). Well-reflecting metals reflect s and p almost alike; dielectrics, and metals that reflect less, do not.

## Examples

`--scene sea --stokes --polarizer 0,90` (800×450, 64 spp, 46 seconds in the RGB development build). With the horizontal polarizing filter (0°), the reflection from the sea surface (horizontally polarized) remains and the stones under water are barely visible. With vertical (90°), the reflection weakens and the stones become clearly visible. This is the well-known effect of using a polarizing filter to remove the glare from a water surface.

In the soap bubble (`scenes/soap_bubble.json --stokes --polarizer 0,90`, 900×900, 96 spp, 35 seconds in the RGB development build), the reflections at the top and bottom of the bubble, where the plane of incidence is vertical, are bright through the horizontal filter and dark through the vertical one (the interference colors of the film stay as they are).

The gold sphere (`--scene gold --stokes`) has a small degree of polarization (almost 0 at the front) that grows toward the edge, where light arrives obliquely. The direction of polarization is along the edge of the sphere (tangential).

## Conventions

- A Stokes vector is expressed with the direction of propagation k and a reference direction x (y = k × x). S1 > 0 is linear polarization along x, and S2 > 0 is 45° between x and y. The image reference is the right direction of the camera.
- The basis of an interface is (s, p, k). s is perpendicular to the plane of incidence, p = k × s. The Fresnel amplitudes are those of Born & Wolf, with r_p = −r_s at normal incidence.
- Each specular event multiplies "the Mueller matrix ÷ its first element" onto the chain from the right. The first element (the unpolarized reflectance or transmittance) is already carried by the weight of the existing path. The rotation of the basis between events is `basisChange`.

## Limits

- The supported materials are dielectrics (smooth glass and water), metals (`RoughConductor`, at any roughness), thin films and measured materials. A metal uses the complex Fresnel amplitudes at the microfacet normal (the half vector), and the multiple-scattering term is mixed in as unpolarized. A thin film (`thinFilm`) builds its Mueller matrix from the complex amplitudes (reflection and transmission) that include the multiple reflections in the film. The spectral build computes this per wavelength of the path, and the RGB build integrates with the wavelength weights of each channel, into a table over (cos θ, thickness) (built the first time it is used). A measured material (`"type": "measured"`, `MeasuredHpBRDF`) returns its measured Mueller matrix (divided by the unpolarized element), rotated into the (s, p, k) basis. Coatings and rough dielectrics are treated as depolarizing.
- Light sources are unpolarized; the environment light is unpolarized unless `--sky-pol` makes it a Rayleigh sky. The sky model is a clear sky of single scattering: over an HDR map with clouds, buildings or trees above the horizon, those are polarized too. Birefringence and polarization inside media are not included.
- The light added by the photon map (caustics), `--lvc` and manifold NEE is unpolarized (also under `--sky-pol`).
