# Polarization (minimal version)

The polarization state of light is followed as a Stokes vector and the action of an interface as a Mueller matrix, and the **polarization of light reflected and refracted by glass, water, metal and thin films** is turned into an image. Images with a polarizing filter in front of the camera can also be produced.

- Code: `include/Core/Polarization.hpp` (Stokes, Mueller, Fresnel amplitudes), `updatePolarization` in `Integrator.cpp`, `Dielectric::scatterMueller`
- Self-check: `GoLD_rayt.exe --scene polarization` (closed-form Fresnel, the Mueller chain checked against a direct 3D electric-field calculation, the Stokes vector of a ray reflected by glass)
- Usage: `--stokes` (`<out>_S0..S3.pfm`, degree of polarization `_dolp.png`, angle of polarization `_aolp.png`), `--polarizer 0,45,90` (images through a polarizing filter, `<out>_pol0.png` ...)

## Approach

The path is traced from the camera, and a Mueller matrix is multiplied in at each event of a material that knows about polarization (reflection and refraction at smooth glass and water, reflection from metal, reflection and transmission of thin films). When the path hits a diffuse surface, a rough dielectric, or scattering in a medium, the light beyond that point is treated as unpolarized. Light that samples a light source directly (NEE) is also counted through the Mueller matrix of that vertex. Unpolarized light has Stokes vector (L, 0, 0, 0), so the polarization seen at the camera is "the first column of the chain × the scalar radiance". The existing scalar path tracing is unchanged; polarization is just accumulated alongside it. Images without `--stokes` do not change.

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
- Light sources and environment light are unpolarized. Polarized light sources, birefringence, and polarization inside media are not included.
- The light added by the photon map (caustics), `--lvc` and manifold NEE is unpolarized.
