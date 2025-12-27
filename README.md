 A lightweight 3D renderer built from scratch.

Aurisidus
A high-performance, physically-based C++ rendering engine dedicated to the rigorous simulation of light transport and material appearance.

🌌 Overview
Aurisidus is a next-generation renderer developed with a core focus on the mathematical rigor of metallic reflections and spectral accuracy. By moving beyond simple RGB-based calculations and embracing spectral power distributions, Aurisidus aims to capture the true optical essence of metals and complex materials.
While currently specialized in advanced material research, Aurisidus is architected to evolve into a versatile, general-purpose production renderer.
🛠 Key Research Themes
• Rigorous Metallic Reflection Models: Implementing advanced microfacet distributions (e.g., GGX, GTR) and wave-optics-inspired models to simulate the complex interplay of light on metallic surfaces.
• Full Spectral Rendering: Shifting from traditional tri-stimulus (RGB) pipelines to a full spectral decomposition to eliminate color shifts and accurately render phenomena like thin-film interference and metamerism.
• High-Performance Architecture: Optimized in C++ with a focus on SIMD utilization, efficient memory layouts, and modern acceleration structures (BVH) to ensure high-speed convergence.
• Algorithmic Innovation: Developing custom importance sampling techniques to reduce noise in complex specular-standard-specular light paths.
🚀 Future Roadmap
• [x] Core Spectral Integrator
• [x] Advanced Metallic BRDFs
• [ ] Real-time Denoising Integration
• [ ] Support for Volumetric Path Tracing
• [ ] Universal Scene Description (USD) Support for General-Purpose Workflows
