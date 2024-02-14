# NetworkPDE.jl: efficient solution and optimization of PDEs on networks via symbolic programming
Note that this package is a work in progress! It relies on a fork of ModelingToolkit.jl (while we wait for a necessary PR to be integrated, https://github.com/cmhyett/ModelingToolkit.jl) to work. 

We utilize symbolic infrastructure in Julia to streamline solving PDEs on networks. Of particular note is the state-of-the-art integration with reverse-mode automatic sensitivity analysis, allowing for efficient calculation of gradients even in high-dimensional parameterizations, and enabling a differentiable-programming approach to PDE-constrained optimization problems on networks.

Networked systems are ubiquitous, and their efficient characterization and control are keys to many of today's challenges: e.g., planning and executing the clean energy transition, reliable delivery of clean water, optimal operation of traffic and public transport, among others. These systems are described by PDEs on the edges, augmented with appropriate (e.g., Kirchoff) boundary conditions coupling the network edges. This boundary-coupling of PDE systems usually requires so-called "book-keeping" in the solver. We utilize symbolic infrastructure in Julia to symbolically pre-compute coupling conditions and eliminate book-keeping, allowing for more efficient memory accesses in solving the forward problem, and integrating natively with Julia's differential equations and AD ecosystem. Of particular note is the full operability of reverse-mode AD, enabling sensitivity calculations that do not depend on the dimensionality of the parameterization - facilitating a new and efficient approach of PDE-constrained optimization problems on networks.