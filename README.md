# TransitionZoneDeformation.jl

TransitionZoneDeformation.jl provides a suite of tools for computing the deformation of railway tracks at transition zones. The track is modelled as an Euler-Bernoulli beam on a piecewise-homogeneous viscoelastic foundation under point forces. The code is only equiped to solve for the effect of one point load at a time, however the system is linear, and so the effect of multiple point forces (or distributed load) can be found by summing the results of multiple computations. More details about the model can be found in my research [1,2]

Below is a sample figure showing a moving point load moving over successive transtion zones, gradually increasing the foundation stiffness and damping: (The red dashed curves indicate the steady-state solutions for the left-most and right-most regions, and the black curve indicates the time-dependent deformation)

<img width="800" height="800" alt="GradiatedTransitionZones_10" src="https://github.com/user-attachments/assets/e43f59ad-0b00-4532-8eb7-36d8b0a15741" />

## Installation
```julia
pkg> add https://github.com/Josiah-Murray/SpatialLaplaceInversion.git
```

## References
[1] Murray, J., Meylan, M. H., Ngo, T., Thamwattana, N., & Indraratna, B. (2025). Analytical Solution for Railway Transition Zones With Abrupt Changes in Elastic Stiffness. Engineering Reports, 7(5). https://doi.org/10.1002/eng2.70106

[2] Murray, J., Pethiyagoda, R., Meylan, M., & Thamwattana, N. (2026). Laplace transform method for the deformation of an Euler–Bernoulli beam on a piecewise homogeneous foundation with moving point load. Applied Mathematics in Science and Engineering, 34(1). https://doi.org/10.1080/27690911.2026.2657615
