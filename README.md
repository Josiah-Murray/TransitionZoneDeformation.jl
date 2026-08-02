# Euler-Bernoulli beams and railway transition zones

This repository contains the code associated with the first half of my PhD. It is primarily concerned with solving for the deformation of an Euler-Bernoulli beam resting on a viscoelastic foundation with abrupt changes in the foundation properties (i.e. transition zones).
The code (as well as this README.md) is an active work in progress as I continue working on my PhD.
Notably, the code depends on a custom module `NILaplace`, in the file `LaplaceInversionImplementations.jl`. I haven't made this module public, however, if you're interested, feel free to contact me.

Below is a sample figure showing a moving point load moving over successive transtion zones, gradually increasing the foundation stiffness and damping: (The red dashed curves indicate the steady-state solutions for the left-most and right-most regions, and the black curve indicates the time-dependent deformation)

<img width="800" height="800" alt="GradiatedTransitionZones_10" src="https://github.com/user-attachments/assets/e43f59ad-0b00-4532-8eb7-36d8b0a15741" />
