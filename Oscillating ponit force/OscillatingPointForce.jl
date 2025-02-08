#|||||||||||||#
#||--Intro--||#
#|||||||||||||#

#=
This file solves for the dynamic deformation of an infinite Euler-Bernoulli
beam on a piecewise constant viscoelastic foundation subject to a moving load.
The model is solved using the method of undetermined coefficients.
=#

include("Functions -  OscillatingPointForce.jl")


#|||||||||||||||||#
#||--Variables--||#
#|||||||||||||||||#

#||--Model parameters--||#
#Things like the beam stiffness, etc.
EI = 2.3*10^3 #Beam stiffness
m = 48 #Mass per unit length of beam
C0 = 10^7 #Foundation damping before tz
C1 = 2*C0 #Foundation damping after tz
k0 = 6.9*10^7 #Foundation stiffness before tz
k1 = 2*k0 #Foundation stiffness after tz

xtz = 1 #Position of transition zone
xp = 0 #Starting position of point force

P = 10^3 #Force conveyed by point load #BUG force has negative effect for some reason. Missing negative probably.
v = 0 #Speed of point force

parameters = [EI, m, C0, C1, k0, k1, xp, xtz, P, v]

#||--Solution parameters--||#
#Things like the time and space values that for which I want an evaluation.
xLeft = -1
xRight = 4
xNum = 100
xVals = LinRange(xLeft,xRight,xNum)

tMin = 0
tMax = 6
tNum = 100
tVals = LinRange(tMin,tMax,tNum)


#||||||||||||||||||||||||||||#
#||--Finding the solution--||#
#||||||||||||||||||||||||||||#

#Calculate the deformation of the beam for given times and x values
deformations = CalcDynamicDeformation(xVals,tVals, parameters)
println("Finished calculating - Beginning graphing:")
##

#||||||||||||||||||||||||||||||||||||||||#
#||--Solution with no transition zone--||#
#||||||||||||||||||||||||||||||||||||||||#
parameters = [EI, m, C0, C0, k0, k0, xp, xtz, P, v]
homogeneousDeformations = CalcDynamicDeformation(xVals,tVals, parameters)


#||||||||||||||||#
#||--Graphing--||#
#||||||||||||||||#

include("Graphing - OscillatingPointForce.jl")



#|||||||||||||||||||||||||||#
#||--Graphing difference--||#
#|||||||||||||||||||||||||||#


differenceDeformation =zeros(size(deformations))

for i in eachindex(deformations)
  differenceDeformation[i] = deformations[i]-homogeneousDeformations[i]
end
