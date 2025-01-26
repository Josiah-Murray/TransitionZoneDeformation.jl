#|||||||||||||#
#||--Intro--||#
#|||||||||||||#

#=
This file solves for the dynamic deformation of an infinite Euler-Bernoulli
beam on a piecewise constant viscoelastic foundation subject to a moving load.
The model is solved using the method of undetermined coefficients.
=#
using Plots
plotlyjs()
include("Functions-MovingPointForce.jl")


#|||||||||||||||||#
#||--Variables--||#
#|||||||||||||||||#

#||--Model parameters--||#
#Things like the beam stiffness, etc.
EI = 2.3*10^3 #Beam stiffness
m = 48 #Mass per unit length of beam
C0 = 10^7 #Foundation damping before tz
C1 = C0 #Foundation damping after tz
k0 = 10^7 #Foundation stiffness before tz
k1 = 2*k0 #Foundation stiffness after tz

xtz = 1 #Position of transition zone
xp = 0 #Starting position of point force

P = -10^4 #Force conveyed by point load
v = 1 #Speed of point force

parameters = [EI, m, C0, C1, k0, k1, xp, xtz, P, v]


#||--Solution parameters--||#
#Things like the time and space values that for which I want an evaluation.
xLeft = -1
xRight = 2
xNum = 100
xVals = LinRange(xLeft,xRight,xNum)

tMin = 0
tMax = 2
tNum = 100
tVals = LinRange(tMin,tMax,tNum)



#||||||||||||||||||||||||||||#
#||--Finding the solution--||#
#||||||||||||||||||||||||||||#

#Calculate the deformation of the beam for given times and x values
deformations = CalcDynamicDeformation(xVals,tVals, parameters)
##

tempDeformations = copy(deformations)
for i in eachindex(tempDeformations)
  tempDeformations[i] = max(-0.001, tempDeformations[i])
  tempDeformations[i] = min(0.001, tempDeformations[i])
end



myPlt = plot(xVals, tVals,  deformations, st=:surface)
display(myPlt)
println("Finished")
##

#||||||||||||||||#
#||--Graphing--||#
#||||||||||||||||#

#TODO Implement graphing of solution.
