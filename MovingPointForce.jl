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
k1 = 4*k0 #Foundation stiffness after tz

xtz = 2 #Position of transition zone
xp = 0 #Starting position of point force

P = 10^4 #Force conveyed by point load
v = 1 #Speed of point force

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
##

tempDeformations = copy(deformations)
for i in eachindex(tempDeformations)
  tempDeformations[i] = max(-0.001, tempDeformations[i])
  tempDeformations[i] = min(0.001, tempDeformations[i])
end




##

#||||||||||||||||#
#||--Graphing--||#
#||||||||||||||||#

#||||----Surface----||||#

myPlt = plot(xVals, tVals,  deformations, st=:surface)
display(myPlt)
println("Finished")


#TODO Add an animation.

@gif for ti in 1:tNum
  p = plot(xVals, deformations[ti,:], ylimits = (-0.01,0.001))
end
