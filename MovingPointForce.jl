#|||||||||||||#
#||--Intro--||#
#|||||||||||||#

#=
This file solves for the dynamic deformation of an infinite Euler-Bernoulli
beam on a piecewise constant viscoelastic foundation subject to a moving load.
The model is solved using the method of undetermined coefficients.
=#

include("Functions-MovingPointForce.jl")


#|||||||||||||||||#
#||--Variables--||#
#|||||||||||||||||#

#||--Model parameters--||#
#Things like the beam stiffness, etc.
EI = 1 #Beam stiffness
m = 1 #Mass per unit length of beam
C0 = 1 #Foundation damping before tz
C1 = 2 #Foundation damping after tz
k0 = 1 #Foundation stiffness before tz
k1 = 2 #Foundation stiffness after tz

xtz = 1 #Position of transition zone

P = 1 #Force conveyed by point load
v = 1 #Speed of point force

parameters = [EI, m, C0, C1, k0, k1, xtz, P, v]


#||--Solution parameters--||#
#Things like the time and space values that for which I want an evaluation.
xLeft = -10
xRight = 10
xNum = 100
xVals = LinRange(xLeft,xRight,xNum)

tMin = 0
tMax = 5
tNum = 101
tVals = LinRange(tMin,tMax,tNum)



#||||||||||||||||||||||||||||#
#||--Finding the solution--||#
#||||||||||||||||||||||||||||#

#Calculate the deformation of the beam for given times and x values
deformations = CalcDynamicDeformation(xVals,tVals, parameters)


##

#||||||||||||||||#
#||--Graphing--||#
#||||||||||||||||#

#TODO Implement graphing of solution.
