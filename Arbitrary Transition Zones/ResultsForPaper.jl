using WAV #To play a sound when computation finished

include("Functions - Inversion schemes.jl")
include("Functions - MPF_ArbitraryTZs.jl")

#|||||||||||||||||#
#||--Variables--||#
#|||||||||||||||||#

#||--Model parameters--||#
#Things like the beam stiffness, etc.
#For constant parameters I use the default values from
  #"Steady-state dynamic response of a Bernoulli-Euler beam on a
  #viscoelastic foundation subject to a platoon of moving loads"
#For changing variables(i.e. k, C, v), I use the ranges from the same paper.

const EI = 2.3*10^3 #Beam stiffness
const m = 48 #Mass per unit length of beam
const P = -10^4 #Force conveyed by point load



const xp = 0 #Starting position of point force




##

#||||||||||||||||||||||||||||||||||||#
#||--Comparing to travelling wave--||#
#||||||||||||||||||||||||||||||||||||#


k = 6.9*10^7
C = 10^7

v_cr = CriticalVelocity(EI, m, k) #~128 for this set of parameters

v_sub = 100
v_sup = 160


leftTZ = TransitionZone(2, k, 2*k, C, 2*C)
xtz_list = [leftTZ]


parameters_sub = [EI, m, xp, xtz_list, P, v_sub]
parameters_sup = [EI, m, xp, xtz_list, P, v_sup]



#||--Solution parameters--||#
#The time and space values for which I want an evaluation.
xLeft = 0
xRight = 100
xNum = 100
xVals = LinRange(xLeft,xRight,xNum)

tMin = 0
tMax = xRight/v_sup
tNum = 100
tVals = LinRange(tMin,tMax,tNum)


#inversionMethod = WeeksMethodImplementation
inversionMethod = directQuadratureMethodImplementation
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters_sup, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals,parameters_sup, k, C)

SteadyStateDeformations = [SteadyDeformation1]



include("Graphing - MPF_ArbitraryTZs.jl")
GifWithFeatures(deformations, xVals, tVals, parameters_sup, SteadyStateDeformations)


#Play sound when computation finished.
y, fs = wavread(raw"C:/Windows/Media/Ring01.wav")
wavplay(y, fs)

println("ENDED")
