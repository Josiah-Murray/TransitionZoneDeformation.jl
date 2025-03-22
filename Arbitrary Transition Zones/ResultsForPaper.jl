using WAV #To play a sound when computation finished
using Plots
include("Functions - Inversion schemes.jl")
include("Functions - MPF_ArbitraryTZs.jl")
include("Graphing - MPF_ArbitraryTZs.jl")


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


const xLeft = 0
const xRight = 1
const xNum = 100
const xVals = LinRange(xLeft,xRight,xNum)

const tMin = 0
const tMax = 1
const tNum = 100
const tVals = LinRange(tMin,tMax,tNum)

colour1 = RGB(0.19215685, 0.27843137, 0.3333333)
colour2 = RGB(0.14901960784313725, 0.6274509803921569, 0.8549019607843137)

colour3 = RGB(0, 0.27450980392156865, 0.4980392156862745)
colour4 = RGB(0.6470588235294118, 0.8, 0.5098039215686274)

redLineColour = RGB(0.5,0.2,0.15)





##

#||||||||||||||||||||||||||||||||||||#
#||--Comparing to travelling wave--||#
#||||||||||||||||||||||||||||||||||||#

#=
#An initial verification of the method

k = 6.9*10^7
C = 10^7
v = 1



leftTZ = TransitionZone(xRight+1, k, k, C, C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]

laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals,parameters, k, C)

SteadyStateDeformations = [SteadyDeformation1]






#GifWithFeatures(deformations, xVals, tVals, parameters, SteadyStateDeformations, false)

GraphTimeEvolution(deformations, xVals, tVals, 1:10:100, parameters, SteadyStateDeformations, false, colour4, colour3)
=#

##

#||||||||||||||||||||||||||||||||||||#
#||--Comparing to travelling wave--||#
#||||||||||||||||||||||||||||||||||||#

#An initial verification of the method

k = 6.9*10^7
C = 10^7
v = 1



leftTZ = TransitionZone(xRight+1, k, k, C, C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]

laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals,parameters, k, C)

SteadyStateDeformations = [SteadyDeformation1]






#GifWithFeatures(deformations, xVals, tVals, parameters, SteadyStateDeformations, false)

GraphTimeEvolution(deformations, xVals, tVals, 1:10:100, parameters, SteadyStateDeformations, false, colour4, colour3)








#Play sound when computation finished.
y, fs = wavread(raw"C:/Windows/Media/Ring01.wav")
wavplay(y, fs)

println("ENDED")
