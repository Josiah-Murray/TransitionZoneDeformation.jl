sqrt(-1) #So that I don't accidentally run all the code at once.
using WAV #To play a sound when computation finished
using Plots
using JLD

figurePath = "C:/Users/joemu/Documents/PhD/Julia/Moving Point Force/Figures/"

include("Functions - Inversion schemes.jl")
include("Functions - MPF_ArbitraryTZs.jl")
include("Functions BigFloat - MPF_ArbitraryTZs.jl")
include("Graphing - MPF_ArbitraryTZs.jl")

#MARK: Variables
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


const xLeft = -1
const xRight = 1
const xNum = 100
const xVals = LinRange(xLeft,xRight,xNum)

const tMin = 0
const tMax = 1
const tNum = 100
const tVals = LinRange(tMin,tMax,tNum)

const tVals10 = LinRange(tMin,tMax,10)

colour1 = RGB(0.19215685, 0.27843137, 0.3333333)
colour2 = RGB(0.14901960784313725, 0.6274509803921569, 0.8549019607843137)

colour3 = RGB(0, 0.27450980392156865, 0.4980392156862745)
colour4 = RGB(0.6470588235294118, 0.8, 0.5098039215686274)

redLineColour = RGB(0.5,0.2,0.15)





##
#MARK: Verification
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
deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k, C)

SteadyStateDeformations = [SteadyDeformation1]
graphName = "ComparingToTravellingWaveDQ"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)

##

#|||||||||||||||||||||||||||||||||||||||||||||#
#||--Weeks method for travelling wave case--||#
#|||||||||||||||||||||||||||||||||||||||||||||#



k = 6.9*10^7
C = 10^7
v = 1



leftTZ = TransitionZone(xRight+1, k, k, C, C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]

SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k, C)
SteadyStateDeformations = [SteadyDeformation1]

#||--N=256--||#
println("Weeks 256")
function N(x)
  return 256
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
graphName = "ComparingToTravellingWaveWeeks256"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)


#||--N=1024--||#
println("Weeks 1024")
function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
graphName = "ComparingToTravellingWaveWeeks1024"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)

#MARK: Transition zones

#||||||||||||||||||||||||||||||||||#
#||--Including a transition zone--||#
#||||||||||||||||||||||||||||||||||#

#TODO: Shade region of higher stiffness

k1 = 6.9*10^7
k2 = 2*k1
C1 = 10^7
C2 = 2*C1
v = 1

#||--Up--||#

leftTZ = TransitionZone(0.4, k1, k2, C1, C2)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]
graphName = "TransitionZonex0p4Weeks1024"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)


#||--And back down--||#


leftTZ = TransitionZone(0.3, k1, k2, C1, C2)
rightTZ = AddConsistentTZ(0.7, leftTZ, k1,C1)
xtz_list = [leftTZ, rightTZ]


parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

graphName = "TransitionZonex0p40p6Weeks1024"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)

#||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#||--Gradual increase with lots of transition zones--||#
#||||||||||||||||||||||||||||||||||||||||||||||||||||||#
k1 = 6.9*10^7
k2 = 2*k1
C1 = 10^7
C2 = 2*C1
v = 1

#||--With Weeks--||#

numTransitionZones = 4
firstTZLocation = 0.3
leftTZ = TransitionZone(firstTZLocation, k1, k1 + (k2-k1)/numTransitionZones, C1, C1 + (C2-C1)/numTransitionZones)
xtz_list = [leftTZ]
for i in 2:numTransitionZones
  global xtz_list = [xtz_list AddConsistentTZ(firstTZLocation + (xRight-firstTZLocation)*(i-1)/(numTransitionZones), xtz_list[end],k1+i*(k2-k1)/numTransitionZones, C1 + i*(C2-C1)/numTransitionZones)]
end



parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

graphName = "GradualTransitionZoneWeeks1024"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)

#||--With direct quadrature--||#

laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)


deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

graphName = "GradualTransitionZoneDQ"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)


#|||||||||||||||||||#
#||--Two bridges--||#
#|||||||||||||||||||#

const xLeftLong = -1
const xRightLong = 2
const xNumLong = 200
const xValsLong = LinRange(xLeftLong,xRightLong,xNumLong)

const tMinLong = 0
const tMaxLong = 2
const tNumLong = 100
const tValsLong = LinRange(tMinLong,tMaxLong,tNumLong)
const tVals10Long = LinRange(tMinLong,tMaxLong, 10)


leftTZ = TransitionZone(0.4, k1, k2, C1, C2)
tz2 = AddConsistentTZ(0.8, leftTZ, k1, C1)
tz3 = AddConsistentTZ(1.2,tz2,k2,C2)
tz4 = AddConsistentTZ(1.6, tz3, k1, C1)
xtz_list = [leftTZ,tz2,tz3,tz4]

#||--With Weeks--||#




parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xValsLong,tVals10Long, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xValsLong,tVals10Long,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xValsLong,tVals10Long,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

graphName = "TwoBridges1024"
gr()
Graph10Times(deformations,xValsLong,tVals10Long,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)

#||--With direct quadrature--||#

laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)


deformations = @time CalcDynamicDeformation(xValsLong,tVals10Long, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xValsLong,tVals10Long,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xValsLong,tVals10Long,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

graphName = "TwoBridgesDQ"
gr()
Graph10Times(deformations,xValsLong,tVals10Long,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, colour1, figurePath, graphName)


#MARK: Numerical stuff

#||||||||||||||||||||||||||||#
#||--Numerical assessment--||#
#||||||||||||||||||||||||||||#

#||--More realistic case--||#
const xLeft_Full = 0
const xRight_Full = 30
const xNum_Full = 100
const xVals_Full = LinRange(xLeft_Full,xRight_Full,xNum_Full)

const tMin_Full = 0
const tMax_Full = 1
const tNum_Full = 100
const tVals_Full = LinRange(tMin_Full,tMax_Full,tNum_Full)

tVals10_Full = LinRange(tMin_Full,tMax_Full,10)

k = 6.9*10^7
C = 10^7
v = 30



leftTZ = TransitionZone(5, k, 2*k, C, 2*C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]


function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
#inversionMethod = WeeksMethodImplementation

#laplaceParameters = [50,100]

#inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)

deformations = @time CalcDynamicDeformation(xVals_Full,tVals10_Full, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals_Full,tVals10_Full,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, 2*k, 2*C)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

graphName = "LongIntervalWeeks1024"
gr()
Graph10Times(deformations,xVals_Full,tVals10_Full,parameters,SteadyStateDeformations, [-0.00005,0.00005], colour2, colour1, figurePath, graphName)
#[-0.000035,0.000015]



#MARK: BigFloat
#||--For comparison with BigFloat stuff--||#

k = 6.9*10^7
C = 10^7
v = 30



leftTZ = TransitionZone(5, k, 2*k, C, 2*C)
rightTZ = AddConsistentTZ(10,leftTZ,2*k,2*C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]



laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)

deformations = @time CalcDynamicDeformation(xVals_Full,tVals10_Full, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals_Full,tVals10_Full,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, 2*k, 2*C)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

graphName = "LongIntervalDQ_50_100"
gr()
Graph10Times(deformations,xVals_Full,tVals10_Full,parameters,SteadyStateDeformations, [-0.00005,0.00005], colour2, colour1, figurePath, graphName)
#[-0.000035,0.000015]

#----Actual BigFloat------#

k = 6.9*10^7
C = 10^7
v = 30



leftTZ = TransitionZoneBF(5, k, 2*k, C, 2*C)
rightTZ = AddConsistentTZBF(10,leftTZ,2*k,2*C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]


laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)

deformations = @time CalcDynamicDeformationBF(xVals_Full,tVals10_Full, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, 2*k, 2*C)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]
graphName = "LongIntervalWeeks_BF_DQ_50_100"
gr()
Graph10Times(deformations,xVals_Full,tVals10_Full,parameters,SteadyStateDeformations, [-0.000035,0.000035], colour2, colour1, figurePath, graphName)
#[-0.000035,0.000015]

#||--BigFloat Weeks--||#

k = 6.9*10^7
C = 10^7
v = 30



leftTZ = TransitionZoneBF(5, k, 2*k, C, 2*C)
rightTZ = AddConsistentTZBF(10,leftTZ,2*k,2*C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 526
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)

deformations = @time CalcDynamicDeformationBF(xVals_Full,tVals10_Full, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, 2*k, 2*C)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]
graphName = "LongIntervalWeeks1024_BF"
gr()
Graph10Times(deformations,xVals_Full,tVals10_Full,parameters,SteadyStateDeformations, [-0.000035,0.000035], colour2, colour1, figurePath, graphName)
#[-0.000035,0.000015]


#||--BigFloat Finer parameters--||#


k = 6.9*10^7
C = 10^7
v = 30



leftTZ = TransitionZoneBF(5, k, 2*k, C, 2*C)
rightTZ = AddConsistentTZBF(10,leftTZ,2*k,2*C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]


laplaceParameters = [200,4000]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)

#----Calculation pre-saved to "LongIntervalWeeks_BF_DQ_200_4000.jld"----#
deformations = load("LongIntervalWeeks_BF_DQ_200_4000.jld")["data"]
#deformations = @time CalcDynamicDeformationBF(xVals_Full,tVals10_Full, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, 2*k, 2*C)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]
graphName = "LongIntervalWeeks_BF_DQ_200_4000"
gr()
Graph10Times(deformations,xVals_Full,tVals10_Full,parameters,SteadyStateDeformations, [-0.000035,0.00002], colour2, colour1, figurePath, graphName)
#[-0.000035,0.000015]

#MARK: Play sound

#Play sound when computation finished.
y, fs = wavread(raw"C:/Windows/Media/Ring01.wav")
wavplay(y, fs)

println("ENDED")
