using WAV #To play a sound when computation finished
using Plots

figurePath = "C:/Users/joemu/Documents/PhD/Julia/Moving Point Force/Figures/"

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


const xLeft = -1
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

#An initial verification of the method

tVals10 = LinRange(tMin,tMax,10)

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
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, figurePath, graphName)

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
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, figurePath, graphName)


#||--N=1024--||#
println("Weeks 1024")
function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xVals,tVals10, parameters, inversionMethod)
graphName = "ComparingToTravellingWaveWeeks1024"
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00055,0.00015], colour2, figurePath, graphName)



#||||||||||||||||||||||||||||||||||#
#||--Including a transition zone--||#
#||||||||||||||||||||||||||||||||||#

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
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

MultiPlotTimeEvolution(deformations, xVals, tVals, [1,10,20,30,40,50,60,70,80,90,100], parameters, SteadyStateDeformations, false, colour4, colour3)

#||--And back down--||#


leftTZ = TransitionZone(0.4, k1, k2, C1, C2)
rightTZ = AddConsistentTZ(0.6, leftTZ, k1,C1)
xtz_list = [leftTZ, rightTZ]


parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

MultiPlotTimeEvolution(deformations, xVals, tVals, [1,10,20,30,40,50,60,70,80,90,100], parameters, SteadyStateDeformations, false, colour4, colour3)

#||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#||--Gradual increase with lots of transition zones--||#
#||||||||||||||||||||||||||||||||||||||||||||||||||||||#


#||--With Weeks--||#

numTransitionZones = 4
firstTZLocation = 0.4
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
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

MultiPlotTimeEvolution(deformations, xVals, tVals, [1,10,20,30,40,50,60,70,80,90,100], parameters, SteadyStateDeformations, false, colour4, colour3)

#||--With direct quadrature--||#

laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)


deformations = @time CalcDynamicDeformation(xVals,tVals, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals,tVals,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

MultiPlotTimeEvolution(deformations, xVals, tVals, [1,10,20,30,40,50,60,70,80,90,100], parameters, SteadyStateDeformations, false, colour4, colour3)



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


leftTZ = TransitionZone(0.4, k1, k2, C1, C2)
tz2 = AddConsistentTZ(0.6, leftTZ, k1, C1)
tz3 = AddConsistentTZ(1.4,tz2,k2,C2)
tz4 = AddConsistentTZ(1.6, tz3, k1, C1)
xtz_list = [leftTZ,tz2,tz3,tz4]

#||--With Weeks--||#




parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
deformations = @time CalcDynamicDeformation(xValsLong,tValsLong, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xValsLong,tValsLong,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xValsLong,tValsLong,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

#TODO Fix wrong point force location
plot_list = MultiPlotTimeEvolution_ReturnPlotList(deformations, xValsLong, tValsLong, [1,10,20,30,40,50,60,70,80,90,100], parameters, SteadyStateDeformations, [-0.001,0.001], colour4, colour3)
plot(plot_list[1:10]..., layout = grid(5,2), size=(1200,1000))
#||--With direct quadrature--||#

laplaceParameters = [50,100]

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)


deformations = @time CalcDynamicDeformation(xValsLong,tValsLong, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xValsLong,tValsLong,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xValsLong,tValsLong,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

MultiPlotTimeEvolution(deformations, xValsLong, tValsLong, [1,10,20,30,40,50,60,70,80,90,100], parameters, SteadyStateDeformations, [-0.001,0.001], colour4, colour3)


#Play sound when computation finished.
y, fs = wavread(raw"C:/Windows/Media/Ring01.wav")
wavplay(y, fs)

println("ENDED")
