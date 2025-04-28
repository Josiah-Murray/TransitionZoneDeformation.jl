using WAV #To play a sound when computation finished
using Plots
using JLD

println(pwd())

DataPath = "Paper/Data"

include(joinpath(functionFolder, "Functions - Inversion schemes.jl"))
include(joinpath(functionFolder, "Functions - MPF_ArbitraryTZs.jl"))
include(joinpath(functionFolder, "Graphing - MPF_ArbitraryTZs.jl"))


#||||||||||||||||||||||/
#MARK: Constants
#||||||||||||||||||||||\

SaveFlag = true
LoadFlag = false

#||--Model constants--||#

const EI = 2.3*10^3 #Beam stiffness
const m = 48 #Mass per unit length of beam
const P = -10^4 #Force conveyed by point load

const xp = 0 #Starting position of point force


#||--Space time setups--||#

#--> Short track

const v_short = 1

const xLeft_short = -1
const xRight_short = 1
const xNum_short = 100
const xVals_short = LinRange(xLeft_short,xRight_short,xNum_short)

const tMin_short = 0
const tMax_short = 1
const tNum_short = 100
const tVals_short = LinRange(tMin_short,tMax_short,tNum_short)
const tVals10_short = LinRange(tMin_short,tMax_short,10)


#--> 90m track
const v_90m = 30

const xLeft_90m = -10
const xRight_90m = 90
const xNum_90m = 100
const xVals_90m = LinRange(xLeft_90m,xRight_90m,xNum_90m)

const tMin_90m = 0
const tMax_90m = 3
const tNum_90m = 100
const tVals_90m = LinRange(tMin_90m,tMax_90m,tNum_90m)
const tVals10_90m = LinRange(tMin_90m,tMax_90m,10)



#||||||||||||||||||||||/
#MARK: Travelling wave
#||||||||||||||||||||||\

k = 6.9*10^7
C = 10^7




#||--Demonstration of the match and shape--||#
graphName = "Demonstration_TravellingWave"

leftTZ = TransitionZone(0.5, k, k, C, C)
xtz_list = [leftTZ]

parameters = [EI, m, xp, xtz_list, P, v_short]

function N(x)
  return 256
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)

if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals_short, tVals10_short, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end


SteadyDeformation1 = SteadyStateTravellingSolution(xVals_short,tVals10_short,parameters, k, C)
SteadyStateDeformations = [SteadyDeformation1]

#BUG
gr()
Graph10Times(deformations,xVals_short,tVals10_short,parameters,SteadyStateDeformations, false, colour2, colour1, "Figures", graphName)



#||--Errors--||#
function CompareToTravellingWave()
  return NaN #TODO
end
#Some investigation of error using varying parameters
#Direct quadrature

#Weeks' method
