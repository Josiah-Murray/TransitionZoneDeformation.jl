using WAV #To play a sound when computation finished
using Plots
using JLD

println(pwd())

DataPath = ""
functionFolder = "Arbitrary Transition Zones"

include(joinpath(functionFolder, "Functions - Inversion schemes.jl"))
include(joinpath(functionFolder, "Functions - MPF_ArbitraryTZs.jl"))
include(joinpath(functionFolder, "Graphing - MPF_ArbitraryTZs.jl"))

#||||||||||||||||||||||/
#MARK: Colours
#||||||||||||||||||||||\

solidLine = RGB(0.5,0.2,0.4)
#solidLine = RGB(1,0,0.4)
#solidLine = RGB(0.8,0.4,1)
#solidLine = RGB(0.9,0.2,0.2)

brokenLine = RGB(0.84901960784313725, 0.8274509803921569, 0.2549019607843137)
brokenLine = RGB(0,0,0)


transitionZoneLine = RGB(0,0,0)
pointForceLine = RGB(0,0,0.25)

background_light = RGB(0.82,0.9,0.91)
background_dark = RGB(0.65,0.75,0.75)

colours = [solidLine,
brokenLine,
transitionZoneLine,
pointForceLine,
background_light,
background_dark]




#||||||||||||||||||||||/
#MARK: Flags
#||||||||||||||||||||||\

#Perform fresh calculations or load from old?
LoadFlag = false
#Save calculations from file?
SaveFlag = false


#||||||||||||||||||||||/
#MARK: Constants
#||||||||||||||||||||||\



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

#--> High speed short track

const v_fast = 30

const xLeft_fast = -1
const xRight_fast = 1
const xNum_fast = 100
const xVals_fast = LinRange(xLeft_fast,xRight_fast,xNum_fast)

const tMin_fast = 0
const tMax_fast = 1/30
const tNum_fast = 100
const tVals_fast = LinRange(tMin_fast,tMax_fast,tNum_fast)
const tVals10_fast = LinRange(tMin_fast,tMax_fast,10)


#MARK: Functions

function CalculateAccelerations(deformations, tVals)
  Δt = tVals[2]-tVals[1]
  central = (y1,y2,y3) -> (y1 -2*y2  + y3)/Δt^2
  forward = (y1,y2,y3,y4) -> (2*y1 - 5*y2 +4*y3 - y4)/Δt^2
  backward = (y1,y2,y3,y4) -> (-y1 + 4*y2 - 5*y3 +2*y4)/Δt^2
  accelerations = zeros(size(deformations))
  for xi in eachindex(deformations[1,:])
    accelerations[1,xi] = forward(deformations[1,xi], deformations[2,xi], deformations[3,xi], deformations[4,xi])
    for ti in 2:length(tVals)-1
      accelerations[ti,xi] = central(deformations[ti-1,xi], deformations[ti,xi], deformations[ti+1,xi])
    end
    accelerations[end,xi] = backward(deformations[end-3,xi], deformations[end-2,xi],deformations[end-1,xi],deformations[end,xi])
  end
  return accelerations
end





#||||||||||||||||||||||/
#||||--Travelling wave
#||||||||||||||||||||||\

k = 6.9*10^7
C = 10^7






#MARK:||COMPARISON||







##
#||--#MARK: Demonstration of the match and shape
graphName = "Test"

leftTZ = TransitionZone(0.5, k, 2k, C, 2C)
xtz_list = [leftTZ]

parameters = [EI, m, xp, xtz_list, P, v_fast]

function N(x)
  return 256
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)

if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals_fast, tVals10_fast, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end


SteadyDeformation1 = SteadyStateTravellingSolution(xVals_fast,tVals10_fast,parameters, k, C)
SteadyStateDeformations = [SteadyDeformation1]

#BUG
gr()
Graph10Times(deformations,xVals_fast,tVals10_fast,parameters,SteadyStateDeformations, false, colours, "", graphName)
