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
#Save calculations to file?
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



#--> High speed short track

v = 30

xLeft = 0
xRight = 4
xNum = 100
xVals = LinRange(xLeft,xRight,xNum)

tMin = 0
tMin = 0.049
tMax = 1/6
tMax = 0.052
tNum = 100
tVals = LinRange(tMin,tMax,tNum)
tVals10 = LinRange(tMin,tMax,10)



k = 6.9*10^7
C = 10^7

graphName = "Example-WideAndFast"

TZ1 = TransitionZone(1, k, 2*k, C, 2*C)
TZ2 = AddConsistentTZ(1.5, TZ1, k, C)
TZ3 = AddConsistentTZ(2.5, TZ1, 2*k, 2*C)
TZ4 = AddConsistentTZ(3, TZ1, k, C)
xtz_list = [TZ1, TZ2, TZ3, TZ4]

parameters = [EI, m, xp, xtz_list, P, v]
M = 104
setprecision(BigFloat,  6*M)

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> GWRImplementation_memo(xVals,tVals,LaplaceSpaceFunction, parameters, M)
#inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> GWRImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, M)

if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals, tVals10, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end


SteadyDeformation1 = SteadyStateTravellingSolution(xVals,tVals10,parameters, k, C)
SteadyStateDeformations = [SteadyDeformation1]

#BUG
gr()
Graph10Times(deformations,xVals,tVals10,parameters,SteadyStateDeformations, [-0.00005, 0.00005], colours, "", graphName)
