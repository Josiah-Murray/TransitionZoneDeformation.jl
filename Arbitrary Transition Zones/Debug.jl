using WAV #To play a sound when computation finished
using Plots
using JLD

figurePath = "C:/Users/joemu/Documents/PhD/Julia/Moving Point Force Storage/Figures"
DataPath = "C:/Users/joemu/Documents/PhD/Julia/Moving Point Force Storage/Data"

include("Functions - Inversion schemes.jl")
include("Functions - MPF_ArbitraryTZs.jl")
include("Functions BigFloat - MPF_ArbitraryTZs.jl")
include("Graphing - MPF_ArbitraryTZs.jl")



#||--More realistic case--||#

const EI = 2.3*10^3 #Beam stiffness
const m = 48 #Mass per unit length of beam
const P = -10^4 #Force conveyed by point load

const xp = 0 #Starting position of point force

const xLeft_Full = 0
const xRight_Full = 30
const xNum_Full = 100
const xVals_Full = LinRange(xLeft_Full,xRight_Full,xNum_Full)

const tMin_Full = 0
const tMax_Full = 1
const tNum_Full = 100
const tVals_Full = LinRange(tMin_Full,tMax_Full,tNum_Full)

tVals10_Full = LinRange(tMin_Full,tMax_Full,10)

colour1 = RGB(0.19215685, 0.27843137, 0.3333333)
colour2 = RGB(0.14901960784313725, 0.6274509803921569, 0.8549019607843137)

colour3 = RGB(0, 0.27450980392156865, 0.4980392156862745)
colour4 = RGB(0.6470588235294118, 0.8, 0.5098039215686274)

redLineColour = RGB(0.5,0.2,0.15)





##

k = 6.9*10^7
C = 10^7
v = 30



leftTZ = TransitionZone(5, k, 2*k, C, 2*C)
xtz_list = [leftTZ]


parameters= [EI, m, xp, xtz_list, P, v]

laplaceParameters = [50,100]
sVals = LinRange(-50,50,100)

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)

deformations = @time CalcDynamicDeformation(xVals_Full,tVals10_Full, parameters, inversionMethod)
graphName = "LongIntervalWeeks1024"
#save(DataPath*"/"*graphName*".jld", "deformations",deformations)

#deformations = load(DataPath*"/"*graphName*".jld")["deformations"]



SteadyDeformation1 = SteadyStateTravellingSolution(xVals_Full,tVals10_Full,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolutionBF(xVals_Full,tVals10_Full,parameters, 2*k, 2*C)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]


gr()
Graph10Times(deformations,xVals_Full,tVals10_Full,parameters,SteadyStateDeformations, [-0.00005,0.00005], colour2, colour1, figurePath, "temp")
