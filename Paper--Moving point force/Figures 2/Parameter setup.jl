using WAV #To play a sound when computation finished
using Plots
using JLD

println(pwd())

DataPath = "Data"
functionFolder = "../../Arbitrary Transition Zones"

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



solidLine = RGB(1,1,1)

brokenLine = RGB(0.84901960784313725, 0.8274509803921569, 0.2549019607843137)
brokenLine = RGB(0.5,0.2,0.4)


transitionZoneLine = RGB(0,0,0)
pointForceLine = RGB(0,0,0.25)

background_light = RGB(0.2,0.2,0.3)
background_dark = RGB(0.05,0.05,0.2)

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


#||--Model constants--||#

const EI = 2.3*10^3 #Beam stiffness
const m = 48 #Mass per unit length of beam
const P = -10^4 #Force conveyed by point load

const xp = 0 #Starting position of point force

const k_default = 6.9*10^7
const C_default = 10^7

#||--Space time setups--||#

#--> Slow track

const v_slow = 1

const xLeft_slow = -1
const xRight_slow = 1
const xNum_slow = 100
const xVals_slow = LinRange(xLeft_slow,xRight_slow,xNum_slow)

const tMin_slow = 0
const tMax_slow = 1
const tNum_slow = 100
const tVals_slow = LinRange(tMin_slow,tMax_slow,tNum_slow)
const tVals10_slow = LinRange(tMin_slow,tMax_slow,10)


#--> Fast track

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
