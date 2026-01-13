using WAV #To play a sound when computation finished
using Plots
using JLD

println(pwd())

DataPath = "Paper--Moving point force/Data"
functionFolder = "../Constant Point Force/"
figureFolder = "Paper--Moving point force/Figures"

include(joinpath(functionFolder, "Functions - Inversion schemes.jl"))
include(joinpath(functionFolder, "Functions - MPF_ArbitraryTZs.jl"))
include(joinpath(functionFolder, "Graphing - MPF_ArbitraryTZs.jl"))

#||||||||||||||||||||||/
#MARK: Colours
#||||||||||||||||||||||\

solidLine = RGB(0.9,0.1,0.4)
solidLine = RGB(0.8,0.2,0.4)
solidLine = RGB(0,0,0)
#solidLine = RGB(1,0,0.4)
#solidLine = RGB(0.8,0.4,1)
#solidLine = RGB(0.9,0.2,0.2)

brokenLine = RGB(0.84901960784313725, 0.8274509803921569, 0.2549019607843137)
brokenLine = RGB(0,0,0)
brokenLine = RGB(0.8,0.2,0.4)
brokenLine = RGB(0.9,0.2,0.4)


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
LoadFlag = true
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



#||||| A version of the CalcDynamicDeformation function which ensures that the
#Laplace space function uses the Float64 data type for its calculations.
#This is for evaluation of the GWR algorithm.

function CalcDynamicDeformationF64(xVals,tVals, parameters, inversionMethod)
  if parameters[end] == 0
    #TODO Implement v=0?
    @error "CalcDynamicDeformation unimplemented for v=0"
    return NaN
  end
  deformations = inversionMethod(xVals, tVals, LaplaceSpaceFunctionMovingPointForce1TZF64, parameters)
  return deformations
end
#||||||||||||||||||||||||||||||||||||||||||||

#See Above.
function LaplaceSpaceFunctionMovingPointForce1TZF64(x, s, parameters; Coeff_solver = CoefficientSolverMovingPointForce1TZ)

  Dt = Float64
  x  = convert(Dt, x)

  #We convert here so that we don't have to know ahead of time what type we'll need and can let the inversion implementation decide.
  EI, m, xp, xtz_list, P, v = parameters
  EI = convert(Dt, EI)
  m  = convert(Dt,  m)
  xp = convert(Dt, xp)
  P  = convert(Dt,  P)
  v  = convert(Dt,  v)
  parameters = [EI, m, xp, xtz_list, P, v]

  #Set the correct parameters.
  #Default to the furthest right parameters, override if we are to the left.
  k = xtz_list[end].k_right
  C = xtz_list[end].C_right
  for tz in xtz_list
    if real(x) < tz.location
      k = tz.k_left
      C = tz.C_left
      break
    end
  end


  #||--Find r values--||#
  #Only used for graphing appropriately
  rVals = RValues(s, EI, m, C, k)
  r1, r2, r3, r4 = rVals


  #||--Find b values--||#
  #Solve for the undetermined coefficients by solving a linear system
  bVals = Coeff_solver(s, parameters)

  #Find correct segment to isolate the coefficients we need for this x value
  segment = 1
  pointForceAdded = false
  for tz in xtz_list
    if ~pointForceAdded && xp < tz.location
      pointForceAdded = true
      if real(x) < xp
        break
      else
        segment += 1
      end
    end
    if real(x) < tz.location
      break
    else
      segment += 1
    end
  end
  if ~pointForceAdded
    segment += 1
  end

  #Pull out correct b values
  b1,b2,b3,b4 =  bVals[4*(segment-1) + 1], bVals[4*(segment-1) + 2], bVals[4*(segment-1) + 3], bVals[4*(segment-1) + 4]


  #||--Construct Laplace-space function--||#

  ŷ = b1*exp(r1*x) + b2*exp(r2*x) + b3*exp(r3*x) + b4*exp(r4*x) + ICResponse(x,s, parameters, C, k, 0) +   (P/abs(v))*(exp(-s*(x-xp)/v)  /  (  ( (EI*s^4)/v^4 ) +m*s^2 + C*s + k )  )*Heaviside((x-xp)/v)

  return ŷ

end





#||||||||||||||||||||||/
#||||--Travelling wave
#||||||||||||||||||||||\

k = 6.9*10^7
C = 10^7




#MARK: Multiple Transition Zones



numberOfTZs = 100

graphName = "GradiatedTransitionZones_"*string(numberOfTZs)

function GradiatedTransitionZones(numberOfTZs, k_start, k_end, C_start, C_end)
  xtz_list = TransitionZone[]

  k_left = k_start
  C_left = C_start

  for i in 1:numberOfTZs
    location = i/(numberOfTZs+1)
    k_right = k_start + (k_end - k_start)*((i)/(numberOfTZs))
    C_right = C_start + (C_end - C_start)*((i)/(numberOfTZs))
    push!(xtz_list, TransitionZone(location, k_left, k_right, C_left, C_right))
    k_left = k_right
    C_left = C_right
  end
  return xtz_list
end



xtz_list = GradiatedTransitionZones(numberOfTZs, k, 2k, C, 2C)

parameters = [EI, m, xp, xtz_list, P, v_short]

function N(x)
  return 1024
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
SteadyDeformation2 = SteadyStateTravellingSolution(xVals_short,tVals10_short,parameters, 2k, 2C)
SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]


gr()
Graph10Times(deformations,xVals_short,tVals10_short,parameters,SteadyStateDeformations, false, colours, figureFolder*"/", graphName, blank = true)







numberOfTZs = 100

graphName = "GWR_GradiatedTransitionZones_"*string(numberOfTZs)

function GradiatedTransitionZones(numberOfTZs, k_start, k_end, C_start, C_end)
  xtz_list = TransitionZone[]

  k_left = k_start
  C_left = C_start

  for i in 1:numberOfTZs
    location = i/(numberOfTZs+1)
    k_right = k_start + (k_end - k_start)*((i)/(numberOfTZs))
    C_right = C_start + (C_end - C_start)*((i)/(numberOfTZs))
    push!(xtz_list, TransitionZone(location, k_left, k_right, C_left, C_right))
    k_left = k_right
    C_left = C_right
  end
  return xtz_list
end



xtz_list = GradiatedTransitionZones(numberOfTZs, k, 2k, C, 2C)

parameters = [EI, m, xp, xtz_list, P, v_short]


M = 40
setprecision(BigFloat, M*3)
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> GWRImplementation_memo(xVals,tVals,LaplaceSpaceFunction, parameters, M)


if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals_short, tVals10_short, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end


SteadyDeformation1 = SteadyStateTravellingSolution(xVals_short,tVals10_short,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals_short,tVals10_short,parameters, 2k, 2C)
SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]


gr()
Graph10Times(deformations,xVals_short,tVals10_short,parameters,SteadyStateDeformations, false, colours, figureFolder*"/", graphName, blank = true)






y, fs = wavread(raw"C:/Windows/Media/Ring01.wav")
#wavplay(y, fs)

println("ENDED")
