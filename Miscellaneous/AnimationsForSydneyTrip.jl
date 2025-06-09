using WAV #To play a sound when computation finished
using Plots
using JLD

println(pwd())

DataPath = "Paper/Data"
functionFolder = "../Arbitrary Transition Zones"

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


k = 6.9*10^7
C = 10^7




#MARK: Single TZ Weeks




graphName = "SingleTZ_t100_Weeks_4096"

leftTZ = TransitionZone(0.3, k, 2*k, C, 2*C)
xtz_list = [leftTZ]

parameters = [EI, m, xp, xtz_list, P, v_short]

function N(x)
  return 4096
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)

if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals_short, tVals_short, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end

SteadyDeformation1 = SteadyStateTravellingSolution(xVals_short,tVals_short,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals_short,tVals_short,parameters, 2*k, 2*C)
SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]



gr()
yLims = [-0.00055, 0.0001]
backgroundColour_min = background_light
backgroundColour_max = background_dark
@gif for ti in eachindex(tVals_short)
  p = plot()

    #Colour in background for transition zones
    colourGraphSegment(p,xtz_list, xLeft_short, xRight_short, yLims, backgroundColour_min, backgroundColour_max)

    #plot lines for transition zones
    for tz in xtz_list
      p = vline!([tz.location],color=colours[3], width=1, framestyle=:box)
    end

    #gradLineColour = CustomGradient(i/length(tIndexes), colour2,colour1)

    #Plot line for moving point force
    p = vline!([v_short*(tVals_short[ti]) + xp], width=3, color=colours[4])




    p = plot!(xVals_short, deformations[ti,:], lc = colours[1], lw=4, ylimits = yLims, xlimits = (xLeft_short,xRight_short))
    if(SteadyStateDeformations != false)
      for j in eachindex(SteadyStateDeformations)
        p = plot!(xVals_short, SteadyStateDeformations[j][ti,:], linestyle = :dash, color = colours[2], lw=2.5)
      end
    end
    p = plot!(framestyle=:box)
    title!(string(round(tVals_short[ti], digits=3))*" s" , titlefontsize = 12,line=-10)
    xlabel!("Beam coordinate (m)", xguidefontsize=10)
    ylabel!("Deformation (m)", yguidefontsize=10)
end



#MARK: Two steps TZ Weeks



graphName = "TwoStepsTZ_t100_Weeks_1024"

TZ1 = TransitionZone(0.1, k, 2k, C, 2C)
TZ2 = AddConsistentTZ(0.3, TZ1, k, C)
TZ3 = AddConsistentTZ(0.5, TZ2, 2k, 2C)
TZ4 = AddConsistentTZ(0.7, TZ3, k, C)


xtz_list = [TZ1,TZ2,TZ3,TZ4]

parameters = [EI, m, xp, xtz_list, P, v_short]

function N(x)
  return 1024
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)

if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals_short, tVals_short, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end


SteadyDeformation1 = SteadyStateTravellingSolution(xVals_short,tVals_short,parameters, k, C)
SteadyDeformation2 = SteadyStateTravellingSolution(xVals_short,tVals_short,parameters, 2k, 2C)
SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]


gr()
yLims = [-0.00055, 0.0001]
backgroundColour_min = background_light
backgroundColour_max = background_dark
@gif for ti in eachindex(tVals_short)
  p = plot()

    #Colour in background for transition zones
    colourGraphSegment(p,xtz_list, xLeft_short, xRight_short, yLims, backgroundColour_min, backgroundColour_max)

    #plot lines for transition zones
    for tz in xtz_list
      p = vline!([tz.location],color=colours[3], width=1, framestyle=:box)
    end

    #gradLineColour = CustomGradient(i/length(tIndexes), colour2,colour1)

    #Plot line for moving point force
    p = vline!([v_short*(tVals_short[ti]) + xp], width=3, color=colours[4])




    p = plot!(xVals_short, deformations[ti,:], lc = colours[1], lw=4, ylimits = yLims, xlimits = (xLeft_short,xRight_short))
    if(SteadyStateDeformations != false)
      for j in eachindex(SteadyStateDeformations)
        p = plot!(xVals_short, SteadyStateDeformations[j][ti,:], linestyle = :dash, color = colours[2], lw=2.5)
      end
    end
    p = plot!(framestyle=:box)
    title!(string(round(tVals_short[ti], digits=3))*" s" , titlefontsize = 12,line=-10)
    xlabel!("Beam coordinate (m)", xguidefontsize=10)
    ylabel!("Deformation (m)", yguidefontsize=10)
end






#||--#MARK: fast Demonstration of the match and shape

graphName = "Demonstration_fast_t100_TravellingWave"

leftTZ = TransitionZone(0.5, k, k, C, C)
xtz_list = [leftTZ]

parameters = [EI, m, xp, xtz_list, P, v_fast]

function N(x)
  return 4096
end
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)

if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals_fast, tVals_fast, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end


SteadyDeformation1 = SteadyStateTravellingSolution(xVals_fast,tVals_fast,parameters, k, C)
SteadyStateDeformations = [SteadyDeformation1]

save(DataPath*"/"*graphName*".jld", "deformations", deformations)


gr()
yLims = [-4*10^-5, 1*10^-5]
backgroundColour_min = background_light
backgroundColour_max = background_dark
@gif for ti in eachindex(tVals_fast)
  p = plot()

    #Colour in background for transition zones
    colourGraphSegment(p,xtz_list, xLeft_fast, xRight_fast, yLims, backgroundColour_min, backgroundColour_max)

    #plot lines for transition zones
    for tz in xtz_list
      p = vline!([tz.location],color=colours[3], width=1, framestyle=:box)
    end

    #gradLineColour = CustomGradient(i/length(tIndexes), colour2,colour1)

    #Plot line for moving point force
    p = vline!([v_fast*(tVals_fast[ti]) + xp], width=3, color=colours[4])




    p = plot!(xVals_fast, deformations[ti,:], lc = colours[1], lw=4, ylimits = yLims, xlimits = (xLeft_fast,xRight_fast))
    if(SteadyStateDeformations != false)
      for j in eachindex(SteadyStateDeformations)
        p = plot!(xVals_fast, SteadyStateDeformations[j][ti,:], linestyle = :dash, color = colours[2], lw=2.5)
      end
    end
    p = plot!(framestyle=:box)
    title!(string(round(tVals_fast[ti], digits=3))*" s" , titlefontsize = 12,line=-10)
    xlabel!("Beam coordinate (m)", xguidefontsize=10)
    ylabel!("Deformation (m)", yguidefontsize=10)
end
