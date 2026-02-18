using Revise
using Plots
include(joinpath(@__DIR__, "../RailwayDeformation.jl"))
import .RailwayDeformation as RD


transitionZones = [RD.RailDataStructures.TransitionZone( k_left = 1, k_right = 1.5, C_left = 1, C_right = 1.5, position = -1 )]
transitionZones = RD.AppendTransitionZone(transitionZones, 1, 1, 1)
#transitionZones = RD.AppendTransitionZone(transitionZones, 3, 2, 2)

beamParameters = RD.RailDataStructures.RailParameters( EI = 1.0, m = 1.0, transitionZones = transitionZones )
v = -1
xp = 0
tVal = 1.0
xVals = LinRange(-10, 10, 100)

LaplaceDomainFunction = s -> RD.DynamicWithTransitionZones.LaplaceDomainFunction(xVals, s, xp, v, beamParameters)

deformation = RD.LaplaceInversionImplementations.GWRImplementation(xVals, [tVal], beamParameters, LaplaceDomainFunction; M = 20, shift_parameter = 0.01)[ :, 1]
#deformation  = abs.(LaplaceDomainFunction(0.5+1im))
fillColours = [RGB(1,0,0) for j in eachindex(xVals)]



kMax =  0
C_Max = 0
kmin = Inf
C_min = Inf
for tz in transitionZones
  if tz.k_left > kMax
    global kMax = tz.k_left
  end
  if tz.k_right > kMax
    global kMax = tz.k_right
  end
  if tz.C_left > C_Max
    global C_Max = tz.C_left
  end
  if tz.C_right > C_Max
    global C_Max = tz.C_right
  end
  if tz.k_left < kmin
    global kmin = tz.k_left
  end
  if tz.k_right < kmin
    global kmin = tz.k_right
  end
  if tz.C_left < C_min
    global C_min = tz.C_left
  end
  if tz.C_right < C_min
    global C_min = tz.C_right
  end
end


colourDark = [0.2, 0.2, 0.18]
colourLight = [0.8, 0.8, 0.85]

gradientFunc_k = k -> (k - kmin)/(kMax - kmin)
gradientFunc_C = C -> (C - C_min)/(C_Max - C_min)

#=
for (i, x) in enumerate(xVals)
  for tz in transitionZones
    if x < tz.position
      fillColours[i] = RGB(
        colourLight[1] - gradientFunc_k(tz.k_left)*(colourLight[1]-colourDark[1]) - gradientFunc_C(tz.C_left)*(colourLight[1]-colourDark[1])
        , colourLight[2] - gradientFunc_k(tz.k_left)*(colourLight[2]-colourDark[2]) - gradientFunc_C(tz.C_left)*(colourLight[2]-colourDark[2])
        , colourLight[3] - gradientFunc_k(tz.k_left)*(colourLight[3]-colourDark[3]) - gradientFunc_C(tz.C_left)*(colourLight[3]-colourDark[3])
      )
      break
    elseif x >= transitionZones[end].position
      fillColours[i] = RGB(
        colourLight[1] - gradientFunc_k(tz.k_right)*(colourLight[1]-colourDark[1]) - gradientFunc_C(tz.C_right)*(colourLight[1]-colourDark[1])
        , colourLight[2] - gradientFunc_k(tz.k_right)*(colourLight[2]-colourDark[2]) - gradientFunc_C(tz.C_right)*(colourLight[2]-colourDark[2])
        , colourLight[3] - gradientFunc_k(tz.k_right)*(colourLight[3]-colourDark[3]) - gradientFunc_C(tz.C_right)*(colourLight[3]-colourDark[3])
      )
    end
  end
end
=#




p = plot(xVals
  , -deformation
  , xlabel = "x"
  , xlims = (-10, 10)
  , ylims = (-0.4, 0.1)
  , ylabel = "Deformation"
  , background_color = RGB(0.9, 0.9, 0.95)
  #, fill = (-0.5, fillColours)
  , linestyle = :solid
  , linewidth = 4
  , color = RGB(0.2, 0.2, 0.1)
  , frame_style = :box
    )
for tz in transitionZones
  p = vline!([tz.position], linestyle = :dash, color = :white, linewidth = 1)
end
display(p)
