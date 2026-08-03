using Plots
import TransitionZoneDeformation as TZD

#Define the list of transition zones, each of which is a custom data type.
transitionZones = [TZD.RailDataStructures.TransitionZone( k_left = 1, k_right = 1.5, C_left = 1, C_right = 1.5, position = -1 )]

#Add another transition zone at position 5 and with k_right  = 0.5 and C_right = 2. The left properties are added automatically based on the right properties of the previous transition zone.
transitionZones = TZD.AppendTransitionZone(transitionZones, 1, 1, 1)

#Create a RailParameters object that contains the properties of the beam and the list of transition zones.
beamParameters = TZD.RailDataStructures.RailParameters( EI = 1.0, m = 1.0, transitionZones = transitionZones )

#Speed and location of point force.
v = 1
xp = -3

#Time and space values for which to find the solution. (Note that the GWR implementation applied later cannot be evaluated for time zero).
tVals = LinRange{BigFloat}(0.01, 6, 20) #BigFloat data type is recommended for the GWR algorithm.
xVals = LinRange(-10, 10, 100)

#Define the Laplace-domain solution for the deformation of the beam. This will be an input to an inversion algorithm which will then return the time-domain solution.
LaplaceDomainFunction = s -> TZD.DynamicWithTransitionZones.LaplaceDomainFunction(xVals, s, xp, v, beamParameters)

#Apply a numerical inversion algorithm to find the time-domain deformation of the system.
M = 40
setprecision(3*M) #Recommended precision for the GWR algorithm.
deformation = TZD.LaplaceInversionImplementations.GWRImplementation(xVals, tVals, beamParameters, LaplaceDomainFunction; M = 20)

#Produce an animation of the solution.
@gif for (ti,t) in enumerate(tVals)
  p = plot(xVals
    , -deformation[ :, ti]
    , xlabel = "x"
    , xlims = (-10, 10)
    , ylims = (-0.4, 0.1)
    , ylabel = "Deformation"
    , background_color = RGB(0.9, 0.9, 0.95)
    , linestyle = :solid
    , linewidth = 4
    , color = RGB(0.2, 0.2, 0.1)
    , frame_style = :box
    , title = "t = $(round(Float64(t), digits = 2))"
    , legend = false
      )
  for tz in transitionZones
    p = vline!([tz.position], linestyle = :dash, color = :black, linewidth = 1)
  end
end
