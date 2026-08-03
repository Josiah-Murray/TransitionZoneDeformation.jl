using Plots
import TransitionZoneDeformation as TZD

#Define the list of transition zones, each of which is a custom data type.
transitionZones = [TZD.RailDataStructures.TransitionZone( k_left = 4, k_right = 1, C_left = 0.5, C_right = 1, position = 0 )]

#Add another transition zone at position 5 and with k_right  = 0.5 and C_right = 2. The left properties are added automatically based on the right properties of the previous transition zone.
transitionZones = TZD.AppendTransitionZone(transitionZones, 5, 0.5, 2)

#Create a RailParameters object that contains the properties of the beam and the list of transition zones.
beamParameters = TZD.RailDataStructures.RailParameters( EI = 1.0, m = 1.0, transitionZones = transitionZones )

#Define the x values for which the solution should be found.
xVals = LinRange(-10, 10, 100)

#Location of point force.
xp = 3

#Calculate the steady-state deformation for the system.
deformation = TZD.SteadyStateWithTransitionZones.CalculateDeformation(xVals, xp, beamParameters)




p = plot(xVals
  , -deformation
  , xlabel = "x"
  , xlims = (-10, 10)
  , ylims = (-0.4, 0.1)
  , ylabel = "Deformation"
  , background_color = RGB(0.9, 0.9, 0.95)
  , linestyle = :solid
  , linewidth = 4
  , color = RGB(0.2, 0.2, 0.1)
  , frame_style = :box
  , legend = false
    )
for tz in transitionZones
  p = vline!([tz.position], linestyle = :dash, color = :black, linewidth = 1)
end
display(p)
