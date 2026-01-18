module RailwayDeformation

#External dependencies
using PolynomialRoots



#Internal dependencies
include(joinpath(@__DIR__, "Submodules/RailDataStructures.jl"))
using .RailDataStructures

include(joinpath(@__DIR__, "Submodules/SteadyStateWithTransitionZones.jl"))
using .SteadyStateWithTransitionZones

include(joinpath(@__DIR__, "Submodules/Utilities.jl"))
using .Utilities

"""
    AppendTransitionZone(transitionZoneArray, position, k_right, C_right)

Append a transition zone to the end of an existing list such that the 'left' foundation properties are consistent with the 'right' foundation properties of the previous transition zone.
Note that the transition zones must be added in order from left to right.

# Examples
```julia-repl
julia> tzList = [RailwayDeformation.TransitionZone(position = 1, k_left = 1, k_right = 2)]
1-element Vector{Main.RailwayDeformation.TransitionZone}:
 Main.RailwayDeformation.TransitionZone(1, 2, 0, 0, 1)

 julia> tzList = RailwayDeformation.AppendTransitionZone(tzList, 1.5, 1, 2)
 Main.RailwayDeformation.TransitionZone(1, 2, 0, 0, 1)
 Main.RailwayDeformation.TransitionZone(2, 1, 0, 2, 1.5)
```
"""
function AppendTransitionZone(transitionZoneArray, position, k_right, C_right)
  if position <= transitionZoneArray[end].position
    @error "Transition zones should be added left to right. Currently `position`= $position is less than or equal to `transitionZoneArray[end].position`= $(transitionZoneArray[end].position)"
    return nothing
  end
  return [transitionZoneArray; RailDataStructures.TransitionZone(transitionZoneArray[end].k_right, k_right, transitionZoneArray[end].C_right, C_right, position)]
end




end
