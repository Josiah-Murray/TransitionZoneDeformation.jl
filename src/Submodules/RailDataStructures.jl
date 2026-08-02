module RailDataStructures



#MARK: Parameters and TZs
struct RailParameters
  EI
  m
  transitionZones
end

"""
    RailParameters(EI = , m = , transitionZones = )

Create a `RailParameters` object which stores the bending stiffness (`EI`) and mass per unit length (`m`) of the beam, as well as a list of transition zones (see [`TransitionZone`](@ref) and [`AppendTransitionZone`](@ref))
"""
function RailParameters(
  ; EI
  , m
  , transitionZones)
  return RailParameters(EI, m, transitionZones)
end

struct TransitionZone
  k_left
  k_right
  C_left
  C_right
  position
end

"""
    TransitionZone(position = , k_left = , k_right = , C_left = 0, C_right = 0)

Create a `TransitionZone` object which stores the x coordinate (`position`) of a transition zone, as well as the track modulus and foundation damping on either side (`k_left`, `k_right`, `C_left`, and `C_right`, respectively).
The damping coefficients default to zero to streamline the use of the function in the steady-state case.
In general, the transition zone objects should be created with the use of [`AppendTransitionZone`](@ref).
"""
function TransitionZone(
  ; k_left
  , k_right
  , C_left = 0
  , C_right = 0
  , position
)
  return TransitionZone(k_left, k_right, C_left,  C_right, position)
end



end
