#|||||||||||||||#
#||--READ ME--||#
#|||||||||||||||#

#=
This file contains functions used for solving for the dynamic deformation of an
infinite Euler-Bernoulli beam on a piecewise constant viscoelastic foundation
subject to a moving mass-spring system.

The methods are based on the paper:
L
  Fărăgău, A. B., Mazilu, T., Metrikine, A. v, Lu, T., & van Dalen, K. N. (2021).
  Transition radiation in an infinite one-dimensional structure interacting
  with a moving oscillator—the Green’s function method.
  Journal of Sound and Vibration, 492, 115804.
 https://doi.org/10.1016/j.jsv.2020.115804
Γ

The function  'CalcDynamicDeformationMM' is the one that is intended to be called by a user.
=#

#MARK: Structs

#Struct which holds information about the transition zones.
#All fields should contain be numerical types (i.e. Int, Float64, BigFloat, ArbFloat, etc.).
struct SmoothTransitionZone
  L<:Real #Length of the transition zone. (Begins at 0, ends at L)
  k_interpolant<:Function
  C_interpolant<:Function
  k_left<:Real
  k_right<:Real
  C_left<:Real
  C_right<:Real
end

#Fills in the left and right foundation parameters based on the interpolant functions.
function constructTransitionZone(L, k_interpolant, C_interpolant)
  k_left = k_interpolant(0)
  k_right = k_interpolant(L)
  C_left = C_interpolant(0)
  C_right = C_interpolant(L)
  return SmoothTransitionZone(L, k_interpolant, C_interpolant, k_left, k_right, C_left, C_right)
end

struct ParameterStruct
  EI<:Real #Flexural rigidity of the beam.
  ρ<:Real #Mass per unit length of the beam, scaled by EI.
  v<:Real #velocity of the moving mass.
  Q0<:Real #Dead weight of the moving mass.
  M<:Real #Mass of the moving mass.
  TZ<:SmoothTransitionZone #The transition zone object.
  u_0<:Real #Initial displacement of the moving mass.
  uDot_0<:Real #Initial velocity of the moving mass.
end




#----------------------------------------------------------------#
#MARK: Main Calculation
#The main function that should be called to perform the relevant calculations.
#Returns a tuple containing:
# 1. The dynamic deformation of the beam in the computational domain w_c(x,t) as a 2D array.
#    The first dimension corresponds to space, the second to time.  i.e. w_c[i,j] = w_c(x_i, t_j)
# 2. The displacement of the moving mass u(t) as a 1D array.
#Inputs:
# xVals: A 1D array of the spatial points at which to calculate the deformation.
# tVals: A 1D array of the time points at which to calculate the deformation.
# parameters: A ParameterStruct object containing the relevant model parameters.
function CalcDynamicDeformationMM(xVals, tVals, parameters)

  #Extract parameters.
  EI = parameters.EI
  ρ = parameters.ρ
  v = parameters.v
  Q0 = parameters.Q0
  M = parameters.M
  TZ = parameters.TZ
  L = TZ.L


  #Initialise the beam deformation array.
  w_c = zeros(eltype(xVals), length(xVals), length(tVals))
  #Initialise the mass displacement array.
  u = zeros(eltype(xVals), length(tVals))

  #Loop over each time step and space step to calculate the deformation.
  for (ti, t) in enumerate(tVals)

    #Setup time discretisation for use with Q_n calculation.
    n = 100 #TODO: Figure out appropriate discretisation of interval [0,t].
    subTVals = LinRange(0, t, n+1)
    Δτ = subTVals[2]-subTVals[1]

    #Calculate space indendent Green’s function values. Used in calculation of {Q_i}.
    g_vt_vτ_tmτ = [GreensFunction(v*t, v*τ, t-τ, parameters) for τ in subTVals]

    #Calculate {Q_i} i=0,...,n
    #Calculate reaction forces at each sub-time step.
    Q = CalculateReactionForces(g_vt_vτ_tmτ, subTVals, parameters)


    for (xi, x) in enumerate(xVals)

      #Calculate space dependent Green’s function values. Used in calculation of w_c(x,t).
      g_x_vτ_tmτ = [GreensFunction(x, v*τ, t-τ, parameters) for τ in subTVals]


      #TODO: Calculate w_c^L
      w_c_L = 0 #TODO: Replace with actual calculation.

      #TODO: Calculate w_c^IC
      w_c_IC = 0 #TODO: Replace with actual calculation.

      #Calculate w_c(x,t)
      w_c[xi, ti] = sum([ g_x_vτ_tmτ[i]*Q[end]*Heaviside(L-vτ)*(Δτ) for (τi, τ) in enumerate(subTVals)]) + w_c_L + w_c_IC



    end
    #Calculate u
    u[ti] = (-Q0/(2*M))*t^2 + u_0 + uDot_0*t + (1/M)*sum([ Q_n[i]*(t-τ)*(Δτ) for (i, τ) in enumerate(subTVals)])

  end

  return w_c, u

end
#----------------------------------------------------------------#




#TODO: Implement.
#MARK: Green's Function
#Calculates the Green’s function for the beam at position x and time t due to a unit impulse at time τ.
#Inputs:
# x: The spatial position at which to calculate the Green’s function.
# t: The time at which to calculate the Green’s function.
# τ: The time at which the impulse occurs.
# parameters: A ParameterStruct object containing the relevant model parameters.
function GreensFunction(x, t, τ, parameters)
  #Extract parameters.
  EI = parameters.EI
  ρ = parameters.ρ
  v = parameters.v
  Q0 = parameters.Q0
  M = parameters.M
  TZ = parameters.TZ
  L = TZ.L

  return 0 #TODO: Implement
end

#MARK: Reaction Forces
#Reaction force calculation.
#Calculates the reaction forces {Q_i} at each sub-time step i=0,...,n.
#Inputs:
# g_vt_vτ_tmτ: A 1D array of the Green's function values g(v*t, v*τ_i, t-τ_i) for i=0,...,n.
# subTVals: A 1D array of the sub-time values τ_i for i=0,...,n.
# parameters: A ParameterStruct object containing the relevant model parameters.
function CalculateReactionForces(g, subTVals, parameters)

  #TODO: Fix up notation with n etc.

  #Extract parameters.
  EI = parameters.EI
  ρ = parameters.ρ
  v = parameters.v
  Q0 = parameters.Q0
  M = parameters.M
  TZ = parameters.TZ
  L = TZ.L
  u_0 = parameters.u_0
  uDot_0 = parameters.uDot_0

  Δt = subTVals[2]-subTVals[1]

  n = length(subTVals)-1

  Q = zeros(eltype(subTVals), n+1)
  Q[1] = Q0

  c_H = 1.1864*10^11 #Hertz constant as per Fărăgău et al. (2021).

  for i in 1:n
    t = subTVals[i+1]
    #Calculate R (a function of Q)
    R = Q_n -> (2*g[0+1] + g[1+1] - (Δt/M)*(3n-1) + t^2)*Q_n + sum( [(g[j-1+1] + 4*g[j+1] + g[j+1+1] - (Δt/M)*(6*n-6*j))*Q[j+1] for j=1:n-1 ] )
     + w_c_IC_vt_t -u_0 + uDot_0*t #TODO: Pass in w_c_IC(v*t, t) properly.

    #Calculate Q_i
    Q[i+1] = NonLinearSolve(Q_n ->  (Q_n/c_H)^(2/3) - R(Q_n)*Heaviside(R(Q_n)), Q[i]) #Initial guess is previous value. #BUG: 2/3 or 3/2?
  end





  return Q
end



#MARK: Utilities

#Heaviside step function.
function Heaviside(x)
  if x < 0
    return 0
  elseif x == 0
    return 0.5
  else
    return 1
  end
end
