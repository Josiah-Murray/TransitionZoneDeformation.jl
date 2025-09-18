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

#Import my Laplace inversion module.
include("LaplaceInversionAdress/ImportNILaplaceModule.jl")



#MARK: Structs

#Struct which holds information about the transition zones.
#All fields should contain be numerical types (i.e. Int, Float64, BigFloat, ArbFloat, etc.).
struct SmoothTransitionZone
  L #<:Real #Length of the transition zone. (Begins at 0, ends at L)
  k_interpolant#<:Function
  C_interpolant#<:Function
  k_left#<:Real
  k_right#<:Real
  C_left#<:Real
  C_right#<:Real
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
  EI#<:Real #Flexural rigidity of the beam.
  ρ#<:Real #Mass per unit length of the beam, scaled by EI.
  v#<:Real #velocity of the moving mass.
  Q0#<:Real #Dead weight of the moving mass.
  M#<:Real #Mass of the moving mass.
  TZ#<:SmoothTransitionZone #The transition zone object.
  u_0#<:Real #Initial displacement of the moving mass.
  uDot_0#<:Real #Initial velocity of the moving mass.
end




#----------------------------------------------------------------#
#MARK: Main Calculation
#The main function that should be called to perform the relevant calculations.
#Returns a tuple containing:
# 1. The dynamic deformation of the beam in the computational domain w_c(x,t) as a 2D array.
#    The first dimension corresponds to space, the second to time.  i.e. w_c[i,j] = w_c(x_i, t_j)
# 2. The displacement of the moving mass u(t) as a 1D array.
#Inputs:
# numXVals: The number of x-points in the computational domain. The computational domain is [0, L] where L is the length of the transition zone.
# parameters: A ParameterStruct object containing the relevant model parameters.
function CalcDynamicDeformationMM(numXVals, parameters)

  #Extract parameters.
  EI = parameters.EI
  ρ = parameters.ρ
  v = parameters.v
  Q0 = parameters.Q0
  M = parameters.M
  TZ = parameters.TZ
  L = TZ.L

  #Define discretisation.
  N = numXVals - 1 #Number of intervals.
  Δx = L/N #Spatial step size.
  timeStepOffset = 3 #Only calculate for every timeStepOffset time steps to reserve x-values for the necessary integrations in terms of τ.
  Δt = timeStepOffset*Δx/v #Time step size.
  xVals = LinRange(0, L, numXVals) #Spatial discretisation.
  tVals = LinRange(0, L/v, numXVals) #Time discretisation.


  #Initialise the beam deformation array.
  w_c = zeros(eltype(xVals), length(xVals), length(tVals))
  #Initialise the mass displacement array.
  u = zeros(eltype(xVals), length(tVals))

  #Loop over each time step and space step to calculate the deformation.
  for (ti, t) in enumerate(tVals)

    #Setup time discretisation for use with Q_n calculation.
    n = timeStepOffset*ti+1
    τVals = LinRange(0, t, n)
    Δτ = τVals[2]-τVals[1]

    #Calculate the Green's function for each x value, each ξ (which are themselves x values), and each t-τ value.
    #TODO: Memoise this?
    GreensFunctions = zeros(Complex{eltype(xVals)}, length(xVals), length(xVals), length(τVals)) #4D array to hold the Green's function values.
    for (τi, τ) in enumerate(τVals)
      GreensFunctions[:,:,τi] =  InvertedGreensFunction(xVals, xVals, t- τ, parameters)
    end
    #GreensFunctions = InvertedGreensFunction(xVals, xVals, t*ones(length(τVals))- τVals, parameters)
    #Calculate space indendent Green’s function values. Used in calculation of {Q_i}.
    #g_vt_vτ_tmτ = [InvertedGreensFunction(v*t, v*τ, t-τ, parameters) for τ in subTVals]
    g_vt_vτ_tmτ = [GreensFunctions[ti*timeStepOffset, τi, τi] for τi in eachindex(τVals)] #BUG: Is this right?
    #Calculate {Q_i} i=0,...,n
    #Calculate reaction forces at each sub-time step.
    Q = CalculateReactionForces(g_vt_vτ_tmτ, τVals, parameters)


    for (xi, x) in enumerate(xVals)

      #Calculate space dependent Green’s function values. Used in calculation of w_c(x,t).
      g_x_vτ_tmτ = [InvertedGreensFunction(x, v*τ, t-τ, parameters) for τ in τVals] #TODO: Change to use `GreensFunctions` array.


      #TODO: Calculate w_c^L
      w_c_L = 0

      #TODO: Calculate w_c^IC
      w_c_IC = 0
      #Calculate w_c(x,t)
      w_c[xi, ti] = sum([ g_x_vτ_tmτ[i]*Q[end]*Heaviside(L-vτ)*(Δτ) for (τi, τ) in enumerate(τVals)]) + w_c_L + w_c_IC
    end
    #Calculate u
    u[ti] = (-Q0/(2*M))*t^2 + u_0 + uDot_0*t + (1/M)*sum([ Q_n[i]*(t-τ)*(Δτ) for (i, τ) in enumerate(τVals)])

  end

  return w_c, u

end
#----------------------------------------------------------------#




#TODO: Implement non-reflective boundary conditions.
#MARK: Green's Function


#Laplace domain Green's function calculation.
#Calculates the Laplace domain Green's function G(x,ξ,s) for the beam at all xVals due to a unit impulse at position ξ.
#Inputs:
# xVals: The positions at which to calculate the Green's function.
# ξ_index: The index of the impulse position in the spatial discretisation. (Note that due to the discretisation, ξ must be one of the x values).
# s: The Laplace domain variable.
# parameters: A ParameterStruct object containing the relevant model parameters.
function LaplaceDomainGreensFunction(xVals,ξ_index,s,parameters)
  #Extract parameters.
  EI = parameters.EI
  ρ = parameters.ρ
  v = parameters.v
  Q0 = parameters.Q0
  M = parameters.M
  TZ = parameters.TZ
  L = TZ.L

  Δx = xVals[2]-xVals[1]

  #Set up the system of equations to solve for the Greens function in the Laplace domain.
  LHSMatrix = zeros(Complex{eltype(xVals)}, length(xVals), length(xVals))

  RHSVector = zeros(Complex{eltype(xVals)}, length(xVals))
  if(ξ_index == 1  || ξ_index == length(xVals) || ξ_index == length(xVals)-1 || ξ_index == 2) #BUG This is designed for clamped boundary conditions.
    #RHSVector[ξ_index] = 2/(EI*Δx) #Impulse at position ξ at boundary. #TODO: How should I actually handle this?
  else
    RHSVector[ξ_index] = 1/(EI*Δx) #Impulse at position ξ.
  end

  #Fill in LHSMatrix.
  for (xi, x) in enumerate(xVals[3:end-2])
    k = tz.k_interpolant(x)
    C = tz.C_interpolant(x)

    LHSMatrix[xi, xi-2:xi+2] = [ 1/(Δx^4) , -4/(Δx^4) , 6/(Δx^4)  + ρ*s^2 + C*s + k, -4/(Δx^4) , 1/(Δx^4)  ]

  end


  #BUG: Clamped boundary conditions. Pretty sure this is correct? Does Laplace mess with this?
  LHSMatrix[1,1] = 1 #w(0) = 0
  LHSMatrix[end,end] = 1 #w(L) = 0
  LHSMatrix[2,1] = 1 #w'(0) = 0
  LHSMatrix[end-1,end] = 1 #w'(L) = 0

  #Solve the system of equations.
  ĝ = LHSMatrix\RHSVector
  return ĝ


end



#TODO: Implement.
#Calculates the time domain inversion of the Laplace domain Green’s function for the beam at position x and time t due to a unit impulse at time τ.
#Inputs:
# x: The spatial position at which to calculate the Green’s function.
# t: The time at which to calculate the Green’s function.
# ξ: The impulse position.
# parameters: A ParameterStruct object containing the relevant model parameters.
function InvertedGreensFunction(xVals, ξVals, t, parameters)
  #Extract parameters.
  EI = parameters.EI
  ρ = parameters.ρ
  v = parameters.v
  Q0 = parameters.Q0
  M = parameters.M
  TZ = parameters.TZ
  L = TZ.L

  #TODO: Implement numerical Laplace inversion.
  s -> LaplaceDomainGreensFunction(xVals, ξVals, s, parameters)


  return 0
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
    Q[i+1] = NonLinearSolve(Q_n ->  ((Q_n/c_H)^(2/3) - R(Q_n)*Heaviside(R(Q_n)) ), Q[i]) #Initial guess is previous value. #BUG: 2/3 or 3/2?
  end





  return Q
end

#TODO: Implement.
function NonLinearSolve(f, initial_guess)
  return initial_guess
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
