#include("Functions - Quadrature implementations.jl")


#|||||||||||||||#
#||--READ ME--||#
#|||||||||||||||#

#=
This file contains functions used for solving for the dynamic deformation of an
infinite Euler-Bernoulli beam on a piecewise constant viscoelastic foundation
subject to a moving load using the method of undetermined coefficients.

The function  'CalcDynamicDeformation' is the one that is intended to be called by a user.
=#


#|||||||||||||||||||||||||||||#
#||--Concerning parameters--||#
#|||||||||||||||||||||||||||||#

#MARK: TZ struct

#Struct which holds information about the transition zones.
#All fields should contain be numerical types (i.e. Int, Float64, BigFloat, ArbFloat, etc.).
struct TransitionZone
  location
  k_left
  k_right
  C_left
  C_right
end

#A function which creates a new transition zone object such that
#its left parameters match with the right parameters of a previous transition zone.
#Mostly for user convenience.
function AddConsistentTZ(location, previousTZ::TransitionZone, k_right, C_right)
  return TransitionZone(location, previousTZ.k_right, k_right, previousTZ.C_right, C_right)
end

#||||||||||||||||||||||||||||||||||||||||||||
#MARK: Call this function
#=
Calculates the deformation of the beam in the model described above.

xVals => A list of x values (e.g. LinRange(-1,1,100))

tVals => A list of time values (e.g. LinRange(0,1,100)) Note that some inversion methods do not support

parameters => A list containing information about the system being modelled.
  In particular parameters = (EI, m, xp, xtz_list, P, v)
  Where EI is the beam's elastic modulus multiplied by the moment of inertia,
  m is the mass of te beam, xp is the x coordinate of the point force at t=0,
  xtz_list is a list of transition zones (with type 'TransitionZone'),
  P is the strength of the point force, and v is the velocity of the point force.

Inversion method => One of the functions from  "Functions - Inversion schemes.jl" aliased with the relevant
  method inputs pre-inputted
  (e.g. inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals, tVals, LaplaceSpaceFunction, parameters, 40))

=#
function CalcDynamicDeformation(xVals,tVals, parameters, inversionMethod)
  if parameters[end] == 0
    #TODO Implement v=0?
    @error "CalcDynamicDeformation unimplemented for v=0"
    return NaN
  end
  deformations = inversionMethod(xVals, tVals, LaplaceSpaceFunctionMovingPointForce1TZ, parameters)
  return deformations
end
#||||||||||||||||||||||||||||||||||||||||||||



#MARK: Laplace-space
#The solution for the deformation in Laplace space.
#Defined here so that it can be inverted elsewhere.
function LaplaceSpaceFunctionMovingPointForce1TZ(x, s, parameters; Coeff_solver = CoefficientSolverMovingPointForce1TZ)

  Dt = typeof(s)
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
    if x < tz.location
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
      if x < xp
        break
      else
        segment += 1
      end
    end
    if x < tz.location
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



#MARK: IC response

#Takes in a position x, Laplace value s, set of parameters, C value, k value, and an order
#Returns the value of the Laplace domain response to the initial conditions with the given variables.
#If N>0, it returns the value of the nth x-derivative of the response.
function ICResponse(x, s, parameters, C, k, N)
  EI, m, xp, xtz_list, P, v = parameters
  Dt = typeof(s)

  k_init = xtz_list[1].k_left
  C_init = xtz_list[1].C_left


  λ = (k_init/(4*EI))^(1/4)
  β = C_init/(2*sqrt(k_init*m))
  α = v/((4*k_init*EI/(m^2))^(1/4))
  β_cr = (1/α)*sqrt(convert(Dt, 2//27))*(2*α^2 + sqrt(3+α^4))*(-α^2+sqrt(3+α^4))





  responseCoeff(w) = 1/(EI*w^4 + m*s^2 + C*s + k)

  η_Coefficients = [- α^2*β^2, 0, (α^4 - 1), 0,  2*α^2, 0, 1]
  η_list = roots(η_Coefficients)

  #Find η in first segment
  η = zero(Dt)+zero(Dt)*1im;
  for i in eachindex(η_list)
      if abs(imag(η_list[i]))<0.000001 && real(η_list[i]) > 0 #BUG The way I'm handling the imag part of η here may be causing issues.
              η = η_list[i]
      else
    end#else
  end#for

  #||--Responses to initial deformation--||#
  wMinus(θ) = (P*λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*
    ( -( (α*β/η + η^2) )*(
        (
          responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^N)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ ) -
          responseCoeff((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^N)*exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ )
        ) /( 2im*η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) +
      (
        responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^N)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  ) +
        responseCoeff((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^N)*exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  )
      )/2
    )



  wPlus(θ) = (P*λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*
    ( -( (α*β/η - η^2) )*(
        (
          responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^N*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ ) -
          responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^N*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ )
        ) /( 2im*η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) ) ) +
      (
        responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^N*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  ) +
        responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^N*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  )
      )/2
    )



  wMinusCritical(θ) = (P*λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*(  -( (α*β/η + η^2) )
    *(
      responseCoeff( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)* ( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^N * exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
      -responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)* ((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^N * exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
    )/  ( 2*η*sqrt(abs(2*α^2 + η^2 - 2*(α*β/η) )) )
    + (
      responseCoeff((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^N* exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
      + responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^N* exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
    )/2)


  #||--Response to initial velocity--||#
  wDotMinus(θ) = -v*(P*λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*
    ( -( (α*β/η + η^2) )*(
      (
        responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(N+1)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ ) -
        responseCoeff((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(N+1)*exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ )
      ) /( 2im*η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) +
      (
        responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(N+1)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  ) +
        responseCoeff((η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(N+1)*exp( (η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  )
      )/2
    )


  wDotPlus(θ) = -v*(P*λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*
  ( -( (α*β/η - η^2) )*(
      (
        responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(N+1)*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ ) -
        responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(N+1)*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ )
      ) /( 2im*η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) ) ) +
    (
      responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(N+1)*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  ) +
      responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(N+1)*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  )
    )/2
  )


  wDotMinusCritical(θ) = -v*(P*λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*( -( (α*β/η + η^2) )
    *(
      responseCoeff((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(N+1) * exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) -
      responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(N+1) * exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )    )/  ( 2*η*sqrt(abs(2*α^2 + η^2 - 2*(α*β/η) )) )
    + (
      responseCoeff((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(N+1)* exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) +
      responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(N+1)* exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) )/2
    )




  θ = x-xp
  if θ<0
    if(β<β_cr)
      Y₀ = wMinus(θ)
      Y₁ = wDotMinus(θ)
    else
      Y₀ = wMinusCritical(θ)
      Y₁ = wDotMinusCritical(θ)
    end
  else
    Y₀ = wPlus(θ)
    Y₁ = wDotPlus(θ)
  end



  icResponse = (m*s + C)*Y₀ + m*Y₁

  return icResponse
end

#Define the Heaviside function.
function Heaviside(x)
  Dt = typeof(x)
  if x<0
    return zero(Dt)
  elseif x==0
    return convert(Dt, 1//2)
  else
    return one(Dt)
  end
end

#MARK: r values

#Roots of the characteristic equation
function RValues(s, EI, m, C,k)
  Dt = typeof(s)
  r⁴ = -(m*s^2 + C*s + k)/EI

  #Always find rj such that it is in the jth quadrant of the complex plane
  ρ = ( abs(r⁴) )^(1/4)
  θ = mod( atan( imag(r⁴), real(r⁴) ), 2*convert(Dt, pi) )/4

  r1 = ρ*exp(θ*1im)
  r2 = r1*1im
  r3 = -r1
  r4 = -r1*1im

  rVals = [r1,r2,r3,r4]
end

#MARK: System

#Solves for the undetermined coefficients and returns the correct set of 4 depending on the value of x
function CoefficientSolverMovingPointForce1TZ(s, parameters)

  Dt = typeof(s)

  #Type of xtz_list specified so that autocomplete works.
  EI, m, xp, xtz_list, P, v = parameters

  #LHS matrix for system enforcing continuity conditions.
  #We need four rows for each continuity condition (1 for each tz and one for pf)
  #We also need four rows to enforce zero deformation at limits.
  numRows = (length(xtz_list)+2)*4


  #Initialise with first row (enforcing zero deformation at inf. on left)
  LHSMatrix = [ [0 1 0 0 ; 0 0 1 0 ]  zeros(Complex{Dt}, 2, numRows - 4) ]



  #RHS Vector for system.
  RHS = zeros(Complex{Dt}, numRows)

  pointForceAdded = false #So we can stop comparing the locations

  segment = 1


  for i in eachindex(xtz_list)
    tz::TransitionZone = xtz_list[i] #Type specified for ease of coding (i.e. autocomplete)
    if(!pointForceAdded)
      #Check if point force is to left of transition zone
      if(xp<tz.location)

        #Add row associated with continuity at point force
        r1, r2, r3, r4 = RValues(s,EI,m, tz.C_left, tz.k_left)
        continuitySubMatrix = [-exp(r1*xp)  -exp(r2*xp) -exp(r3*xp) -exp(r4*xp) exp(r1*xp)  exp(r2*xp) exp(r3*xp) exp(r4*xp)]
        for i = 1:3
          continuitySubMatrix = [continuitySubMatrix; -r1^i*1exp(r1*xp)  -r2^i*exp(r2*xp) -r3^i*exp(r3*xp) -r4^i*exp(r4*xp) r1^i*exp(r1*xp)  r2^i*exp(r2*xp) r3^i*exp(r3*xp) r4^i*exp(r4*xp) ]
        end
        LHSMatrix = [LHSMatrix;     zeros(Complex{Dt}, 4, 4*(segment-1))      continuitySubMatrix     zeros(Complex{Dt}, 4, numRows-4*(segment-1)-8)  ]

        #RHS Vector associated with point force
        for i = 1:4
          RHS[2 + (segment-1)*4 + i] = -(-s/v)^(i-1)*(  ( P/abs(v)  )/( ( (EI*s^4)/v^4  ) +m*s^2 + tz.C_left*s + tz.k_left   )   ) - ( ICResponse(xp+eps(), s, parameters, tz.C_left, tz.k_left, i-1) - ICResponse(xp-eps(), s, parameters, tz.C_left, tz.k_left, i-1) )#BUG Will epsilon be too large here?
        end


        pointForceAdded = true



        segment += 1
      end
    end

    #----Add row associated with continuity at transition zone
    lr1, lr2, lr3, lr4 = RValues(s,EI,m, tz.C_left, tz.k_left)
    rr1, rr2, rr3, rr4 = RValues(s,EI,m, tz.C_right, tz.k_right)
    continuitySubMatrix = [-exp(lr1*tz.location)  -exp(lr2*tz.location) -exp(lr3*tz.location) -exp(lr4*tz.location) exp(rr1*tz.location)  exp(rr2*tz.location) exp(rr3*tz.location) exp(rr4*tz.location)]
    for i = 1:3
      continuitySubMatrix = [continuitySubMatrix; -lr1^i*1exp(lr1*tz.location)  -lr2^i*exp(lr2*tz.location) -lr3^i*exp(lr3*tz.location) -lr4^i*exp(lr4*tz.location) rr1^i*exp(rr1*tz.location)  rr2^i*exp(rr2*tz.location) rr3^i*exp(rr3*tz.location) rr4^i*exp(rr4*tz.location) ]
    end
    LHSMatrix = [LHSMatrix;     zeros(Complex{Dt}, 4, 4*(segment-1))      continuitySubMatrix     zeros(Complex{Dt}, 4, numRows-4*(segment-1)-8)  ]

    #RHS Vector
    for i = 1:4
      RHS[2 + (segment-1)*4 + i] = (-s/v)^(i-1)*( P*exp(-s*tz.location/v)/abs(v)  )*(  1/( ( (EI*s^4)/v^4  ) +m*s^2 + tz.C_left*s + tz.k_left   ) - 1/( ( (EI*s^4)/v^4  ) +m*s^2 + tz.C_right*s + tz.k_right   )   )*Heaviside((tz.location-xp)/v) - ( ICResponse(tz.location, s, parameters, tz.C_right, tz.k_right, i-1) - ICResponse(tz.location, s, parameters, tz.C_left, tz.k_left, i-1) )
    end


    segment += 1
  end

  #add last row corresponding to zero deformation at inf. on right
  LHSMatrix = [LHSMatrix; zeros(Complex{Dt},2, numRows - 4) [ 1 0 0 0 ; 0 0 0 1] ]



  bVals = try
    LHSMatrix\RHS
  catch
    @error "LHSMatrix contains NaN/inf"
    println("s value: ", s)
  end




  bVals[2] = zero(Dt)
  bVals[3] = zero(Dt)
  bVals[end] = zero(Dt)
  bVals[end-3] = zero(Dt)







  return bVals

end






#|||||||||||||||||\\
#MARK: Steady-state
#|||||||||||||||||//

#=

The functions in this section relate to steady-state (i.e. travelling wave)
solutions for a point force traversing a beam on a homogeneous foundation.

=#



using PolynomialRoots
#This solution is found using the paper "Analytical Solutions for Euler-Bernoulli
#Beam on Viscoelastic Foundation Subjected to Moving Load".
function SteadyStateTravellingSolution(xVals, tVals, parameters, k, C)
  EI, m, xp, ~, P, v = parameters


  λ = (k/(4*EI))^(1/4)
  β = C/(2*sqrt(k*m))
  α = v/((4*k*EI/(m^2))^(1/4))
  β_cr = (1/α)*sqrt(2/27)*(2*α^2 + sqrt(3+α^4))*(-α^2+sqrt(3+α^4))

  η_Coefficients = [- α^2*β^2, 0, (α^4 - 1), 0,  2*α^2, 0, 1]
  η_list = roots(η_Coefficients)

  #Find η in first segment
  η = 0+0im;
  for i in eachindex(η_list)
      if abs(imag(η_list[i]))<0.000001 && real(η_list[i]) > 0
              η = η_list[i]
      else
    end#else
  end#for

  wMinus(θ) = (P*λ/(2*k))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*
  ( -( (α*β/η + η^2) )*(
  (exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ ) -
  exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ )
  ) /( 2im*η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) +
  (exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  ) +
  exp( (η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  )
  )/2
  )


  wPlus(θ) = (P*λ/(2*k))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*
  ( -( (α*β/η - η^2) )*(
  (exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ ) -
  exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ )
  ) /( 2im*η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) ) ) +
  (exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  ) +
  exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  )
  )/2
  )


  wMinusCritical(θ) = (P*λ/(2*k))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*( -( (α*β/η + η^2) )*(  exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) - exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )    )/  ( 2*η*sqrt(abs(2*α^2 + η^2 - 2*(α*β/η) )) )  + ( exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) + exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) )/2)

  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]
      θ = x-xp-v*t
      if θ<0
        if(β<β_cr)
          deformations[ti,xi] = real(wMinus(θ))
        else
          deformations[ti,xi] = real(wMinusCritical(θ))
        end

      else
        deformations[ti,xi] = real(wPlus(θ))
      end

    end
  end
  return deformations
end

#This solution is found using the paper "Analytical Solutions for Euler-Bernoulli
#Beam on Viscoelastic Foundation Subjected to Moving Load".
function DerivativeSteadyStateTravellingSolution(xVals, tVals, parameters, ks, C)
  EI, m, xp, xtz_list, P, v = parameters


  λ = (ks/(4*EI))^(1/4)
  β = C/(2*sqrt(ks*m))
  α = v/((4*ks*EI/(m^2))^(1/4))
  β_cr = (1/α)*sqrt(2/27)*(2*α^2 + sqrt(3+α^4))*(-α^2+sqrt(3+α^4))

  η_Coefficients = [- α^2*β^2, 0, (α^4 - 1), 0,  2*α^2, 0, 1]
  η_list = roots(η_Coefficients)

  #Find η in first segment
  η = 0+0im;
  for i in eachindex(η_list)
      if abs(imag(η_list[i]))<0.000001 && real(η_list[i]) > 0
              η = η_list[i]
      else
    end#else
  end#for

  wMinusDash(θ) = (P*λ/(2*ks))*(  ( -v*η^2*λ*exp( η* λ*θ  ) ) /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*( -( (α*β/η + η^2) )*(sin(sqrt(2*α^2 + η^2 - 2*(α*β/η) )*λ*θ ) /  ( η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) + cos( sqrt(2*α^2 + η^2 - 2*(α*β/η) )*λ*θ  )) -v*sqrt(2*α^2 + η^2 - 2*(α*β/η) )*λ*(P*λ/(2*ks))*(  ( η*exp( η* λ*θ  ) ) /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*( -( (α*β/η + η^2) )*(cos(sqrt(2*α^2 + η^2 - 2*(α*β/η) )*λ*θ ) /  ( η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) - sin( sqrt(2*α^2 + η^2 - 2*(α*β/η) )*λ*θ  ))

  wPlusDash(θ) = (P*λ/(2*ks))*(  ( v*η^2*λ*exp( -η* λ*θ  ) ) /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*( -( (α*β/η - η^2) )*(sin(sqrt(2*α^2 + η^2 + 2*(α*β/η) )*λ*θ ) /  ( η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) ) ) + cos( sqrt(2*α^2 + η^2 + 2*(α*β/η) )*λ*θ  )) -  v*sqrt(2*α^2 + η^2 + 2*(α*β/η) )*λ*  (P*λ/(2*ks))*(  ( η*exp( -η* λ*θ  ) ) /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*( -( (α*β/η - η^2) )*(cos(sqrt(2*α^2 + η^2 + 2*(α*β/η) )*λ*θ ) /  ( η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) ) ) - sin( sqrt(2*α^2 + η^2 + 2*(α*β/η) )*λ*θ  ))


  wMinusCriticalDash(θ) = (P*λ/(2*ks))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*( -( (α*β/η + η^2) )*( -v*(η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) + v*(η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )    )/  ( 2*η*sqrt(abs(2*α^2 + η^2 - 2*(α*β/η) )) )  + ( -v*(η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) -v*(η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ* exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) )/2)



  wDotMinus(θ) = -v*(P*λ/(2*k))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*
    ( -( (α*β/η + η^2) )*(
        (
         ((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ ) -
          ((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ )
        ) /( 2im*η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) +
      (
        ((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  ) +
        ((η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*exp( (η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  )
      )/2
    )


  wDotPlus(θ) = -v*(P*λ/(2*k))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*
    ( -( (α*β/η - η^2) )*(
        (
          ((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ ) -
          ((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ )
        ) /( 2im*η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) )
      ) +
      (
        ((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  ) +
        ((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  )
      )/2
    )


  wDotMinusCritical(θ) = -v*(P*λ/(2*k))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*( -( (α*β/η + η^2) )
    *(    ( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ) * exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
    - ((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ) * exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )    )/  ( 2*η*sqrt(abs(2*α^2 + η^2 - 2*(α*β/η) )) )
    + ( ((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)* exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
    + ((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)* exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) )/2)



  deformations = zeros(length(tVals), length(xVals))
  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]
      θ = x-xp-v*t
      if θ<0
        if(β<β_cr)
          deformations[ti,xi] = real(wMinusDash(θ))
        else
          deformations[ti,xi] = real(wMinusCriticalDash(θ))
        end

      else
        deformations[ti,xi] = real(wPlusDash(θ))
      end

    end
  end
  return deformations
end





#Returns the critical velocity for a constant point force moving
#with a fixed velocity on an purely elastic (i.e. undamped) foundation.
#Kenney, J. T. (1954).
#Steady-State Vibrations of Beam on Elastic Foundation for Moving Load.
#Journal of Applied Mechanics, 21(4), 359–364.
#https://doi.org/10.1115/1.4010934
function CriticalVelocity(EI, m , k)
  v_cr = (4*k*EI/(m^2))^(1/4)
  return v_cr
end



#These two functions implement the damping ratio calculations described in
#Basu, D., & Rao, N. S. V. K. (2013).
#Analytical solutions for Euler–Bernoulli beam on
# visco‐elastic foundation subjected to moving load.
#International Journal for Numerical and Analytical Methods in Geomechanics, 37(8), 945–960.
#https://doi.org/10.1002/nag.1135
function CriticalDampingRatio(EI,m,k,v)
  α = v/((4*k*EI/(m^2))^(1/4))
  β_cr = (1/α)*sqrt(2/27)*(2*α^2 + sqrt(3+α^4))*sqrt(-α^2+sqrt(3+α^4))
  return β_cr
end
function DampingRatio(m,C,k)
  return C/(2*sqrt(k*m))
end
