module SteadyStateWithTransitionZones
#include(joinpath(@__DIR__, "RailDataStructures.jl"))

using ..RailDataStructures
using ..Utilities


"""
    LaplaceDomainFunction(xVals, s, xp, v, beamParameters)

Return the Laplace domain function for an Euler-Bernoulli beam on a viscoelastic foundation with parameters `beamParameters` (see #TODO: reference), for a Laplace domain variable `s`, and each x value in `xVals`.
The beam is forced by a point force starting at position `xp` and moving with velocity `v`. This function is primarily intended to be passed as an function of the form `s -> LaplaceDomainFunction(xVals, s, beamParameters)` into a Laplace inversion routine. See #TODO: reference to Laplace inversion code.
"""
function LaplaceDomainFunction(xVals, s, xp, v, beamParameters::RailDataStructures.RailParameters)

  Dt = real(typeof(s))
  xVals  = convert.(Dt, xVals)



  tzList = beamParameters.transitionZones

  #||--Find b values--||#
  #Solve for the undetermined coefficients by solving a linear system
  bVals = Coeff_solver(s, xp, v, beamParameters)
  ŷ = zeros(Complex{Dt}, length(xVals))

  segment = 1 #For keeping track of which segment we are in, so that correct coefficient values can be pulled out.
  pointForceAdded = false #We 'insert' an extra segment to account for the point force location.
  for (xi, x) in enumerate(xVals)

    #Set the correct parameters.
    #Default to the furthest right parameters, override if we are to the left.
    k = tzList[end].k_right
    C = tzList[end].C_right
    for tz in tzList
      if real(x) < tz.location
        k = tz.k_left
        C = tz.C_left
        break
      end
    end


    #||--Find r values--||#
    #Only used for graphing appropriately
    rVals = RValues(s, beamParameters.EI, beamParameters.m, C, k)
    r1, r2, r3, r4 = rVals




    #TODO: Can probably combine with code for k and C above.
    #Find correct segment to isolate the coefficients we need for this x value
    for tz in tzList
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


    #TODO: Implement option for zero velocity.
    ŷ[xi] = b1*exp(r1*x) + b2*exp(r2*x) + b3*exp(r3*x) + b4*exp(r4*x) + ICResponse(x, s, v, beamParameters, C, k) +   (1/abs(v))*(exp(-s*(x-xp)/v)  /  (  ( (beamParameters.EI*s^4)/v^4 ) +beamParameters.m*s^2 + C*s + k )  )*Utilities.Heaviside((x-xp)/v)
  end

  return ŷ
end

"""
    RValues(s, EI, m, C, k)

Calculate the characteristic roots for the Laplace-domain Euler-Bernoulli beam on a viscoelastic foundation with bending stiffness `EI`, mass per unit length `m`, foundation damping `C`, and track modulus `k`, for a Laplace domain variable `s`.

This is an internal function and not intended to be called directly by the user.
"""
function RValues(s, EI, m, C,k)#TODO Update root finding to use PolynomialRoots.jl and then a sorting function.
  Dt = real(typeof(s))
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

"""
    Coeff_solver(s, xp, v, parameters::RailDataStructures.RailParameters)

Return the unknown coefficents for the Laplace-domain Euler-Bernoulli beam on a viscoelastic foundation with transition zones and for a Laplace domain variable `s`.

This is an internal function and not intended to be called directly by the user.
"""
function CoeffSolver(s, xp, v, parameters::RailDataStructures.RailParameters) #TODO: Update to allow for xp to be anywhere and not enforce zero in limits through system.

  Dt = real(typeof(s))

  #Type of xtz_list specified so that autocomplete works.
  xtz_list = parameters.transitionZones
  EI, m, v = parameters.EI, parameters.m,  v

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
          RHS[2 + (segment-1)*4 + i] = -(-s/v)^(i-1)*(  ( 1/abs(v)  )/( ( (EI*s^4)/v^4  ) +m*s^2 + tz.C_left*s + tz.k_left   )   ) - ( ICResponse(xp+eps(Dt), s, parameters, tz.C_left, tz.k_left, i-1) - ICResponse(xp-eps(Dt), s, parameters, tz.C_left, tz.k_left, i-1) )
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
      RHS[2 + (segment-1)*4 + i] = (-s/v)^(i-1)*( 1*exp(-s*tz.location/v)/abs(v)  )*(  1/( ( (EI*s^4)/v^4  ) +m*s^2 + tz.C_left*s + tz.k_left   ) - 1/( ( (EI*s^4)/v^4  ) +m*s^2 + tz.C_right*s + tz.k_right   )   )*Heaviside((tz.location-xp)/v) - ( ICResponse(tz.location, s, parameters, tz.C_right, tz.k_right, i-1) - ICResponse(tz.location, s, parameters, tz.C_left, tz.k_left, i-1) )
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


"""
    ICResponse(x, s, v, parameters::RailDataStructures.RailParameters, C, k; derivativeOrder::Int = 0)

Returns the response of the Laplace-domain rail to the initial conditions at a position `x`, Laplace variable `s`, for rail parameters `parameters` (see #TODO: reference), foundation damping `C`, track modulus `k`.
Changing the optional keyword parameter `derivativeOrder` will instead return the derivative of the response of order `derivativeOrder` with respect to time.
By default, it uses a travelling wave associated with the left-most beam segment. To change this, use the optional keyword parameters `C_init` and `k_init` to set the foundation damping and track modulus of the travelling wave solution.

This is an internal function and not intended to be called directly by the user.
"""
function ICResponse(x, s, v, parameters::RailDataStructures.RailParameters, C, k; derivativeOrder::Int = 0, C_init = nothing, k_init = nothing)
  EI, m = parameters.EI, parameters.m
  Dt = real(typeof(s))

  if isnothing(C_init)
    C_init = xtz_list[1].C_left
  end
  if isnothing(k_init)
    k_init = xtz_list[1].k_left
  end



  λ = (k_init/(4*EI))^(1/4)
  β = C_init/(2*sqrt(k_init*m))
  α = v/((4*k_init*EI/(m^2))^(1/4))
  β_cr = (1/α)*sqrt(convert(Dt, 2//27))*(2*α^2 + sqrt(3+α^4))*(-α^2+sqrt(3+α^4))





  responseCoeff(w) = 1/(EI*w^4 + m*s^2 + C*s + k)

  η_Coefficients = [- α^2*β^2, 0, (α^4 - 1), 0,  2*α^2, 0, 1]
  η_list = roots(η_Coefficients)

  #Find η in first segment
  η = zero(Dt)+Inf*1im; #Preallocate to have something to compare to in first round of for loop

  for i in eachindex(η_list)
      if abs(imag(η_list[i]))<abs(imag(η)) && real(η_list[i]) > 0
              η = η_list[i]
      else
    end#else
  end#for

  #||--Responses to initial deformation--||#
  wMinus(θ) = (λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + 0.5 * (α*β/η)^2 ))*
    ( -( (α*β/η + η^2) )*(
        (
          responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^derivativeOrder)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ ) -
          responseCoeff((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^derivativeOrder)*exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ )
        ) /( 2im*η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) +
      (
        responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^derivativeOrder)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  ) +
        responseCoeff((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*(((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^derivativeOrder)*exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  )
      )/2
    )



  wPlus(θ) = (λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*
    ( -( (α*β/η - η^2) )*(
        (
          responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^derivativeOrder*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ ) -
          responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^derivativeOrder*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ )
        ) /( 2im*η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) ) ) +
      (
        responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^derivativeOrder*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  ) +
        responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^derivativeOrder*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  )
      )/2
    )



  wMinusCritical(θ) = (λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*(  -( (α*β/η + η^2) )
    *(
      responseCoeff( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)* ( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^derivativeOrder * exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
      -responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)* ((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^derivativeOrder * exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
    )/  ( 2*η*sqrt(abs(2*α^2 + η^2 - 2*(α*β/η) )) )
    + (
      responseCoeff((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^derivativeOrder* exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
      + responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^derivativeOrder* exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )
    )/2)


  #||--Response to initial velocity--||#
  wDotMinus(θ) = -v*(λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*
    ( -( (α*β/η + η^2) )*(
      (
        responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ ) -
        responseCoeff((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (η - 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ )
      ) /( 2im*η*sqrt(2*α^2 + η^2 - 2*(α*β/η) ) ) ) +
      (
        responseCoeff((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (η + 1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  ) +
        responseCoeff((η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)*((η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (η -1im*sqrt(2*α^2 + η^2 - 2*(α*β/η) ))*λ*θ  )
      )/2
    )


  wDotPlus(θ) = -v*(λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*
  ( -( (α*β/η - η^2) )*(
      (
        responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ ) -
        responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ )
      ) /( 2im*η*sqrt(2*α^2 + η^2 + 2*(α*β/η) ) ) ) +
    (
      responseCoeff((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (-η + im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  ) +
      responseCoeff((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)*((-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ)^(derivativeOrder+1)*exp( (-η - im*sqrt(2*α^2 + η^2 + 2*(α*β/η) ))*λ*θ  )
    )/2
  )


  wDotMinusCritical(θ) = -v*(λ/(2*k_init))*(  η /  (η^4 + α^2 * η^2 + convert(Dt, 1//2) * (α*β/η)^2 ))*( -( (α*β/η + η^2) )
    *(
      responseCoeff((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(derivativeOrder+1) * exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) -
      responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(derivativeOrder+1) * exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  )    )/  ( 2*η*sqrt(abs(2*α^2 + η^2 - 2*(α*β/η) )) )
    + (
      responseCoeff((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(derivativeOrder+1)* exp( (η + sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) +
      responseCoeff((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)*((η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ)^(derivativeOrder+1)* exp( (η - sqrt(abs(2*α^2 + η^2 - 2*(α*β/η))   ))*λ*θ  ) )/2
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



end
