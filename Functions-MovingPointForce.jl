#|||||||||||||||#
#||--READ ME--||#
#|||||||||||||||#

#=
This file contains functions used for solving for the dynamic deformation of an
infinite Euler-Bernoulli beam on a piecewise constant viscoelastic foundation
subject to a moving load using the method of undetermined coefficients.

In particular, is it designed with 'MovingPointForce.jl' in mind.
=#


#Calculates the deformation of the beam in the model described above.
function CalcDynamicDeformation(xVals,tVals, parameters)


  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]


      #Solve for the deformation of the beam at each time and x value.
      deformation = TalbotMethod(tVal, s-> LaplaceSpaceFunctionMovingPointForce1TZ(x,s,parameters))


      #Update the deformations matrix.
      deformations[ti,xi] = deformation

    end
  end

  return deformations
end


#The solution for the deformation in Laplace space.
#Defined here so that it can be inverted elsewhere.
function LaplaceSpaceFunctionMovingPointForce1TZ(x,s, parameters)
  EI, m, C0, C1, k0, k1, xtz, P, v = parameters

  #||--Find r values--||#
  if(x<xtz)
    r⁴ = -(m*s^2 + C0*s + k0)/EI
  else
    r⁴ = -(m*s^2 + C1*s + k1)/EI
  end
  #Always find rj such that it is in the jth quadrant of the complex plane
  ρ = ( abs(r⁴) )^(1/4)
  theta = mod( atan( imag(r⁴), real(r⁴) ), 2*pi )/4

  r1 = rho*exp(theta*1im)
  r2 = r1*1im
  r3 = -r1
  r4 = -r1*1im

  rVals = [r1,r2,r3,r4]


  #||--Find b values--||#
  #Solve for the undetermined coefficients by solving a linear system
  b1,b2,b3,b3 = CoefficientSolverMovingPointForce1TZ(rVals,parameters)

  #||--Construct Laplace-space function--||#
  #Need to choose the correct Laplace-domain function
  #based on where our x value is relative to the transition zone.
  if(x<xtz)
    #Solve for y hat 0, i.e. left of transition zone.
    ŷ =  b1*exp(r1*x) + b4*exp(r4*x) + (P/abs(v))*(exp(-s*x/v)  /  (  ( (EI*s^4)/v^4 ) +m*s^2 + C0*s+k0 )  )
  else
    #Solve for y hat 1, i.e. Right of transition zone.
    ŷ =  b2*exp(r2*x) + b3*exp(r3*x) + (P/abs(v))*(exp(-s*x/v)  /  (  ( (EI*s^4)/v^4 ) +m*s^2 + C1*s+k1 )  )
  end

  return ŷ
end

#Solves for the undetermined coefficients
function CoefficientSolverMovingPointForce1TZ(rVals, parameters)
  EI, m, C0, C1, k0, k1, xtz, P, v = parameters

  #LHS matrix for system enforcing continuity conditions.
  LHSMatrix = zeros(4,4)

  #RHS Vector for system.
  RHS = zeros(4)

  for column = 1:4
    for row  = 1:4
      if row == 1 || row == 4
        sign = -1
      end
      LHSMatrix[row, column] = sign*rVals[row]^(column-1)*exp(rVals[row]*xtz)
    end
    RHS[column] = (-s/v)^(column-1)*( P*exp(-s*xtz/v)/abs(v)  )*(  1/( ( (EI*s^4)/v^4  ) +m*s^2 + C0*s + k0   ) - 1/( ( (EI*s^4)/v^4  ) +m*s^2 + C1*s + k1   )   )
  end

  bVals = LHSMatrix\RHS

  return bVals

end

function TalbotMethod(tVal, LaplaceFunction)
  return "Not yet implemented." #//TODO Implement Talbot method
end
