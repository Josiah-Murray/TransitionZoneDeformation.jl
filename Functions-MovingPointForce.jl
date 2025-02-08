using InverseLaplace
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
      deformation = real( DirectQuadratureMethod(t, s-> LaplaceSpaceFunctionMovingPointForce1TZ(x,s,parameters)) )


      #Update the deformations matrix.
      deformations[ti,xi] = deformation

    end
  end

  return deformations
end


#The solution for the deformation in Laplace space.
#Defined here so that it can be inverted elsewhere.
function LaplaceSpaceFunctionMovingPointForce1TZ(x,s, parameters)
  EI, m, C0, C1, k0, k1, xp, xtz, P, v = parameters

  #||--Find r values--||#
  rVals = RValues(x,s, parameters)
  r1, r2, r3, r4 = rVals


  #||--Find b values--||#
  #Solve for the undetermined coefficients by solving a linear system
  b1,b2,b3,b4,b5,b6,b7,b8 = CoefficientSolverMovingPointForce1TZ(s,rVals,parameters)

  #||--Construct Laplace-space function--||#
  #Need to choose the correct Laplace-domain function
  #based on where our x value is relative to the transition zone and starting point of the point force.

  if(x<xp)
    #leftmost segment
    ŷ = b1*exp(r1*x) + b2*exp(r4*x)

  elseif (x<xtz)
    #Inner segment
    ŷ = b3*exp(r1*x) + b4*exp(r2*x) + b5*exp(r3*x) + b6*exp(r4*x) + (P/abs(v))*(exp(-s*x/v)  /  (  ( (EI*s^4)/v^4 ) +m*s^2 + C0*s + k0 )  )

  else
    #Rightmost segment
    ŷ = b7*exp(r2*x) + b8*exp(r3*x) + (P/abs(v))*(exp(-s*x/v)  /  (  ( (EI*s^4)/v^4 ) +m*s^2 + C1*s+k1 )  )

  end



  return ŷ

end

function RValues(x,s, parameters)
  EI, m, C0, C1, k0, k1, xp, xtz, P, v = parameters
  if(x<xtz)
    r⁴ = -(m*s^2 + C0*s + k0)/EI
  else
    r⁴ = -(m*s^2 + C1*s + k1)/EI
  end
  #Always find rj such that it is in the jth quadrant of the complex plane
  ρ = ( abs(r⁴) )^(1/4)
  θ = mod( atan( imag(r⁴), real(r⁴) ), 2*pi )/4

  r1 = ρ*exp(θ*1im)
  r2 = r1*1im
  r3 = -r1
  r4 = -r1*1im

  rVals = [r1,r2,r3,r4]
end


#Solves for the undetermined coefficients
function CoefficientSolverMovingPointForce1TZ(s, rVals, parameters)
  EI, m, C0, C1, k0, k1, xp, xtz, P, v = parameters

  #LHS matrix for system enforcing continuity conditions.
  LHSMatrix = zeros(Complex{Float64}, 8,8)

  #RHS Vector for system.
  RHS = zeros(Complex{Float64},8)

  leftRValues = RValues(xtz-2*eps(),s, parameters)
  rightRValues = RValues(xtz+2*eps(),s, parameters)

  #||--Setting up LHSMatrix--||#

  #First block (associated with left most beam section at xp)
  for row = 1:4
    LHSMatrix[row, 1] = -leftRValues[1]^(row-1)*exp(leftRValues[1]*xp)
    LHSMatrix[row, 2] = -leftRValues[4]^(row-1)*exp(leftRValues[4]*xp)
  end

  #Upper middle block (associated with middle beam section at the point force origin)
  for row = 1:4
    for column = 3:6
      LHSMatrix[row, column] = leftRValues[column-2]^(row-1)*exp(leftRValues[column-2]*xp)
    end
  end


  #Lower middle block (associated with the middle beam section at the transition zone)
  for row = 5:8
    for column = 3:6
      LHSMatrix[row, column] = -leftRValues[column-2]^(row-4-1)*exp(leftRValues[column-2]*xtz)
    end
  end

  #Bottom right block (associated with the rightmost beam section at the transition zone)
  for row = 5:8
    LHSMatrix[row, 7] = rightRValues[2]^(row-4-1)*exp(rightRValues[2]*xtz )
    LHSMatrix[row, 8] = rightRValues[3]^(row-4-1)*exp(rightRValues[3]*xtz )
  end


  #||--RHS vector--||#
  #Top section (associated with the point force origin)
  for row = 1:4
    RHS[row] = -(-s/v)^(row-1)*(  ( P*exp(-s*xp/v)/abs(v)  )/( ( (EI*s^4)/v^4  ) +m*s^2 + C0*s + k0   )   )
  end


  #Bottom section (associated with transition zone)
  for row = 5:8
    RHS[row] = (-s/v)^(row-4-1)*( P*exp(-s*xtz/v)/abs(v)  )*(  1/( ( (EI*s^4)/v^4  ) +m*s^2 + C0*s + k0   ) - 1/( ( (EI*s^4)/v^4  ) +m*s^2 + C1*s + k1   )   )
  end


  bVals = LHSMatrix\RHS



  return bVals

end

function TalbotMethod(tVal, LaplaceFunction)
  return "Not yet implemented." #//TODO Implement Talbot method
end

#Inverts the Laplace transform using a basic quadrature method.
function DirectQuadratureMethod(t, LaplaceFunction)
  sRadius = -10
  sNum = 10
  sStep = 2*sRadius/sNum
  sVals = 0.01*ones(sNum) + LinRange(-sRadius,sRadius,sNum)*1im

  f = 0
  for s in sVals
    f += LaplaceFunction(s)*exp( s*t  )*sStep
  end

  return f
end
