#|||||||||||||||#
#||--READ ME--||#
#|||||||||||||||#
#=
This file contains functions used for solving for the dynamic deformation of an
infinite Euler-Bernoulli beam on a piecewise constant viscoelastic foundation
subject to a moving load using the method of undetermined coefficients.

In particular, is it designed with 'MovingPointForce.jl' in mind.
=#



function CalcDynamicDeformation(xVals,tVals, parameters)

  #DEBUG: forming matrix of all the x and t values
  tempArray = zeros(Float64, length(xVals), length(tVals))
  for xi in eachindex(xVals)
    for ti in eachindex(tVals)

      tempArray[xi, ti] = xVals[xi]

    end
  end

  return "Not yet implemented."
end
