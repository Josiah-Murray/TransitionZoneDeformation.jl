module SteadyStateWithTransitionZones
#include(joinpath(@__DIR__, "RailDataStructures.jl"))

using ..RailDataStructures


"""
    LaplaceDomainFunction(xVals, s, beamParameters)

Return the Laplace domain function for an Euler-Bernoulli beam on a viscoelastic foundation with parameters `beamParameters` (see #TODO: reference), for a Laplace domain variable `s`, and each x value in `xVals`.
This function is primarily intended to be passed as an function of the form `s -> LaplaceDomainFunction(xVals, s, beamParameters)` into a Laplace inversion routine. See #TODO: reference to Laplace inversion code.
"""
function LaplaceDomainFunction(xVals, s, beamParameters::RailDataStructures.RailParameters)

  #TODO Convert to use xVals instead of single x value.
  #TODO Currently only coppied from old code: update.

  Dt = real(typeof(s))
  xVals  = convert.(Dt, xVals)

  #TODO: Do we convert?
  #We convert here so that we don't have to know ahead of time what type we'll need and can let the inversion implementation decide.
  EI, m, xp, tzList, P, v = parameters
  EI = convert(Dt, EI)
  m  = convert(Dt,  m)
  xp = convert(Dt, xp)
  P  = convert(Dt,  P)
  v  = convert(Dt,  v)
  parameters = [EI, m, xp, tzList, P, v]

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
  rVals = RValues(s, EI, m, C, k)
  r1, r2, r3, r4 = rVals


  #||--Find b values--||#
  #Solve for the undetermined coefficients by solving a linear system
  bVals = Coeff_solver(s, parameters)

  #Find correct segment to isolate the coefficients we need for this x value
  segment = 1
  pointForceAdded = false
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

  ŷ = b1*exp(r1*x) + b2*exp(r2*x) + b3*exp(r3*x) + b4*exp(r4*x) + ICResponse(x,s, parameters, C, k, 0) +   (P/abs(v))*(exp(-s*(x-xp)/v)  /  (  ( (EI*s^4)/v^4 ) +m*s^2 + C*s + k )  )*Heaviside((x-xp)/v)

  return ŷ
end

end
