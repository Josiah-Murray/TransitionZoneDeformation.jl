module SteadyStateWithTransitionZones
#include(joinpath(@__DIR__, "RailDataStructures.jl"))

using ..RailDataStructures
using ..PolynomialRoots:roots



"""
    CalculateDeformation(xVals, x_p, beamParameters)

Calculate the deformation of an Euler-Bernoulli beam, with parameters `beamParameters` (see #TODO: reference), at each of a list of points, `xVals`, and subject to a point force with unit strength at `x_p`.
If a point force with strength P is desired, then the output can be simply multiplied by P due to the linearity of the system.
A similar argument can be applied for multiple point forces.
"""
function CalculateDeformation(xVals, x_p, beamParameters::RailDataStructures.RailParameters)

  mType = typeof(xVals[1])
  mType = mType <: Complex ? mType : Complex{mType}
  EI = beamParameters.EI
  #m = beamParameters.m
  tzList = beamParameters.transitionZones

  x_p_index = 0 #Stores the index of x_p in tz_list
  for (tzi, tz) in enumerate(tzList)
    if x_p < tz.position
      x_p_index = tzi
      k = tz.k_left
      C = tz.C_left
      insert!(tzList, tzi, RailDataStructures.TransitionZone(position = x_p, k_left = k, k_right = k, C_left = C, C_right = C))
      break
    elseif x_p == tz.position
      x_p_index = tzi
      break
    end
  end

  if x_p_index == 0
    k = tzList[end].k_right
    C = tzList[end].C_right
    append!(tzList, RailwayDeformation.TransitionZone(position = x_p, k_left = k, k_right = k, C_left = C, C_right = C))
  end

  #----Get characteristic roots for each segment and remove growing terms in left- and right-most segments.
  angleMeasure = cNum -> mod(angle(cNum)+pi/2, 2*pi)#For sorting the roots anti-clockwise from -pi/2

  LHS = zeros(mType, 4*length(tzList), 4*length(tzList))
  RHS = zeros(mType, 4*length(tzList))
  #Construct linear system for coefficients of the solution in each segment
  for (tzi, tz) in enumerate(tzList)
    if tzi == 1
      roots([tz.k_left, 0, 0, 0, EI])
      leftRoots  = sort(roots([tz.k_left ,0,0,0,EI]), by = angleMeasure)
      rightRoots = sort(roots([tz.k_right,0,0,0,EI]), by = angleMeasure)
      for i = 1:4
        LHS[i, 1:6] = [[  -leftRoots[j]^(i-1)*exp(leftRoots[j]*tz.position)   for j in 1:2]; [rightRoots[j]^(i-1)*exp(rightRoots[j]*tz.position)   for j in 1:4] ]
      end
    elseif tzi == length(tzList)
      for i in 1:4
        leftRoots  = sort(roots([tz.k_left ,0,0,0,EI]), by = angleMeasure)
        rightRoots = sort(roots([tz.k_right,0,0,0,EI]), by = angleMeasure)
        LHS[end-4+i,   end-5:end] = [ [  -leftRoots[j]^(i-1)*exp(leftRoots[j]*tz.position)   for j in 1:4]; [  rightRoots[j]^(i-1)*exp(rightRoots[j]*tz.position)   for j in 3:4] ]
      end
    else
      for i in 1:4
        leftRoots  = sort(roots([tz.k_left ,0,0,0,EI]), by = angleMeasure)
        rightRoots = sort(roots([tz.k_right,0,0,0,EI]), by = angleMeasure)
        LHS[4*(tzi-1)+i,   2+4*(tzi-2)+1:2+4*(tzi)] = [[ -leftRoots[j]^(i-1)*exp(leftRoots[j]*tz.position)   for j in 1:4]; [  rightRoots[j]^(i-1)*exp(rightRoots[j]*tz.position)   for j in 1:4] ]
      end
    end

    if x_p == tz.position
      RHS[4*(tzi)] = 1/EI
    end
  end



  coefficients = LHS\RHS


  #----Calculate deformation at each xVal
  deformation = zeros(mType, length(xVals))
  segmentIndex = 1
  charRoots = sort(roots([tzList[1].k_left ,0,0,0,EI]), by = angleMeasure)
  for (xi, x) in enumerate(xVals) #Assumes that xVals is are in ascending order (say from a call to LinRange).

    while segmentIndex <= length(tzList) && x > tzList[segmentIndex].position
      println(segmentIndex)#BUG remove
      segmentIndex += 1
    end


    if segmentIndex > length(tzList)
      charRoots = sort(roots([tzList[segmentIndex-1].k_right ,0,0,0,EI]), by = angleMeasure)
    else
      charRoots = sort(roots([tzList[segmentIndex].k_left ,0,0,0,EI]), by = angleMeasure)
    end



    if segmentIndex == 1
      deformation[xi] = sum([coefficients[j]*exp(charRoots[j]*x) for j in 1:2])
    elseif segmentIndex > length(tzList)
      deformation[xi] = sum([coefficients[end-j+1]*exp(charRoots[5-j]*x) for j in 1:2])
    else
      deformation[xi] = sum([coefficients[2+4*(segmentIndex-2)+j]*exp(charRoots[j]*x) for j in 1:4])
    end
  end
  return real(deformation)

end



end
