using InverseLaplace


#Implementation of the GWR algorithm from iLaplace package
function GaverStehfestImplementation(xVals,tVals,LaplaceSpaceFunction, parameters)

  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]


      #deformation = real( DirectQuadratureMethod(t, s-> LaplaceSpaceFunctionMovingPointForce1TZ(x,s,parameters)) )

      ft = GWR( s -> real(LaplaceSpaceFunction(x,s,parameters)), 8 )
      deformation = real( ft(t) )


      #Update the deformations matrix.
      deformations[ti,xi] = deformation

    end
  end


  return deformations
end



#Implementation of Weeks method using the iLaplace package
function WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters)
  deformations = zeros(length(tVals), length(xVals))

  sigma0 = 0 #Right-most pole

  #Weeks' parameters
  #sigma = 0.01 + sigma0
  #b = 2*(sigma-sigma0)
  #N=1024


  eps0 = 0.0002;
  epsT = 0.001;
  maxTime = tVals[end]

  sigma = sigma0 + (1/maxTime)*log(epsT/eps0)
  b = 2*(sigma-sigma0)


  for xi in eachindex(xVals)
    x=xVals[xi]

    #Different x Values require different values of N.
    b = 0.4 #Controls drop off at edges. Drop off lower for higher b
    N_max = 1024 #Controls maximum number of terms in approximation

    N = ( N_max*exp(-b*(x-  (xVals[1]+xVals[end])/2  )^2)  )

    N = Int(round(N))

    ft = Weeks( s -> LaplaceSpaceFunction(x,s,parameters), Int(N))




    for ti in eachindex(tVals)
      t=tVals[ti]


      deformation = real( ft(t) )


      #Update the deformations matrix.
      deformations[ti,xi] = deformation

    end
  end


  return deformations
end


#Implementation of an unaccelerated direct quad. method
function directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters)
  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]


      deformation = real( DirectQuadratureMethod(t, s-> LaplaceSpaceFunction(x,s,parameters)) )


      #Update the deformations matrix.
      deformations[ti,xi] = deformation

    end
  end


  return deformations
end

#Inverts the Laplace transform using a basic quadrature method.
function DirectQuadratureMethod(t, LaplaceFunction)
  sRadius = 50
  sNum = 100
  sStep = 2*sRadius/sNum
  sVals = 0.01*ones(sNum) + LinRange(-sRadius,sRadius,sNum)*1im

  f = 0
  for s in sVals
    f += 1/(2*pi)*LaplaceFunction(s)*exp( s*t  )*sStep #Should there be an i here?
  end

  return f
end
