using InverseLaplace


#Implementation of the GWR algorithm from iLaplace package
function GWRImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)


  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]

      if t==0
        deformations[ti,xi] = 0

      else


        #deformation = real( DirectQuadratureMethod(t, s-> LaplaceSpaceFunctionMovingPointForce1TZ(x,s,parameters)) )

        ft = GWR( s -> real(LaplaceSpaceFunction(x,s,parameters)), N)
        deformation = real( ft(t) )


        #Update the deformations matrix.
        deformations[ti,xi] = deformation

      end

    end
  end


  return deformations
end



#Implementation of Weeks method using the iLaplace package
function WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters)

  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]



    ft = Weeks( s -> LaplaceSpaceFunction(x,s,parameters))




    for ti in eachindex(tVals)
      t=tVals[ti]


      deformation = real( ft(t) )


      #Update the deformations matrix.
      deformations[ti,xi] = deformation

    end
  end


  return deformations
end

#Implementation of Weeks method using the iLaplace package.
#Adds the ability to specify a function N(x) which specifies the number of terms
#to use in the Laguerre sum.
function WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)
  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]

    #Different x Values require different values of N.
    #b = 0.4 #Controls drop off at edges. Drop off lower for higher b
    #N_max = 1024 #Controls maximum number of terms in approximation

    #N = ( N_max*exp(-b*(x-  (xVals[1]+xVals[end])/2  )^2)  )

    #N = Int(round(N))

    ft = Weeks( s -> LaplaceSpaceFunction(x,s,parameters), Int(N(x)))




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

#Implementation of an unaccelerated direct quad. with ability to specify quadrature details.
function directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)
  deformations = zeros(length(tVals), length(xVals))

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]


      deformation = real( DirectQuadratureMethod(t, s-> LaplaceSpaceFunction(x,s,parameters), laplaceParameters) )


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

#Inverts the Laplace transform using a basic quadrature method with ability to specify
#quadrature details.
function DirectQuadratureMethod(t, LaplaceFunction, laplaceParameters)
  sRadius, sNum = laplaceParameters
  sStep = 2*sRadius/sNum
  sVals = 0.01*ones(sNum) + LinRange(-sRadius,sRadius,sNum)*1im

  f = 0
  for s in sVals
    f += 1/(2*pi)*LaplaceFunction(s)*exp( s*t  )*sStep #Should there be an i here?
  end

  return f
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
