using InverseLaplace


#For creating a version of the coefficient solver that is 'memoized',
#i.e. it stores calculated coefficients so it only has to perform certain calculations once.
#See "Hands on design patterns and best practices with Julia" by Tom Kwong
function memoize(f)
  memo = Dict()
  return (x, parameters) -> begin
    x_hash = hash(x)
    if haskey(memo, x_hash)
      return memo[x_hash]
    else
      value = f(x, parameters)
      memo[x_hash] = value
      return value
    end
  end
end

#Implementation of the GWR algorithm from iLaplace package with memoisation of the undetermined coefficients.
function GWRImplementation_memo(xVals,tVals,LaplaceSpaceFunction, parameters, N)


  deformations = zeros(length(tVals), length(xVals))

  #TODO: I don't particularly like that this means this file and Functions - MPF_ArbitraryTZs.jl are interdependent.
  #Perhaps extract CoefficientSolverMovingPointForce1TZ to its own file? or maybe not necessary.
  coeff_memo = memoize(CoefficientSolverMovingPointForce1TZ)

  tzero_warning = false #So you only get the related warning once

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]

      if t==0
        if ~tzero_warning
          @warn "GWR cannot be performed for t=0. Defaulting NaN"
          tzero_warning = true
        end
        deformations[ti,xi] = NaN
      else


        #deformation = real( DirectQuadratureMethod(t, s-> LaplaceSpaceFunctionMovingPointForce1TZ(x,s,parameters)) )


        ft = GWR( s -> real(LaplaceSpaceFunction(x,s,parameters; Coeff_solver = coeff_memo)), N)
        deformation = real( ft(t) )


        #Note: this is not defined in this repository yet.
        #deformation = my_gwr(s -> LaplaceSpaceFunction(x,s,parameters), t, N)

        #Update the deformations matrix.
        deformations[ti,xi] = deformation
      end

    end
  end


  return deformations
end

#Implementation of the GWR algorithm from iLaplace package
function GWRImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)


  deformations = zeros(length(tVals), length(xVals))
  tzero_warning = false

  for xi in eachindex(xVals)
    x=xVals[xi]
    for ti in eachindex(tVals)
      t=tVals[ti]

      if t==0
        if ~tzero_warning
          @warn "GWR cannot be performed for t=0. Defaulting NaN"
          tzero_warning = true
        end
        deformations[ti,xi] = NaN

      else


        #deformation = real( DirectQuadratureMethod(t, s-> LaplaceSpaceFunctionMovingPointForce1TZ(x,s,parameters)) )


        ft = GWR( s -> real(LaplaceSpaceFunction(x,s,parameters)), N)
        deformation = real( ft(t) )


        #Note: this is not defined in this repository yet.
        #deformation = my_gwr(s -> LaplaceSpaceFunction(x,s,parameters), t, N)

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
