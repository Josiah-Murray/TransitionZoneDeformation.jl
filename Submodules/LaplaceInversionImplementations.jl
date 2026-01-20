module LaplaceInversionImplementations


"""
    GWRImplementation(xVals, tVals, parameters, LaplaceFunction; M = 20, shift_parameter = 0)

Return an approximation to the inversion of the Laplace-domain function `LaplaceFunction` using the Gaver-Wynn Rho algorithm at each time in `tVals`.
`LaplaceFunction` should be a function of a single complex variable `s` which returns an array of the same length as `xVals`, representing the Laplace-domain solution at each spatial location in `xVals`.
The optional parameters `M` and `shift_parameter` correspond to the parameters of the Gaver-Wynn Rho algorithm.
"""
function GWRImplementation(xVals, tVals, parameters, LaplaceFunction; M = 20, shift_parameter = 0)
  deformations = zeros(typeof(tVals[1]), length(xVals), length(tVals))
  for (ti, t) in enumerate(tVals)
    deformations[:,ti] = real(NILaplace.GaverWynnRho.GWR_array(LaplaceFunction, t, M, shift_parameter = shift_parameter) )
  end

  return deformations
end

"""
    GWRImplementation(; M = 20, shift_parameter = 0)

Returns a function that behaves like `GWRImplementation` with the parameters `M` and `shift_parameter` prespecified.
"""
function GWRImplementation(; M = 20, shift_parameter = 0)
  return (xVals, tVals, parameters, LaplaceFunction) -> GWRImplementation(xVals, tVals, parameters, LaplaceFunction; M = M, shift_parameter = shift_parameter)
end


"""
    WeeksImplementation(xVals, tVals, parameters, LaplaceFunction; M = 20, σ = 1/2, b = 1)

Return an approximation to the inversion of the Laplace-domain function `LaplaceFunction` using Weeks' method at each time in `tVals`.
`LaplaceFunction` should be a function of a single complex variable `s` which returns an array of the same length as `xVals`, representing the Laplace-domain solution at each spatial location in `xVals`.
The optional parameters `M`, `σ`, and `b` correspond to the parameters of Weeks' method.
"""
function WeeksImplementation(xVals, tVals, parameters, LaplaceFunction; M = 20, σ = 1/2, b = 1)

  deformation_WeeksStruct = NILaplace.Weeks.GenerateWeeksApproximation(
    LaplaceFunction,
    M,
    σ,
    b
  )
  deformations = zeros(typeof(tVals[1]), length(xVals), length(tVals))
  for (ti, t) in enumerate(tVals)
    deformations[:,ti] = real(NILaplace.Weeks.EvalWeeks(deformation_WeeksStruct, t))
  end

  return deformations
end

"""
    WeeksImplementation(; M = 20, σ = 1/2, b = 1)

Returns a function that behaves like `WeeksImplementation` with the parameters `M`, `σ`, and `b` prespecified.
"""
function WeeksImplementation(;M, σ, b)
  return (xVals, tVals, parameters, LaplaceFunction) -> WeeksImplementation(xVals, tVals, parameters, LaplaceFunction; M = M, σ = σ, b = b)
end








end
