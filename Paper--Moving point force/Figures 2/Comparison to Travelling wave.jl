#=
Produces comparisons to travelling wave for various methods.
=#

include("Parameter setup.jl")


#MARK: Naive implementations

#|||||||||||||||||||||||\
#||| Direct quadrature|||\
#|||||||||||||||||||||||||\
graphName = "Comp_DQ"

leftTZ = TransitionZone(0.5, k_default, k_default, C_default, C_default)
xtz_list = [leftTZ]

parameters = [EI, m, xp, xtz_list, P, v_slow]


laplaceParameters = [50,100]
inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)


if !LoadFlag
  deformations = @time CalcDynamicDeformation(xVals_slow, tVals10_slow, parameters, inversionMethod)
  if SaveFlag
    save(DataPath*"/"*graphName*".jld", "deformations", deformations)
  end
else
  deformations = load(DataPath*"/"*graphName*".jld")["deformations"]
end


SteadyDeformation1 = SteadyStateTravellingSolution(xVals_slow,tVals10_slow,parameters, k_default, C_default)
SteadyStateDeformations = [SteadyDeformation1]

gr()
Graph10Times(deformations,xVals_slow,tVals10_slow,parameters,SteadyStateDeformations, [-0.0006, 0.0002], colours, "", graphName)
#|||||||||||||||||||||||||/
#||||||||||||||||||||||||/





#||||||||||||||||||\
#||| Weeks method|||\
#||||||||||||||||||||\





#||||||||||||||||||||/
#|||||||||||||||||||/
