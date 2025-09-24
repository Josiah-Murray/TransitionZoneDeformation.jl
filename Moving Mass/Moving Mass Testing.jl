transitionZone = constructTransitionZone(3, x -> 1, x-> 1)

parameters = ParameterStruct(
  EI =1,
  ρ = 3.0,
  v = 1.0,
  Q0 = 1,
  M = 10,
  TZ = transitionZone,
  u_0 = 0,
  uDot_0 = 0
)

deformation = CalcDynamicDeformationMM(40, parameters)
