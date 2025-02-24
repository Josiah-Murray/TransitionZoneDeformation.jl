

include("Functions - MPF_ArbitraryTZs.jl")


#|||||||||||||||||#
#||--Variables--||#
#|||||||||||||||||#

#||--Model parameters--||#
#Things like the beam stiffness, etc.
const EI = 2.3*10^3 #Beam stiffness
const m = 48 #Mass per unit length of beam

const xp = 0 #Starting position of point force

const P = 10^4 #Force conveyed by point load #BUG force has negative effect for some reason. Missing negative probably.
const v = 1 #Speed of point force


#||||||||||||||||||||||||#
#||--Transition zones--||#
#||||||||||||||||||||||||#



const leftTZ = TransitionZone(2, 6.9*10^7, 2*6.9*10^7, 10^7, 2*10^7 )
#rightTZ = AddConsistentTZ(4, leftTZ, 6.9*10^7, 10^7)
#xtz_list = [leftTZ, rightTZ]
xtz_list = [leftTZ]


parameters = [EI, m, xp, xtz_list, P, v]



#||--Solution parameters--||#
#The time and space values for which I want an evaluation.
xLeft = -1
xRight = 4
xNum = 20
xVals = LinRange(xLeft,xRight,xNum)

tMin = 0
tMax = 6
tNum = 20
tVals = LinRange(tMin,tMax,tNum)


#||||||||||||||||||||||||||||#
#||--Finding the solution--||#
#||||||||||||||||||||||||||||#

#Calculate the deformation of the beam for given times and x values
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters)
println("Finished calculating - Beginning graphing:")
##


#||||||||||||||||#
#||--Graphing--||#
#||||||||||||||||#

include("Graphing - MPF_ArbitraryTZs.jl")
println("DEBUG")
