using WAV #To play a sound when computation finished

include("Functions - MPF_ArbitraryTZs.jl")


#|||||||||||||||||#
#||--Variables--||#
#|||||||||||||||||#

#||--Model parameters--||#
#Things like the beam stiffness, etc.
const EI = 2.3*10^3 #Beam stiffness
const m = 48 #Mass per unit length of beam

const xp = 0 #Starting position of point force

const P = -10^4 #Force conveyed by point load
const v = 1 #Speed of point force


#||||||||||||||||||||||||#
#||--Transition zones--||#
#||||||||||||||||||||||||#



leftTZ = TransitionZone(2, 6.9*10^7, 2*6.9*10^7, 10^7, 10^7 )
rightTZ = AddConsistentTZ(3, leftTZ, 6.9*10^7, 10^7)
#anotherRight = AddConsistentTZ(4, rightTZ, 6.9*10^7, 2*10^7)
xtz_list = [leftTZ, rightTZ]
#xtz_list = [leftTZ]


parameters = [EI, m, xp, xtz_list, P, v]



#||--Solution parameters--||#
#The time and space values for which I want an evaluation.
xLeft = -1
xRight = 5
xNum = 100
xVals = LinRange(xLeft,xRight,xNum)

tMin = 0
tMax = 5
tNum = 100
tVals = LinRange(tMin,tMax,tNum)


#||||||||||||||||||||||||||||#
#||--Finding the solution--||#
#||||||||||||||||||||||||||||#

#NOTE - To change the underlying Laplace inversion method,
#change the implementation used in 'invertLaplace()'.


#Calculate the deformation of the beam for given times and x values
deformations = @time CalcDynamicDeformation(xVals,tVals, parameters)
println("Finished calculating - Beginning graphing:")
##


#||||||||||||||||#
#||--Graphing--||#
#||||||||||||||||#

include("Graphing - MPF_ArbitraryTZs.jl")



#Play sound when computation finished.
y, fs = wavread(raw"C:/Windows/Media/Ring01.wav")
wavplay(y, fs)

println("ENDED")
