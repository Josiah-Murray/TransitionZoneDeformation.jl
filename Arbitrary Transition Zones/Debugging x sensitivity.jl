using WAV #To play a sound when computation finished
using Plots
include("Functions - Inversion schemes.jl")
include("Functions - MPF_ArbitraryTZs.jl")
include("Graphing - MPF_ArbitraryTZs.jl")


#|||||||||||||||||#
#||--Variables--||#
#|||||||||||||||||#

#||--Model parameters--||#
#Things like the beam stiffness, etc.
#For constant parameters I use the default values from
  #"Steady-state dynamic response of a Bernoulli-Euler beam on a
  #viscoelastic foundation subject to a platoon of moving loads"
#For changing variables(i.e. k, C, v), I use the ranges from the same paper.

const EI = 2.3*10^3 #Beam stiffness
const m = 48 #Mass per unit length of beam
const P = -10^4 #Force conveyed by point load



const xp = 0 #Starting position of point force



const xLeft = 0
const xRight = 3
const xNum = 100
const xVals = LinRange(xLeft,xRight,xNum)

const tMin = 0
const tMax = 2
const tNum = 100
const tVals = LinRange(tMin,tMax,tNum)

colour1 = RGB(0.19215685, 0.27843137, 0.3333333)
colour2 = RGB(0.14901960784313725, 0.6274509803921569, 0.8549019607843137)

colour3 = RGB(0, 0.27450980392156865, 0.4980392156862745)
colour4 = RGB(0.6470588235294118, 0.8, 0.5098039215686274)

redLineColour = RGB(0.5,0.2,0.15)




k1 = 6.9*10^7
k2 = 2*k1
C1 = 10^7
C2 = 2*C1
v = 1




leftTZ = TransitionZone(0.4, k1, k2, C1, C2)
tz2 = AddConsistentTZ(0.6, leftTZ, k1, C1)
tz3 = AddConsistentTZ(1.4,tz2,k2,C2)
tz4 = AddConsistentTZ(1.6, tz3, k1, C1)
xtz_list = [leftTZ,tz2,tz3,tz4]


#||--With Weeks--||#


parameters= [EI, m, xp, xtz_list, P, v]

function N(x)
  return 1024
end
#xValsTest = xVals[1:50]
#xValsTest = LinRange(0.31,0.325, 400)
xValsTest = xVals

inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> WeeksMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, N)


#laplaceParameters = [50,100]
#inversionMethod = (xVals, tVals,LaplaceSpaceFunction, parameters) -> directQuadratureMethodImplementation(xVals,tVals,LaplaceSpaceFunction, parameters, laplaceParameters)


deformations = @time CalcDynamicDeformation(xValsTest,tVals, parameters, inversionMethod)
SteadyDeformation1 = SteadyStateTravellingSolution(xValsTest,tVals,parameters, k1, C1)
SteadyDeformation2 = SteadyStateTravellingSolution(xValsTest,tVals,parameters, k2, C2)

SteadyStateDeformations = [SteadyDeformation1, SteadyDeformation2]

MultiPlotTimeEvolution(deformations, xValsTest, tVals, [1,10,20,30,40,50,60,70,80,90,100], parameters, SteadyStateDeformations, [-0.001,0.001], colour4, colour3)


#Play sound when computation finished.
y, fs = wavread(raw"C:/Windows/Media/Ring01.wav")
wavplay(y, fs)

println("ENDED")




#=
#||--Testing against finite difference--||#

#returns the central difference approximation of order 2 for the fourth derivative.
#Takes in a 1D array and a centre index i
function FDM_Centre4(array, i)
  ∂⁴_x = array[i-2] - 4*array[i-1] + 6*array[i] - 4*array[i+1] + array[i+2]
  return ∂⁴_x
end

function FDM_Forward1(array,i)
  ∂_t =(-3/2)*array[i] + 2*array[i+1]-(1\2)*array[i+2]
  return ∂_t
end

function FDM_Forward2(array, i)
  ∂²_t = 2*array[i] - 5*array[i+1] +4*array[i+2] - array[i+4]
end

function approximateFit(deformations, parameters, k, C)
  EI, m, xp, xtz_list, P, v = parameters

  xIndexes = 3:length(deformations[1,:])-3

  tIndexes = 1:length(deformations[:,1])-4

  L = zeros(length(tIndexes), length(xIndexes))
  for i in eachindex(xIndexes)
    xi = xIndexes[i]
    for j in eachindex(tIndexes)
      ti = tIndexes[j]

      L[j,i] = EI*FDM_Centre4(deformations[ti,:],xi) + m*FDM_Forward2(deformations[:,xi],ti) + C*FDM_Forward1(deformations[:,xi],ti) + k*deformations[ti,xi]



    end
  end

  return L
end


L = approximateFit(deformations, parameters, k1, C1)

xIndexes = 3:length(deformations[1,:])-3
tIndexes = 1:length(deformations[:,1])-4


plot(tVals[tIndexes], xVals[xIndexes],L, st=:surface)

plot(L, st=:surface)
=#
