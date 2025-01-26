
using Plots
default(legend=false)
plotlyjs()

#||||----Surface plot----||||#
myPlt = plot(xVals, tVals,  deformations, st=:surface)
display(myPlt)

#||||----Animation----||||#
@gif for ti in 1:tNum
  p = plot(xVals, deformations[ti,:], ylimits = (-0.01,0.001))
end
