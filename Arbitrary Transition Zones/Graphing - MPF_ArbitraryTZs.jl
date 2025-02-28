
using Plots
default(legend=false)
plotlyjs()

#||||----Surface plot----||||#
myPlt = plot(xVals, tVals,  deformations, st=:surface)
display(myPlt)

#||||----Animation----||||#

function PlotSteadyStates(p, ti, SteadyStateDeformations)

  for i in eachindex(SteadyStateDeformations)
    p = plot!(xVals, SteadyStateDeformations[i][ti,:], linestyle = :dash, color = RGB(1,1,1))
  end

  return p
end

@gif for ti in 1:tNum
  p=plot()
  for tz in xtz_list
    p = vline!([tz.location],color=RGB(0.1,0.1,0.2), width=3, foreground_color=RGB(0.7,0.7,0.7), background_color=RGB(0.05,0.05,0.1), framestyle=:box)
  end

  p = vline!([v*(ti-1)*(tMax-tMin)/tNum + xp], width=3, color=RGB(0.3,0.3,0.4))
  p = plot!(xVals, deformations[ti,:], lc=RGB(0.5,0.4,0.9), lw=3, ylimits = (-0.001,0.001), xlimits = (xLeft,xRight))
  p = PlotSteadyStates(p, ti, SteadyStateDeformations)


  xlabel!("Beam coordinate (m)")
  ylabel!("Deformation (m)")
end

println("Finished Graphing")


##
