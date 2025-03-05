
using Plots
default(legend=false)
plotlyjs()

#||||----Surface plot----||||#
#myPlt = plot(xVals, tVals,  deformations, st=:surface)
#display(myPlt)

#||||----Animation----||||#

function PlotSteadyStates(p, ti, SteadyStateDeformations)

  for i in eachindex(SteadyStateDeformations)
    p = plot!(xVals, SteadyStateDeformations[i][ti,:], linestyle = :dash, color = RGB(1,1,1))
  end

  return p
end

#Data should be an array with data[t, x]
function GifWithFeatures(data, xVals, tVals, parameters, steadyStateDeformations = false, yLims = false)

  ~, ~, xp, tz_list, ~, v = parameters
  tNum = length(tVals)
  xLeft = xVals[1]
  xRight = xVals[end]

  @gif for ti in 1:tNum
    p=plot()
    #plot lines for transition zones
    for tz in tz_list
      p = vline!([tz.location],color=RGB(0.1,0.1,0.2), width=3, foreground_color=RGB(0.7,0.7,0.7), background_color=RGB(0.05,0.05,0.1), framestyle=:box)
    end

    #Plot line for moving point force
    p = vline!([v*(ti-1)*(tMax-tMin)/tNum + xp], width=3, color=RGB(0.3,0.3,0.4))

    #Plot graph, optionally specifying ylimits
    if(ylims!=false)
      p = plot!(xVals, data[ti,:], lc=RGB(0.5,0.4,0.9), lw=3, ylimits = (-0.0001,0.0001), xlimits = (xLeft,xRight))
    else
      p = plot!(xVals, data[ti,:], lc=RGB(0.5,0.4,0.9), lw=3,  xlimits = (xLeft,xRight))
    end

    #Optionally plot steady state solutions
    if(steadyStateDeformations != false)
      p = PlotSteadyStates(p, ti, SteadyStateDeformations)
    end


    xlabel!("Beam coordinate (m)")
    ylabel!("Deformation (m)")
  end

end

println("Finished Graphing")


##
#=
#||||----Difference Graph----||||#

difference1 = SteadyStateDeformations[1] - deformations
difference2 = SteadyStateDeformations[2] - deformations

@gif for ti in 1:tNum
  p=plot()
  for tz in xtz_list
    p = vline!([tz.location],color=RGB(0.1,0.1,0.2), width=3, foreground_color=RGB(0.7,0.7,0.7), background_color=RGB(0.05,0.05,0.1), framestyle=:box)
  end

  p = vline!([v*(ti-1)*(tMax-tMin)/tNum + xp], width=3, color=RGB(0.3,0.3,0.4))
  p = plot!(xVals, difference1[ti,:], lc=RGB(0.5,0.4,0.9), lw=3, ylimits = (-0.0001,0.0001), xlimits = (xLeft,xRight))
  p = plot!(xVals, difference2[ti,:], lc=RGB(0.5,0.9,0.4), lw=3, ylimits = (-0.0001,0.0001), xlimits = (xLeft,xRight))



  xlabel!("Beam coordinate (m)")
  ylabel!("Deformation (m)")
end
=#
