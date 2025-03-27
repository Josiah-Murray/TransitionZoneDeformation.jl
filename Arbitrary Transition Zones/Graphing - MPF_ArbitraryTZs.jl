
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

  #If y limits aren't specified, come up with some
  if (yLims==false)
    println("Automatic y limits:")
    if (steadyStateDeformations!=false)
      yLims = [Inf, -Inf]
      for s in SteadyStateDeformations
        s_max  = maximum(s)
        s_min = minimum(s)
        if s_min < yLims[1]
          yLims[1] = maximum(s)
        end
        if s_max > yLims[2]
          yLims[2] = s_max
        end
      end
      def_max = maximum(deformations)
      def_min = minimum(deformations)

      if def_min < yLims[1]
        yLims[1] = def_min
      end
      if def_max > yLims[2]
        yLims[2] = def_max
      end
    else
      yLims = [minimum(deformations), maximum(deformations)]
    end

    #Add some fraction of the height to top and bottom
    height = yLims[2]-yLims[1]
    fraction = 0.1
    yLims = [yLims[1]-height*fraction, yLims[2] + height*fraction]

    println(yLims)
   end

  @gif for ti in 1:tNum
    p=plot()
    #plot lines for transition zones
    for tz in tz_list
      p = vline!([tz.location],color=RGB(0.1,0.1,0.2), width=3, foreground_color=RGB(0.7,0.7,0.7), background_color=RGB(0.05,0.05,0.1), framestyle=:box)
    end

    #Plot line for moving point force
    p = vline!([v*(ti-1)*(tMax-tMin)/tNum + xp], width=3, color=RGB(0.3,0.3,0.4))

    #Plot graph
    p = plot!(xVals, data[ti,:], lc=RGB(0.5,0.4,0.9), lw=3, ylimits = yLims, xlimits = (xLeft,xRight))

    #Optionally plot steady state solutions
    if(steadyStateDeformations != false)
      p = PlotSteadyStates(p, ti, SteadyStateDeformations)
    end


    xlabel!("Beam coordinate (m)")
    ylabel!("Deformation (m)")
  end

end


#Takes an input between 0 and 1 and outputs a linearly interpolated (in rgb space)
#colour between the two.
function CustomGradient(input, colour1,colour2)

  return colour1 + input*(colour2-colour1)

end

#For a set of tIndexes (should probably be no more than about 5?), plot each of the
#relevant deformation profiles on the same graph.
#If no desired ylims or steady states, set these values to false.
function GraphTimeEvolution(data, xVals, tVals, tIndexes, parameters, steadyStateDeformations, yLims, colour1, colour2)

  ~, ~, xp, tz_list, ~, v = parameters
  tNum = length(tVals)
  xLeft = xVals[1]
  xRight = xVals[end]

  #If y limits aren't specified, come up with some
  if (yLims==false)
    println("Automatic y limits:")
    if (steadyStateDeformations!=false)
      yLims = [Inf, -Inf]
      for s in SteadyStateDeformations
        s_max  = maximum(s)
        s_min = minimum(s)
        if s_min < yLims[1]
          yLims[1] = maximum(s)
        end
        if s_max > yLims[2]
          yLims[2] = s_max
        end
      end
      def_max = maximum(deformations)
      def_min = minimum(deformations)

      if def_min < yLims[1]
        yLims[1] = def_min
      end
      if def_max > yLims[2]
        yLims[2] = def_max
      end
    else
      yLims = [minimum(deformations), maximum(deformations)]
    end

    #Add some fraction of the height to top and bottom
    height = yLims[2]-yLims[1]
    fraction = 0.1
    yLims = [yLims[1]-height*fraction, yLims[2] + height*fraction]

    println(yLims)
  end

  p = plot()

  #plot lines for transition zones
  for tz in tz_list
    p = vline!([tz.location],color=RGB(0.8,0.8,0.9), width=3, framestyle=:box)
  end


  for i in eachindex(tIndexes)


    ti = reverse(tIndexes)[i]

    gradLineColour = CustomGradient(i/length(tIndexes), colour1,colour2)

    #Plot line for moving point force
    p = vline!([v*(ti-1)*(tMax-tMin)/tNum + xp], width=1, color=gradLineColour+0*(RGB(1,1,1)-gradLineColour))




    p = plot!(xVals, data[ti,:], lc = gradLineColour, lw=4, ylimits = yLims, xlimits = (xLeft,xRight))
    if(steadyStateDeformations != false)
      for j in eachindex(SteadyStateDeformations)
        p = plot!(xVals, SteadyStateDeformations[j][ti,:], linestyle = :dash, color = RGB(0.5,0.2,0.15), lw=3.5)
      end
    end
  end
  p = plot!(framestyle=:box)
   display(p)
end


#For a set of tIndexes (should probably be no more than about 5?), plot each of the
#relevant deformation profiles on the same graph.
#If no desired ylims or steady states, set these values to false.
function MultiPlotTimeEvolution(data, xVals, tVals, tIndexes, parameters, steadyStateDeformations, yLims, colour1, colour2)

  ~, ~, xp, tz_list, ~, v = parameters
  tNum = length(tVals)
  xLeft = xVals[1]
  xRight = xVals[end]

  #If y limits aren't specified, come up with some
  if (yLims==false)
    println("Automatic y limits:")
    if (steadyStateDeformations!=false)
      yLims = [Inf, -Inf]
      for s in SteadyStateDeformations
        s_max  = maximum(s)
        s_min = minimum(s)
        if s_min < yLims[1]
          yLims[1] = maximum(s)
        end
        if s_max > yLims[2]
          yLims[2] = s_max
        end
      end
      def_max = maximum(deformations)
      def_min = minimum(deformations)

      if def_min < yLims[1]
        yLims[1] = def_min
      end
      if def_max > yLims[2]
        yLims[2] = def_max
      end
    else
      yLims = [minimum(deformations), maximum(deformations)]
    end

    #Add some fraction of the height to top and bottom
    height = yLims[2]-yLims[1]
    fraction = 0.1
    yLims = [yLims[1]-height*fraction, yLims[2] + height*fraction]

    println(yLims)
  end






  for i in eachindex(tIndexes)
    p = plot()
    #plot lines for transition zones
    for tz in tz_list
      p = vline!([tz.location],color=RGB(0.8,0.8,0.9), width=3, framestyle=:box)
    end

    ti = tIndexes[i]

    gradLineColour = CustomGradient(i/length(tIndexes), colour2,colour1)

    #Plot line for moving point force
    p = vline!([v*(ti-1)*(tMax-tMin)/tNum + xp], width=1, color=gradLineColour+0*(RGB(1,1,1)-gradLineColour))




    p = plot!(xVals, data[ti,:], lc = gradLineColour, lw=4, ylimits = yLims, xlimits = (xLeft,xRight))
    if(steadyStateDeformations != false)
      for j in eachindex(SteadyStateDeformations)
        p = plot!(xVals, SteadyStateDeformations[j][ti,:], linestyle = :dash, color = RGB(0.5,0.2,0.15), lw=3.5)
      end
    end
    p = plot!(framestyle=:box)
    title!(string(tVals[ti])) #TODO: Update how this works
    display(p)
  end

end


#For a set of tIndexes (should probably be no more than about 5?), plot each of the
#relevant deformation profiles on the same graph.
#If no desired ylims or steady states, set these values to false.
function MultiPlotTimeEvolution_ReturnPlotList(data, xVals, tVals, tIndexes, parameters, steadyStateDeformations, yLims, colour1, colour2)

  ~, ~, xp, tz_list, ~, v = parameters
  tNum = length(tVals)
  xLeft = xVals[1]
  xRight = xVals[end]

  #If y limits aren't specified, come up with some
  if (yLims==false)
    println("Automatic y limits:")
    if (steadyStateDeformations!=false)
      yLims = [Inf, -Inf]
      for s in SteadyStateDeformations
        s_max  = maximum(s)
        s_min = minimum(s)
        if s_min < yLims[1]
          yLims[1] = maximum(s)
        end
        if s_max > yLims[2]
          yLims[2] = s_max
        end
      end
      def_max = maximum(deformations)
      def_min = minimum(deformations)

      if def_min < yLims[1]
        yLims[1] = def_min
      end
      if def_max > yLims[2]
        yLims[2] = def_max
      end
    else
      yLims = [minimum(deformations), maximum(deformations)]
    end

    #Add some fraction of the height to top and bottom
    height = yLims[2]-yLims[1]
    fraction = 0.1
    yLims = [yLims[1]-height*fraction, yLims[2] + height*fraction]

    println(yLims)
  end





  p_list = []

  for i in eachindex(tIndexes)
    p = plot()
    #plot lines for transition zones
    for tz in tz_list
      p = vline!([tz.location],color=RGB(0.8,0.8,0.9), width=3, framestyle=:box)
    end

    ti = tIndexes[i]

    gradLineColour = CustomGradient(i/length(tIndexes), colour2,colour1)

    #Plot line for moving point force
    p = vline!([v*(tVals[i]) + xp], width=3, color=colour2)




    p = plot!(xVals, data[ti,:], lc = colour1, lw=4, ylimits = yLims, xlimits = (xLeft,xRight))
    if(steadyStateDeformations != false)
      for j in eachindex(SteadyStateDeformations)
        p = plot!(xVals, SteadyStateDeformations[j][ti,:], linestyle = :dash, color = RGB(0.5,0.2,0.15), lw=3.5)
      end
    end
    p = plot!(framestyle=:box)
    title!(string(round(tVals[ti], digits=3))*" s" , titlefontsize = 12,line=-10)
    xlabel!("Beam coordinate (m)", xguidefontsize=10)
    ylabel!("Deformation (m)", yguidefontsize=10)

    p_list = [p_list; p]
  end

  return p_list

end

#Takes in a set of deformations defined for 10 time values and plots each time in a 5 by 2 grid.
#Note that, when called in a file, the line `gr()` should be added first, else the formatting will be wrong.
#For whatever reason, adding it to the function itself, doesn't help.
#It also needs to be added before each call, not just once in the document.
#As I have come to say often in my PhD, 'Everything is broken and I don't know why'.
function Graph10Times(deformations, xVals, tVals, parameters, SteadyStateDeformations, ylims, colour1, colour2, storageFolderPath, graphName)
  p_list = MultiPlotTimeEvolution_ReturnPlotList(deformations, xVals, tVals, 1:10, parameters, SteadyStateDeformations, ylims, colour1, colour2)
  p1 = plot(p_list[1:5]..., layout = grid(5,1), size = (1200,600))
  p2 = plot(p_list[6:10]..., layout = grid(5,1), size = (1200,600))
  p = plot(p1,p2, layout = grid(1,2), size = (1200,1200),leftmargin=10Plots.mm)
  display(p)
  savefig(p, storageFolderPath*graphName*".png")
  savefig(p, storageFolderPath*graphName*".pdf")
end
