#I am experiencing some problems with Weeks inversion for my moving load problem.
#It has occurred to me that this is probably because Weeks' method has trouble dealing
#with the long down time before the point force reaches it.
#This code aims to confirm that this is indeed a problem with Weeks' method.

#It turns out I was wrong :P

using Plots
using SpecialFunctions
using InverseLaplace

t_shift = 5

function approximateWithShift(t_shift)
  tVals = range(0,10,50)

  F(s) = (sqrt(π)/2)exp(s^2/4)erfc(s/2)*exp(-s*t_shift)

  f_e(t) = exp(-(t-t_shift)^2)


  f_t = Weeks(F, 16 )
  f_approx = zeros(length(tVals),1)
  for i in eachindex(tVals)
    f_approx[i] = f_t(tVals[i])
  end



  p = plot(tVals, f_e, xlims=[0,10], ylims = [-2,2])
  p = plot!(tVals, f_approx)
  display(p)
end



for i in range(0,10,10)
  approximateWithShift(i)
end
