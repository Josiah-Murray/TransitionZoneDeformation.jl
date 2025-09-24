

using GLMakie
using Polynomials

fig = Figure()
ax1 = fig[1,1] = Axis(fig, xlabel="Real", ylabel="Imaginary")

lsgrid = SliderGrid(
    fig[2,1],
    (label = "A", range = 0:1:100, startvalue = 50)
    , (label = "B", range = 0:1:100, startvalue = 50)
    , (label = "C", range = -100:1:0, startvalue = -50))

sliderobservables = [s.value for s in lsgrid.sliders]

imaginary_root_indicator = Observable(:blue)

function root_indicator(x)
  if x<50
    return :red
  else
    return :blue
  end
end


A = lsgrid.sliders[1].value

function compute_roots(A,B,C)
  f = Polynomial([C, 0, B, 0, A, 0, 1])
  f_roots = roots(f)
  return f_roots
end


points = lift(sliderobservables...) do A, B, C
    compute_roots(A, B, C)
end



scatterplot = scatter!(ax1, @lift(real($points)), @lift(imag($points)), markersize=8, color=@lift(root_indicator($A)))

fig
