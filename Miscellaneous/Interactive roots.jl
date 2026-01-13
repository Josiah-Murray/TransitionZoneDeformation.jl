

using GLMakie
using Polynomials

fig = Figure()
ax1 = fig[1,1] = Axis(fig, xlabel="Real", ylabel="Imaginary")

lsgrid = SliderGrid(
    fig[2,1],
    (label = "D", range = 0:1:1000, startvalue = 50)
    , (label = "Q", range = 0:1:1000, startvalue = 50)
    , (label = "γ", range = -10:0.1:10, startvalue = 0)
    , (label = "s", range = 0:1:1000, startvalue = 50)
    , (label = "c", range = 0:1:1000, startvalue = 50))

sliderobservables = [s.value for s in lsgrid.sliders]

imaginary_root_indicator = Observable(:blue)

function root_indicator(x)
  if x<0
    return :red
  else
    return :blue
  end
end


A = lsgrid.sliders[1].value

function compute_roots(D, Q, γ, s, c)
  f = Polynomial([-s^2/c^2 , 0, (γ*s^2 +1), 0, Q, 0, D])
  f_roots = Polynomials.roots(f)
  return f_roots
end


points = lift(sliderobservables...) do D, Q, γ, s, c
    compute_roots(D, Q, γ, s, c)
end

discriminant = lift(sliderobservables...) do D, Q, γ, s, c
    Δ = (γ*s^2 + 1)^2*Q^2 - 4*(γ*s^2 + 1)^3*D -4*(-s^2/c^2)*Q^3 - 27*(-s^2/c^2)^2*D^2 + 18*(-s^2/c^2)*(γ*s^2 + 1)*Q*D
    return Δ
end



scatterplot = GLMakie.scatter!(ax1, @lift(real($points)), @lift(imag($points)), markersize=8, color=@lift(root_indicator($discriminant)))

fig
