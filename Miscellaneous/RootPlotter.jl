using Plots
using Polynomials

AVals = LinRange(0, 100, 40)
BVals = LinRange(0, 100, 40)
CVals = LinRange(-100, 0, 40)


p = plot()
f_roots = []
for A in AVals, B in BVals, C in CVals
  f = Polynomial([C, 0, B, 0, A, 0, 1])
  global f_roots = vcat(f_roots, roots(f))
end


p = scatter!(real(f_roots), imag(f_roots), xlabel="Real", ylabel="Imaginary", legend=false)
