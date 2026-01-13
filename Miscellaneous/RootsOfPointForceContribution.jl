using GLMakie



using GLMakie
#using Polynomials


DataPath = "Paper--Moving point force/Data"
functionFolder = "../Constant Point Force/"
figureFolder = "Paper--Moving point force/Figures"

include(joinpath(functionFolder, "Functions - Inversion schemes.jl"))
include(joinpath(functionFolder, "Functions - MPF_ArbitraryTZs.jl"))
include(joinpath(functionFolder, "Graphing - MPF_ArbitraryTZs.jl"))



fig = Figure()
ax1 = fig[1,1] = Axis(fig, xlabel="Real", ylabel="Imaginary")
ax2 = fig[1,2] = Axis3(fig, xlabel= "Real", ylabel="Imaginary", zlabel="Magnitude")

lsgrid = SliderGrid(
    fig[2,1],
    (label = "EI", range = 0:1:1000, startvalue = 50)
    , (label = "m", range = 0:1:1000, startvalue = 50)
    , (label = "C", range = 0:01:1000, startvalue = 1)
    , (label = "k", range = 0:1:1000, startvalue = 50)
    , (label = "v", range = 0:1:10, startvalue = 1)
)

sliderobservables = [s.value for s in lsgrid.sliders]


A = lsgrid.sliders[1].value

function compute_roots(EI, m, C, k, v)
  f = [k, C, m, 0, EI/v^4]
  f_roots = roots(f)
  return f_roots
end


points = lift(sliderobservables...) do  EIv4, m, C, k, v
    compute_roots(EIv4, m, C, k, v)
end

sImVls = LinRange(-1,1,50)
sReVls = LinRange(-1,1,50)


laplaceSpaceSolution = lift(sliderobservables...) do  EI, m, C, k,v
    parameters = [EI, m, 0, [TransitionZone(1, k,k,C,C)], 1, v]

    [LaplaceSpaceFunctionMovingPointForce1TZ(-1,(s_r + 1im*s_i), parameters) for s_r in sReVls, s_i in sImVls]
end

function ComplexColour(z)
  θ = angle(z)
  scaled = (θ + π)
  shift = π/2
  index = 1.5
  dimness = 0.5
  return RGBf(
    ((sin(scaled        + shift)+1)/2)^(index*log(abs(z)+1))*dimness,
    ((sin(scaled + 2π/3 + shift)+1)/2)^(index*log(abs(z)+1))*dimness,
    ((sin(scaled + 4π/3 + shift)+1)/2)^(index*log(abs(z)+1))*dimness
    )
end
colorMapping = lift(laplaceSpaceSolution) do lss
    ComplexColour.(lss)
end



scatterplot = GLMakie.scatter!(ax1, @lift(real($points)), @lift(imag($points)), markersize=8)
p = GLMakie.surface!(ax2, sReVls, sImVls, @lift(abs.($laplaceSpaceSolution)) )
fig
