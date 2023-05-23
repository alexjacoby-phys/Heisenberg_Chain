include("H_Maker.jl")
include("qi_tools.jl")
import Plots


L = 8


@time H = make_HH(L)
#@time H += 2. *disorder_HH(L)
#@time KrylovKit.eigsolve(H)

PS1 = fill(1, L)
#PS1[1] = PS1[2] = 2
@time psi0 = ps_states(PS1, L)









δ = 0.2
τ = 0.:δ:3.



dat = zeros(Float64, 0)
for t in τ
    psi, info = KrylovKit.exponentiate(H, im * t, psi0)
    rho = Partial_Trace_Opstate(psi, L÷2, L-(L÷2), 0, L)
    append!(dat, SvN_Small(Matrix(rho + 10e-10 * LinearAlgebra.I)))
    println("time=", t)
end

Plots.plot(τ, dat)
#Plots.savefig("growth_L=10v3.pdf")


derivative = zeros(Float64,length(τ),length(τ)) + LinearAlgebra.I;
for i in 1:(length(τ)-1)
    derivative[i+1,i] = -1.
end
derivative * exp.(dat) / δ
slope = sum((derivative * exp.(dat) / δ)[2:7])/6


using LaTeXStrings
Plots.plot(slope*τ, exp.(dat), xlabel= L"Time, Units of $\left(\lim_{t\to 0 }\partial_{t}e^{\mathcal{S}_{\rm op}}\right)^{-1}$", ylabel=L"{ e^{\mathcal{S}_{\rm op}} }", #=seriestype=:scatter , =# color=:black, ytickfontsize=10, xtickfontsize=10, xguidefontsize=15, yguidefontsize=15, label = L"L=10"#=legend=:none , markersize=8=#)




# rho == rho'
# LinearAlgebra.tr(rho)
#LinearAlgebra.eigvals(Matrix(rho))