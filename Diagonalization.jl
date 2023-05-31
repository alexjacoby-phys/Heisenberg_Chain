include("H_Maker.jl")
include("qi_tools.jl")
import Plots

A = [ 0 1 ; 0 1]


for i in enumerate(A)
    print(i)
end
L = 14

H = make_H(L)


LinearAlgebra.eigen(Matrix(H))


# Matrix(H)



#  H = Matrix(make_H(L))


# psi  = KrylovKit.exponentiate(H, im * δ, psi)
# KrylovKit.exponentiate(H, im * δ, psi0)



@time H = make_HH(L)
#@time H = make_HH(L, [0.5, 0.5, 1.0])
#@time H += .1 *disorder_HH(L)
#@time KrylovKit.eigsolve(H)

PS1 = fill(1, L)
PS1[1] = PS1[2] = 2
@time psi0 = ps_states(PS1, L)

maxT = 1.0
δ = 0.01
τ = 0.0:δ:maxT


# dat = zeros(Float64, 0)
# for t in τ
#     psi, info = KrylovKit.exponentiate(H, im * t, psi0)
#     rho = Partial_Trace_Opstate(psi, L ÷ 2, L - (L ÷ 2), 0, L)
#     append!(dat, SvN_Small(Matrix(rho + 10e-10 * LinearAlgebra.I)))
#     println("time=", t)
# end


dat2 = zeros(Float64, 0)
psi = psi0
for t in τ
    psi, info = KrylovKit.exponentiate(H, im * δ, psi)
    append!(dat2, SvN(psi,L÷2,L-(L÷2),4))
    println("time=", t)
end


#Plots.plot(τ, dat)
Plots.plot(τ, dat2)
#Plots.savefig("growth_L=10v3.pdf")

dat = dat2
# derivative = zeros(Float64, length(τ), length(τ)) + LinearAlgebra.I;
# for i in 1:(length(τ)-1)
#     derivative[i+1, i] = -1.0
# end
# derivative * exp.(dat) / δ
# slope = sum((derivative*exp.(dat)/δ)[2:7]) / 6


using LaTeXStrings
# Plots.plot( slope *τ, exp.(dat), xlabel=L"Time, Units of $\left(\lim_{t\to 0 }\partial_{t}e^{\mathcal{S}_{\rm op}}\right)^{-1}$", ylabel=L"{ e^{\mathcal{S}_{\rm op}} }", color=:black, ytickfontsize=10, xtickfontsize=10, xguidefontsize=15, yguidefontsize=15, label=L"L=9") #=legend=:none , markersize=8=#

Plots.plot(τ, exp.(dat), xlabel=L"Time, No \  Units", ylabel=L"{ e^{\mathcal{S}_{\rm op}} }", color=:black, ytickfontsize=10, xtickfontsize=10, xguidefontsize=15, yguidefontsize=15, label=L"L=6") #=legend=:none , markersize=8=#




"Saves last figure with inputs L, δ, and τ"
function save(L::Int64, δ::Float64, τ::Float64)
    fn = string("Length=", L, "_", "delta=", δ, "_", "tau=", τ)
    Plots.savefig(string(fn, ".pdf"))
end
save(L, δ, maxT)


# rho == rho'
# LinearAlgebra.tr(rho)
#LinearAlgebra.eigvals(Matrix(rho))