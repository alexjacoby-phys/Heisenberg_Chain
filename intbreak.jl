include("H_Maker.jl")
include("qi_tools.jl")


H = Matrix(make_H(14))
LinearAlgebra.eigen(H)

"Enter in order L, T, δT, disorder strengths (format as a tuple), and then provide a file name."
function run_me(L::Int64, T::Float64, δ::Float64, disorder_strengths::Vector{Float64}, ps::Vector{Int64}, filename::String; ϵ::Float64=10e-15)
    τ = 0.0:δ:T
    dat_array = zeros(Float64, length(disorder_strengths), length(τ))
    H0 = make_HH(L)
    H1 = disorder_HH(L)
    psi0 = ps_states(ps, L)
    for (n, disorder) in enumerate(disorder_strengths)
        H = H0 + disorder * H1
        psi = psi0
        data = zeros(Float64, length(τ))
        for (k, t) in enumerate(τ)
            psi, convergence = KrylovKit.exponentiate(H, im * δ, psi)
            rho = Partial_Trace_Opstate(psi, L ÷ 2, L - (L ÷ 2), 0, L)
            data[k] = SvN_Small(Matrix(rho + ϵ * LinearAlgebra.I))
            println("time= ", t)
        end
        dat_array[n, :] = data
    end
    home = pwd()
    mkdir(filename)
    cd(string(home, "/", filename))
    info = string("info_", filename, ".txt")
    dat = string("dat_", filename, ".txt")
    touch(info)
    touch(dat)
    io = open(dat, "a")
    DelimitedFiles.writedlm(io, dat_array)
    close(io)
    io = open(info, "w")
    println(io, string("System size = ", L))
    println(io, string("Max Time = ", T))
    println(io, string("Time Step = ", δ))
    println(io, string("Disorder Strengths = ", disorder_strengths))
    close(io)
    cd(home)
    return dat_array
end



L = 10
T = 5.0
δ = 0.02
τ = 0.0:δ:T
disorder_strengths = [0.0, 0.05, 0.1, 0.25, 0.5]
PS1 = fill(1, L)
PS1[1] = PS1[2] = 2


data = run_me(L, T, δ, disorder_strengths, PS1, "Trial3")




# dat_array = zeros(length(disorder_strengths),length(τ))
# H = make_HH(L)
# HD = disorder_HH(L)
# touch("testfile.txt")
# io = open("testfile.txt", "a")
# println(io, "more test stuff")
# println(io, "more test stuff")
# println(io, "more test stuff")
# close(io)


