include("Functions.jl")


home = pwd()

fob = "OP_ENT_DAT_Vtest"
mkdir(fob)
cd(fob)

L_min = 10
L_max = 15
L_RNG = L_min:1:L_max

T_min = 0.
T_max = 2.
δ = 0.05
τ = T_min : δ : T_max

touch("Time.txt")
DelimitedFiles.writedlm("Time.txt",Vector(τ))

for L in L_RNG
    fn = string("L=", L, ".txt")
    touch(fn)

    A = L ÷ 2
    B = L - A
    Q = L ÷ 2
    partitions = Charge_Partitions(L, Q, A, B)

    cvec = basis_converter(L, Q)
    opinit = LinearAlgebra.normalize(kron(coupling(1, 2, spZ, spZ, L)...)[cvec, cvec])

    eigdat = LinearAlgebra.eigen(Matrix(H_Q(L, Q)))
    Evec = eigdat.values
    S = eigdat.vectors

    dat_mat = [Vector{Float64}([]) for i in τ]

    for (t_index, t) in pairs(τ)
        psi = Op_time_evolution(Evec, S, t, opinit)
        dat_vec = []
        for (QA, QB) in partitions, (QAtilde, QBtilde) in partitions
            CTC = QtoSVD_Basis_Tuple(L, Q, A, B, QA, QAtilde)
            svdmat = [psi[entry...] for entry in CTC]
            append!(dat_vec, LinearAlgebra.svdvals(svdmat)...)
        end
        append!(dat_mat[t_index], (dat_vec .^2)...)

    end
    DelimitedFiles.writedlm(fn, dat_mat)
    println("L = ",L, " Done")
end