
# using Pkg
# Pkg.activate("jj8144/home/.julia/LOE_Code")
# Pkg.instantiate()

include("Functions.jl")


home = pwd()

fob = "/Users/alexjacoby/Documents/Research_Code/Heisenberg_Chain/OP_ENT_DELLAV2/OP_DATA"
index_fob = "/Users/alexjacoby/Documents/Research_Code/Heisenberg_Chain/OP_ENT_DELLAV2/INDEX_DATA"
# mkdir(fob)
cd(fob)

L_min = 4
L_max = 4
L_RNG = L_min:1:L_max

T_min = 400.
T_max = 500.
δ = 5.
τ = T_min : δ : T_max

touch("Time.txt")
DelimitedFiles.writedlm("Time.txt",Vector(τ))

for L in L_RNG
    fn = string("L=", L, ".txt")
    touch(fn)

    A = L ÷ 2
    B = L - A
    Q = L ÷ 2
    L_Index_fob = string(index_fob, "/L=", L)
    partitions = Int64.(DelimitedFiles.readdlm(string(L_Index_fob, "/Partition_Double.txt")))
    partition_double = [partitions[i, :] for i in 1:size(partitions)[1]]

    cvec = Int64.(dropdims(DelimitedFiles.readdlm(string(L_Index_fob, "/QBasis.txt")), dims=(2)))

    opinit = LinearAlgebra.normalize(kron(coupling(1, 2, spZ, spZ, L)...)[cvec, cvec])

    eigdat = LinearAlgebra.eigen(Matrix(make_H(L)[cvec, cvec]))
    Evec = eigdat.values
    S = eigdat.vectors

    dat_mat = [Vector{Float64}([]) for i in τ]


    for (t_index, t) in pairs(τ)
        psi = Op_time_evolution(Evec, S, t, opinit)
        dat_vec = []
        for (p_no, (QA, QB, QAtilde, QBtilde)) in pairs(partition_double)
            CTC_L_FILE = string(L_Index_fob, "/L_Partition_number=", p_no, ".txt")
            CTC_R_FILE = string(L_Index_fob, "/R_Partition_number=", p_no, ".txt")
            CTCL = Int64.(DelimitedFiles.readdlm(CTC_L_FILE))
            CTCR = Int64.(DelimitedFiles.readdlm(CTC_R_FILE))
            CTC = [(CTCL[i, j], CTCR[i, j]) for i in 1:size(CTCL)[1], j in 1:size(CTCL)[2]]
            cd(string(fob, "/",))
            svdmat = [psi[entry...] for entry in CTC]
            append!(dat_vec, LinearAlgebra.svdvals(svdmat)...)
            CTC = nothing
            svdmat = nothing
            GC.gc()
        end
        append!(dat_mat[t_index], (dat_vec .^ 2)...)
        psi = nothing
        data_vec = nothing
    end
    DelimitedFiles.writedlm(fn, dat_mat)
    println("L = ", L, " Done")
end
cd(home)







