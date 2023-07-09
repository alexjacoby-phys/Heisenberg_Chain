include("Functions.jl")



L=14

T_min = 400.
T_max = 420.
δ = 5.
τ = T_min : δ : T_max



A = L ÷ 2
B = L - A
Q = L ÷ 2
partitions = Charge_Partitions(L, Q, A, B)

cvec = basis_converter(L, Q)
opinit = LinearAlgebra.normalize(kron(coupling(1, 2, spZ, spZ, L)...)[cvec, cvec])

@time eigdat = LinearAlgebra.eigen(Matrix(H_Q(L, Q)))
Evec = eigdat.values;
S = eigdat.vectors;

dat_mat = [Vector{Float64}([]) for i in τ]

@time for (t_index, t) in pairs(τ)
    @time psi = Op_time_evolution(Evec, S, t, opinit);
    dat_vec = []
    for (QA, QB) in partitions, (QAtilde, QBtilde) in partitions
        CTC = QtoSVD_Basis_Tuple(L, Q, A, B, QA, QAtilde)
        svdmat = [psi[entry...] for entry in CTC]
        append!(dat_vec, LinearAlgebra.svdvals(svdmat)...);
        CTC = nothing
        svdmat = nothing
        GC.gc()d
    end
    append!(dat_mat[t_index], (dat_vec .^2)...)
    psi = nothing
    data_vec = nothing
    @time GC.gc()
end

eigdat = nothing