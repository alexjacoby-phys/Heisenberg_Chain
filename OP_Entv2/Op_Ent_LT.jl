
include("Functions.jl")


L = 10
Q = 5
A = 4 
B = 6

T = 50.
δ = 0.02
τ = 0.0:δ:T






eigdat = LinearAlgebra.eigen(Matrix(H_Q(L,Q)))

Evec = eigdat.values
S = eigdat.vectors

dat = []
cvec = basis_converter(L, Q)
opinit = LinearAlgebra.normalize(kron(coupling(1, 2, spZ, spZ, L)...)[cvec, cvec])
psi = Op_time_evolution(Evec, S, 10., opinit)
for (QA,QB) in partitions, (QAtilde,QBtilde) in partitions
    CTC = QtoSVD_Basis_Tuple(L, Q, A, B, QA, QAtilde)
    svdmat = [psi[entry...] for entry in CTC]
    append!(dat,LinearAlgebra.svdvals(svdmat)...)
end
dat
sum(dat .^2)
# QA = 3
# QAtilde = 2
for (tind, t ) in pairs(τ)
    println(tind, "     ", t)
end


entropies = []
partitions = Charge_Partitions(L,Q,A,B)
for (tind, t) in pairs(τ)
    dat = []
    opinit = LinearAlgebra.normalize(kron(coupling(1, 2, spZ, spZ, L)...)[cvec, cvec])
    psi = Op_time_evolution(Evec, S, t, opinit)
    for (QA,QB) in partitions, (QAtilde,QBtilde) in partitions
        CTC = QtoSVD_Basis_Tuple(L, Q, A, B, QA, QAtilde)
        svdmat = [psi[entry...] for entry in CTC]
        schmidt_vals = LinearAlgebra.svdvals(svdmat)
        es = schmidt_vals .^2 + ones(Float64,length(schmidt_vals))*10e-15
        append!(dat,es)
    end
    svn_vec = - dat .* log.(dat)
    entropy = sum(svn_vec)
    append!(entropies, entropy)
    println(t)
    println(string("normalized to ", sum(dat)))
end

import Plots
Plots.plot(τ,entropies)
Plots.savefig("U(1)Proof_of_concept.pdf")


# dat = zeros(Float64,length(τ),binomial(L,Q))
##something quite weird here. Dat comes out uniform in time so time evolution not doing anything. Seems like an indexing error.


dat = fill([],length(τ))
Charge_Partitions(L,Q,A,B)







HQ = Matrix(H_Q(L,Q)) 

for (itind,t) in pairs(τ)
    println(string("time = ",t))
    psi = Op_time_evolution(Evec,S,t,opinit)
    for (QA,QAtilde) in Charge_Partitions(L,Q,A,B)
        CTC  =  QtoSVD_Basis_Tuple(L, Q, A, B, QA , QAtilde)
        opqstate = [psi[entry...] for entry in CTC]
        append!(dat[itind],LinearAlgebra.svdvals(opqstate))
    end
end
dat

data = [ (entry .^2 )/(sum( entry .^2)) for entry in dat]

-LinearAlgebra.dot(data[1], log.(data[1]))

Y = [-LinearAlgebra.dot(entry, log.(entry)) for entry in data]











# function get_index_very_specific(X::Tuple{Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}},Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}}}; length::Int64)
#     bitindex_Qspace = vcat(X[1][1], X[2][1])
#     bitindex_Qtildespace = vcat(X[1][2], X[2][2])
#     index_Qspace = getbitindex(bitindex_Qspace, length)
#     index_Qtildespace = getbitindex(bitindex_Qtildespace, length)
#     return (index_Qspace, index_Qtildespace)
# end


# function get_index_very_specific(X::Tuple{Tuple{Vector{Int64},Vector{Int64}},Tuple{Vector{Int64},Vector{Int64}}}; length::Int64)
#     bitindex_Qspace = vcat(X[1][1], X[2][1])
#     bitindex_Qtildespace = vcat(X[1][2], X[2][2])
#     index_Qspace = getbitindex(bitindex_Qspace, length)
#     index_Qtildespace = getbitindex(bitindex_Qtildespace, length)
#     return (index_Qspace, index_Qtildespace)
# end





# function indx_shuffle(X::Tuple{Tuple{Vector{Int64},Vector{Int64}},Tuple{Vector{Int64},Vector{Int64}}}; length::Int64)
#     return (vcat(X[1][1], X[2][1]), vcat(X[1][2], X[2][2]))
# end

# QB = Q - QA
# QBtilde = Q - QAtilde

# Full_basis = basis_bitstrings(L, Q)
# basis_dict = Dict([(getbitindex(Full_basis[i, :], L), i) for i in 1:size(Full_basis, 1)])

# BasisA = Matrix(basis_bitstrings(A, QA))
# BasisB = Matrix(basis_bitstrings(B, QB))
# BasisAtilde = Matrix(basis_bitstrings(A, QAtilde))
# BasisBtilde = Matrix(basis_bitstrings(B, QBtilde))

# nba = [BasisA[i, :] for i in 1:size(BasisA, 1)]
# nbatilde = [BasisAtilde[i, :] for i in 1:size(BasisAtilde, 1)]
# nbb = [BasisB[i, :] for i in 1:size(BasisB, 1)]
# nbbtilde = [BasisBtilde[i, :] for i in 1:size(BasisBtilde, 1)]

# nbadouble = vec(collect(Iterators.product(nba, nbatilde)))
# nbbdouble = vec(collect(Iterators.product(nbb, nbbtilde)))

# J = collect(Iterators.product(nbadouble, nbbdouble))
# display(J)
# indx_shuffle.(J,length =L)
# #something wrong here should be symmetric
# P = get_index_very_specific.(J, length=L)

# [(basis_dict[entry[1]], basis_dict[entry[2]]) for entry in P]





function Q_op(L::Int64)
    Q_op_init = SparseArrays.spzeros(ComplexF64, 2^L, 2^L)
    Q_op_init[LinearAlgebra.diagind(Q_op_init)] = [L - sum(digits(i, base=2)) for i in 0:((2^L)-1)]
    return Q_op_init
end



L=10
Q=4

cvec = basis_converter(L,Q)

Q_MAT = Q_op(L)


[ Q_MAT[i,i] for i in cvec]











begin
    # # This is useful if you ever forget which way around the change of basis should go.
    # LinearAlgebra.norm(S* LinearAlgebra.diagm(eigdat.values)* S' - H)
    # LinearAlgebra.norm(S' * LinearAlgebra.diagm(eigdat.values) * S - H)
end




