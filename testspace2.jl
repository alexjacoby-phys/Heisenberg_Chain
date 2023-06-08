include("subspace.jl")


L = 16
Q = 8

A = 7
B = 9
QA = 4
QAtilde = 3

QB = Q - QA
QBtilde = Q - QAtilde


@time Full_basis = basis_bitstrings(L, Q);

BasisA = basis_bitstrings(A, QA)
BasisB = basis_bitstrings(B, QB)
BasisAtilde = basis_bitstrings(A, QAtilde)
BasisBtilde = basis_bitstrings(B, QBtilde)

nba = [BasisA[i, :] for i in 1:size(BasisA, 1)];
nbatilde = [BasisAtilde[i, :] for i in 1:size(BasisAtilde, 1)];
nbb = [BasisB[i, :] for i in 1:size(BasisB, 1)];
nbbtilde = [BasisBtilde[i, :] for i in 1:size(BasisBtilde, 1)];

nbadouble = vec(collect(Iterators.product(nba, nbatilde)));
nbbdouble = vec(collect(Iterators.product(nbb, nbbtilde)));

J = collect(Iterators.product(nbadouble, nbbdouble));

function get_index_very_specific(X::Tuple{Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}},Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}}}; length::Int64)
    return (getbitindex(vcat(X[1][1], X[2][1]), length), getbitindex(vcat(X[1][2], X[2][2]), length))
end
@time P = get_index_very_specific.(J,length=L);

basis_dict = Dict([(getbitindex(Full_basis[i,:], L), i) for i in 1:size(Full_basis,1)])

#Check that these two matrices correspond to each other correctly
P
QtoAB = [ (basis_dict[entry[1]],basis_dict[entry[2]]) for entry in P]














function get_index_very_specific(X::Tuple{Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}},Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}}}; length::Int64)
    return (getbitindex(vcat(X[1][1], X[2][1]), length), getbitindex(vcat(X[1][2], X[2][2]), length))
end

function QtoSVD_Basis(; L::Int64, Q::Int64, A::Int64, B::Int64, QA::Int64, QAtilde::Int64)
    Full_basis = basis_bitstrings(L, Q)
    basis_dict = Dict([(getbitindex(Full_basis[i, :], L), i) for i in 1:size(Full_basis, 1)])

    BasisA = basis_bitstrings(A, QA)
    BasisB = basis_bitstrings(B, QB)
    BasisAtilde = basis_bitstrings(A, QAtilde)
    BasisBtilde = basis_bitstrings(B, QBtilde)

    nba = [BasisA[i, :] for i in 1:size(BasisA, 1)]
    nbatilde = [BasisAtilde[i, :] for i in 1:size(BasisAtilde, 1)]
    nbb = [BasisB[i, :] for i in 1:size(BasisB, 1)]
    nbbtilde = [BasisBtilde[i, :] for i in 1:size(BasisBtilde, 1)]

    nbadouble = vec(collect(Iterators.product(nba, nbatilde)))
    nbbdouble = vec(collect(Iterators.product(nbb, nbbtilde)))

    J = collect(Iterators.product(nbadouble, nbbdouble))

    P = get_index_very_specific.(J, length=L)

    return [(basis_dict[entry[1]], basis_dict[entry[2]]) for entry in P]
end
























NO = [H[entry...] for entry in QtoAB]
H = LinearAlgebra.normalize(H_Q(L,Q))

snogs = LinearAlgebra.normalize(spX)

SparseArrays.sparse(NO)

sum(LinearAlgebra.svdvals(NO))