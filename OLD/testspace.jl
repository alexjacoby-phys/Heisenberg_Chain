
L = 16
Q = 8

A = 8
B = 8
QA = 4 
QAtilde = 4

QB = Q - QA
QBtilde = Q - QAtilde


@time Full_basis = basis_bitstrings(L, Q)

BasisA = basis_bitstrings(A, QA)
BasisB = basis_bitstrings(B, QB)
BasisAtilde = basis_bitstrings(A, QAtilde)
BasisBtilde = basis_bitstrings(B, QBtilde)

nba = [BasisA[i,:] for i in 1:size(BasisA,1)];
nbatilde = [BasisAtilde[i, :] for i in 1:size(BasisAtilde, 1)];
nbb = [BasisB[i, :] for i in 1:size(BasisB, 1)];
nbbtilde = [BasisBtilde[i, :] for i in 1:size(BasisBtilde, 1)];

nbadouble = vec(collect(Iterators.product(nba,nbatilde)));
nbbdouble = vec(collect(Iterators.product(nbb, nbbtilde)));

J =  collect(Iterators.product(nbadouble,nbbdouble));
J

nbadouble[1][2]
snorgles = vcat(nbadouble[1][1],nbbdouble[1][1],nbadouble[1][2],nbbdouble[1][2])



function get_index_very_specific(X::Tuple{Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}},Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}}};L::Int64)
    return (getbitindex(vcat(X[1][1], X[2][1]), L), getbitindex(vcat(X[1][2], X[2][2]), L))
end

get_index_very_specific(J[1,1],L=16)



dimA_double = size(BasisA, 1) * size(BasisAtilde, 1)
dimB_double = size(BasisB, 1) * size(BasisBtilde, 1)
basis_init = SparseArrays.spzeros(Int64,2*L)
double_basis_mat = fill(basis_init,size(BasisA,1)*size(BasisAtilde,1),size(BasisB,1)*size(BasisBtilde,1))
test = [ vcat(nbadouble[i][1],nbbdouble[j][1],nbadouble[i][2],nbbdouble[j][2]) for j in 1:size(double_basis_mat,2), i in 1:size(double_basis_mat,1)];
test2 = [ getbitindex(vcat(nbadouble[i][1], nbbdouble[j][1], nbadouble[i][2], nbbdouble[j][2]),2*L) for j in 1:dimB_double, i in 1:dimA_double];

double_basis_mat[:,:][1][1][1:A];

doubled_basis_A = vec([vcat(BasisA[i, :]..., BasisAtilde[j, :]...) for j in 1:size(BasisAtilde, 1), i in 1:size(BasisA, 1)])
doubled_basis_B = vec([vcat(BasisB[i, :]..., BasisBtilde[j, :]...) for j in 1:size(BasisBtilde, 1), i in 1:size(BasisB, 1)])




map_mat = collect(Iterators.product(doubled_basis_A, doubled_basis_B))




H = SparseArrays.spzeros(ComplexF64,2^L,2^L)
cvec = basis_converter(L,Q)
H[cvec, cvec] = H_Q(L, Q)


fullspace_index_mat = get_Fullspace_index.(A, B, map_mat, L)
new_indmat = [Qspace_dict(fullspace_index_mat[i, j], basis_dict) for i in 1:size(fullspace_index_mat, 1), j in 1:size(fullspace_index_mat, 2)]
