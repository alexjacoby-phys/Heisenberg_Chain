include("H_Maker.jl")



##Please just generally check this code is correct because I feel like it is doing weird shit
function bitstrings(N::Int64)
    possibilities = collect(Iterators.product(Iterators.repeated((true, false), N)...))
    return reverse.(reshape(possibilities, (2^N)))
end

function basis_bitstrings(L::Int64, Q::Int64; sparse::Bool=true)
    if sparse == true
        states = SparseArrays.spzeros(Int64, binomial(L, Q), L)
    else
        states = zeros(binomial(L, Q), L)
    end
    n = 1
    for k in bitstrings(L)
        if sum(k) == Q
            states[n, :] = collect(Int64.(k))
            n += 1
        end
    end
    return states
end


function getbitindex(bitstring::Vector{Int64}, L::Int64)
    return sum(bitstring .* (2 .^ (0:(L-1))))
end

function getbitindex(bitstring::SparseArrays.SparseVector{Int64,Int64}, L::Int64)
    return sum(bitstring .* (2 .^ (0:(L-1))))
end



function basis_converter(L::Int64, Q::Int64)
    basis = basis_bitstrings(L, Q)
    D = binomial(L, Q)
    cvec = zeros(Int64, D)
    for i in 1:D
        cvec[i] = getbitindex(basis[i, :], L)
    end
    return cvec
end


function H_Q(L::Int64, Q::Int64)
    cvec = basis_converter(L, Q)
    H0 = make_H(L)
    return H0[cvec, cvec]
end


function get_Fullspace_index(A::Int64, B::Int64, svd_index::Tuple{Vector{Int64},Vector{Int64}}, L::Int64)
    ind_vec = vcat(svd_index[1][1:A]..., svd_index[2][1:B]..., svd_index[1][(A+1):2A]..., svd_index[2][(B+1):2B]...)
    return getbitindex(ind_vec, 2 * L)
end



function Qspace_dict(full_entry::Int64, basis_map::Dict{Int64,Int64})
    return basis_map[full_entry]
end


# "This function is outdated and shouldn't be used. It is correct but it will crash julia for even reasonably sized chains."
# function operator_reshape_ind_mat(;L::Int64, A::Int64,B::Int64,Q::Int64, QA::Int64, QAtilde::Int64)
#     if (A + B) != L
#         return "You are a big dumb idiot"
#     end
#     QB = Q - QA
#     QBtilde = Q - QAtilde


#     Full_basis = basis_bitstrings(L, Q)
#     Full_double_basis = vec([vcat(Full_basis[i, :]..., Full_basis[j, :]...) for j in 1:size(Full_basis, 1), i in 1:size(Full_basis, 1)])
#     basis_dict = Dict([(getbitindex(Full_double_basis[i], 2 * L), i) for i in 1:length(Full_double_basis)])

#     BasisA = basis_bitstrings(A, QA)
#     BasisB = basis_bitstrings(B, QB)
#     BasisAtilde = basis_bitstrings(A, QAtilde)
#     BasisBtilde = basis_bitstrings(B, QBtilde)

#     doubled_basis_A = vec([vcat(BasisA[i, :]..., BasisAtilde[j, :]...) for j in 1:size(BasisAtilde, 1), i in 1:size(BasisA, 1)])
#     doubled_basis_B = vec([vcat(BasisB[i, :]..., BasisBtilde[j, :]...) for j in 1:size(BasisBtilde, 1), i in 1:size(BasisB, 1)])

#     map_mat = collect(Iterators.product(doubled_basis_A, doubled_basis_B))
#     fullspace_index_mat = get_Fullspace_index.(A, B, map_mat, L)
#     new_indmat = [Qspace_dict(fullspace_index_mat[i, j], basis_dict) for i in 1:size(fullspace_index_mat, 1), j in 1:size(fullspace_index_mat, 2)]
#     return new_indmat
# end


function get_index_very_specific(X::Tuple{Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}},Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}}}; length::Int64)
    return (getbitindex(vcat(X[1][1], X[2][1]), length), getbitindex(vcat(X[1][2], X[2][2]), length))
end

function tupledex(R::Tuple{Int64,Int64}, N::Int64)
    return R[1] + (R[2] - 1) * N
end







#"I think this is wrong????"
# function QtoSVD_Basis(; L::Int64, Q::Int64, A::Int64, B::Int64, QA::Int64, QAtilde::Int64)
#     QB = Q - QA
#     QBtilde = Q - QAtilde

#     Full_basis = basis_bitstrings(L, Q)
#     basis_dict = Dict([(getbitindex(Full_basis[i, :], L), i) for i in 1:size(Full_basis, 1)])

#     BasisA = basis_bitstrings(A, QA)
#     BasisB = basis_bitstrings(B, QB)
#     BasisAtilde = basis_bitstrings(A, QAtilde)
#     BasisBtilde = basis_bitstrings(B, QBtilde)

#     nba = [BasisA[i, :] for i in 1:size(BasisA, 1)]
#     nbatilde = [BasisAtilde[i, :] for i in 1:size(BasisAtilde, 1)]
#     nbb = [BasisB[i, :] for i in 1:size(BasisB, 1)]
#     nbbtilde = [BasisBtilde[i, :] for i in 1:size(BasisBtilde, 1)]

#     nbadouble = vec(collect(Iterators.product(nba, nbatilde)))
#     nbbdouble = vec(collect(Iterators.product(nbb, nbbtilde)))

#     J = collect(Iterators.product(nbadouble, nbbdouble))

#     P = get_index_very_specific.(J, length=L)

#     return tupledex.([(basis_dict[entry[1]], basis_dict[entry[2]]) for entry in P],length(Full_basis))
# end


"To transform an operator use [op[entry...] for entry in output_of_this_function]"
function QtoSVD_Basis_Tuple(; L::Int64, Q::Int64, A::Int64, B::Int64, QA::Int64, QAtilde::Int64)
    QB = Q - QA
    QBtilde = Q - QAtilde

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


# N= 973
# Avec = Vector{Int64}(1:N^2)
# A = reshape(Avec, N,N)
# B = [ (i,j) for i in 1:N, j in 1:N]

# function tupledex(R::Tuple{Int64,Int64},N::Int64)
#     return R[1] + (R[2]-1)*N 
# end

# tupledex.(B,N) == A

















begin
    #Work in Q = 3 subspace on A and Q = 1 subspace on Atilde
    # L = 7
    # Q = 4

    # A = 3
    # B = 4


    # QA = 2
    # QB = Q - QA

    # QAtilde = 1
    # QBtilde = Q - QAtilde



    # Full_basis =  basis_bitstrings(L,Q)
    # Full_double_basis = vec([vcat(Full_basis[i, :]..., Full_basis[j, :]...) for j in 1:size(Full_basis, 1), i in 1:size(Full_basis, 1)])

    # basis_dict = Dict([(getbitindex(Full_double_basis[i], 2*L), i) for i in 1:length(Full_double_basis)])
    # #basis_dict = Dict([ (getbitindex(Full_basis[i,:],L),i) for i in 1:size(Full_basis,1)])




    # fullspace_index_mat = get_Fullspace_index.(A,B,map_mat)







    # new_indmat = [Qspace_dict(fullspace_index_mat[i, j], basis_dict) for j in 1:size(fullspace_index_mat, 2), i in 1:size(fullspace_index_mat, 1)]




    # LinearAlgebra.svdvals(nm)




    # eigendat = LinearAlgebra.eigen(Matrix(H_Q(16, 7)))

    # S = eigendat.vectors'
    # t = 1.0
    # U = S' * LinearAlgebra.diagm(exp.(-im * t * eigendat.values)) * S





    # function H_Qv1(L::Int64, Q::Int64; sparse::Bool=false)
    #     D = binomial(L, Q)
    #     H = SparseArrays.spzeros(D, D)
    #     H0 = make_H(L)
    #     cvec = basis_converter(L, Q)
    #     println("Starting Constructor")
    #     R = D^2
    #     for i in 1:D, j in 1:D
    #         H[i, j] = H0[cvec[i], cvec[j]]
    #     end
    #     return H
    # end

    # function H_Qv2(L::Int64, Q::Int64)
    #     Dim = binomial(L, Q)
    #     cvec = basis_converter(L, Q)
    #     indmat = collect(Iterators.product(cvec, cvec))
    #     H0 = make_H(L)
    #     H = SparseArrays.spzeros(ComplexF64, Dim, Dim)
    #     for (Qind, fullind) in pairs(indmat)
    #         H[Qind] = H0[fullind...]
    #     end
    #     return H
    # end

    # function H_Qv3(L::Int64,Q::Int64)
    #     cvec = basis_converter(L,Q)
    #     H0 = make_H(L)
    #     return H0[cvec,cvec]
    # end
end
