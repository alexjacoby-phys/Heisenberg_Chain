include("H_Maker.jl")
import SparseArrays


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


#Work in Q = 3 subspace on A and Q = 1 subspace on Atilde
L = 7
Q = 4

A = 3
B = 4


QA = 2
QB = Q - QA

QAtilde = 1
QBtilde = Q - QAtilde





BasisA = basis_bitstrings(A,QA)
BasisB = basis_bitstrings(B,QB)
BasisAtilde = basis_bitstrings(A,QAtilde)
BasisBtilde = basis_bitstrings(B,QBtilde)

doubled_basis_A = vec([vcat(BasisA[i,:]...,BasisAtilde[j,:]...) for j in 1:size(BasisAtilde,1), i in 1:size(BasisA,1)])
doubled_basis_B = vec([vcat(BasisB[i, :]..., BasisBtilde[j, :]...) for j in 1:size(BasisBtilde, 1), i in 1:size(BasisB, 1)])

map_mat = collect(Iterators.product(doubled_basis_A,doubled_basis_B))

test = map_
ind_vec = vcat(test[1][1:A]...,test[2][1:B]...,test[1][(A+1):2A]...,test[2][(B+1):2B]...)
getbitindex(ind_vec,2*L)

function get_H_Q_index(A::Int64,B::Int64,svd_index::Tuple{Vector{Int64},Vector{Int64}})
    ind_vec = vcat(svd_index[1][1:A]...,svd_index[2][1:B]...,svd_index[1][(A+1):2A]...,svd_index[2][(B+1):2B]...)
    return getbitindex(ind_vec,2*L)
end

get_H_Q_index.(A,B,map_mat)