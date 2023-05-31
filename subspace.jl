include("H_Maker.jl")
import SparseArrays


function bitstrings(N::Int64)
    possibilities = collect(Iterators.product(Iterators.repeated((true, false), N)...))
    return reverse.(reshape(possibilities, (2^N)))
end


function basis_bitstrings(L::Int64,Q::Int64; sparse::Bool = true)
    if sparse == true
        states = SparseArrays.spzeros(Int64,binomial(L,Q),L)
    else
        states = zeros(binomial(L,Q),L)
    end
    n = 1
    for k in bitstrings(L)
        if sum(k) == Q
            states[n,:] = collect(Int64.(k))
            n += 1
        end
    end
    return states
end


function getbitindex(bitstring::Vector{Int64},L::Int64)
    return sum(bitstring .* (2 .^(0:(L-1))))
end
function getbitindex(bitstring::SparseArrays.SparseVector{Int64, Int64}, L::Int64)
    return sum(bitstring .* (2 .^ (0:(L-1))))
end

basis_bitstrings(10,3)

function basis_converter(L::Int64, Q::Int64)
    basis = basis_bitstrings(L,Q)
    D = binomial(L,Q)
    cvec = zeros(Int64,D)
    for i in 1:D
        cvec[i] = getbitindex(basis[i,:],L)
    end
    return cvec
end




function H_Q(Q::Int64,L::Int64)
    D = binomial(L,Q)
    H = SparseArrays.spzeros(D,D)
    H0 = make_H(L)
    cvec = basis_converter(L,Q)
    println("Starting Constructor")
    R = D^2
    for i in 1:D, j in 1:D
        H[i,j] = H0[cvec[i],cvec[j]]
    end
    return H 
end
H_Q(9,18)