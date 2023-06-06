using LinearAlgebra, SparseArrays



A = (1,2,3,4)

A[4]
pairs(spX)


for (index,vals) in pairs(Y)
    println(index[1],"   ", index[2], "   ", vals)
end





function embed(dig::Vector{Int64}, L::Int64)
    append!(dig, zeros(Int64, L - length(dig)))
    return dig
end



basis_bitstrings(10,3)

basis_converter(10,3)

# function Q_reshape(L::Int64,A::Int64, B::Int64, Q::Int64, QA::Int64,QB::Int64,op::Matrix{ComplexF64})
#     cvec = basis_converter(L,Q)
#     svdmat = SparseArrays.spzeros(ComplexF64, )


L = 10 
psi = LinearAlgebra.normalize(rand(ComplexF64,2^L))
sites = (1, 2 , 4) 

sort([4,2,1])


function svd_state_reshape(state:Vector{ComplexF64},sites_prelim::Vector{Int64} ;Hdim = 2)
    sites = 
    basis_bitstrings(length(sites))

end



function Q_state_reshape()
    
end




