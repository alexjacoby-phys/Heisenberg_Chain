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

function Q_reshape(L::Int64,A::Int64, B::Int64, Q::Int64, QA::Int64,QB::Int64,op::Matrix{ComplexF64})
    cvec = basis_converter(L,Q)
    svdmat = SparseArrays.spzeros(ComplexF64, )




function svd_op_reshape(op::Matrix{ComplexF64},sites::Tuple{Int64},L::Int64; Hdim::Int64 = 2, ϵ = 10e-15)
    opstate_svd = zeros(Float64,Hdim^(2*length(sites)), Hdim^(2L-2*length(sites)))
    for (ind, val) in pairs(op)
        a = embed(digits(ind[1]-1, base = Hdim),L)
        b = embed(digits(ind[2]-1, base = Hdim),L)
        for n in 1:L for m in 1:L
            
        end
    end
    return
end


in()