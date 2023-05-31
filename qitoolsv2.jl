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



testtuple =(1,2,3,4,5)
length(testtuple)



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