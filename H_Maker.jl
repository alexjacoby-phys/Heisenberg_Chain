import LinearAlgebra, SparseArrays, KrylovKit

begin
    Id = Array{ComplexF64}([1 0; 0 1])
    X = Array{ComplexF64}([0 1; 1 0])
    Y = Array{ComplexF64}([0 -im; im 0])
    Z = Array{ComplexF64}([1 0; 0 -1])
    spId = SparseArrays.sparse(Id)
    spX = SparseArrays.sparse(X)
    spY = SparseArrays.sparse(Y)
    spZ = SparseArrays.sparse(Z)
    σ = [spId, spX, spY, spZ]

    function coupling(i::Int64, j::Int64, A::Array{ComplexF64}, B::Array{ComplexF64}, L::Int64)
        oparray = fill(spId, L)
        oparray[i] = A
        oparray[j] = B
        return oparray
    end


    function coupling(i::Int64, j::Int64, A::SparseArrays.SparseMatrixCSC{ComplexF64,Int64}, B::SparseArrays.SparseMatrixCSC{ComplexF64,Int64}, L::Int64)
        oparray = fill(spId, L)
        oparray[i] = A
        oparray[j] = B
        return oparray
    end


    function OS(i::Int64, A::Array{ComplexF64}, L::Int64)
        oparray = fill(spId, L)
        oparray[i] = A
        return oparray
    end


    function OS(i::Int64, A::SparseArrays.SparseMatrixCSC{ComplexF64,Int64}, L::Int64)
        oparray = fill(spId, L)
        oparray[i] = A
        return oparray
    end

    function disorder_mat()
        B = 2 .* rand(Float64,3) .- 1.
        return B[1]*spX + B[2]*spY + B[3] *spZ
    end


    function bilin(i::Int64, j::Int64, L::Int64)
        optens = [coupling(i, j, spX, spX, L), coupling(i, j, spY, spY, L), coupling(i, j, spZ, spZ, L)]
        return optens
    end


    function kron_array(oparray::Vector{Array{ComplexF64}})
        return kron(oparray...)
    end


    function kron_array(oparray::Vector{SparseArrays.SparseMatrixCSC{ComplexF64,Int64}})
        return kron(oparray...)
    end


    function full_coupling(i::Int64, j::Int64, L::Int64)
        return sum(kron_array.(bilin(i, j, L)))
    end


    function lattice_mod(x::Int64, L::Int64)
        if x ÷ L == 0
            return x
        elseif x ÷ L != 0
            return mod(x - 1, L) + 1
        end
    end

### SOMETHING WRONG HERE__ YOUR SPECTRUM IS NOT NOT INVERSION SYMMETRIC!!!!!!!
    function make_H(L::Int64)
        H = SparseArrays.spzeros(2^L, 2^L)
        for i in 1:L
            H += full_coupling(lattice_mod(i, L), lattice_mod(i + 1, L), L)
        end
        return H
    end


    #you might want to check all HH functions once again just to make sure they work correctly; your indexing was previously quite incorrect. 
    #Make_H is definitely correct since it replicated a gs energy. This replicates 2xgs energy, which I think it correct.
    function make_HH(L::Int64)
        H = SparseArrays.spzeros(2^(2 * L), 2^(2 * L))
        for i in 1:L
            H += (full_coupling(lattice_mod(2i, 2 * L), lattice_mod(2i + 2, 2 * L), 2L) - full_coupling(lattice_mod(2i + 1, 2 * L), lattice_mod(2i + 3, 2 * L), 2L))
        end
        return H
    end
    
    function disorder_H(L::Int64)
        H = SparseArrays.spzeros(2^L, 2^L)
        for i in 1:L
            H += kron_array(OS(i,disorder_mat(),L))
        end
        return H
    end

    function disorder_HH(L::Int64)
        H = SparseArrays.spzeros(2^(2L), 2^(2L))
        for i in 1:L
            dm = disorder_mat()
            H +=  kron_array(OS(lattice_mod(2i,2L),dm,2L)) + kron_array(OS(lattice_mod(2i+1,2L),dm,2L))
        end
        return H
    end
end