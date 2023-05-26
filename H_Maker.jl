import LinearAlgebra, SparseArrays, KrylovKit, Plots, LaTeXStrings, DelimitedFiles
# import Pkg
# Pkg.add("LinearAlgebra")
# Pkg.add("SparseArrays")
# Pkg.add("KrylovKit")
# Pkg.add("LaTeXStrings")
# Pkg.add("Plots")
# Pkg.add("DelimitedFiles")

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
        B = 2 .* rand(Float64, 3) .- 1.0
        return B[1] * spX + B[2] * spY + B[3] * spZ
    end

    "Creates tensor of all the operators needed for a bilinear spin exchange. One can pass an anisotropy parameter-- a vector-- to get an anisotropic coupling. The components of the vector are the coefficient on the XX, YY, and ZZ couplings."
    function bilin(i::Int64, j::Int64, L::Int64)
        optens = [coupling(i, j, spX, spX, L), coupling(i, j, spY, spY, L), coupling(i, j, spZ, spZ, L)]
        return optens
    end
    function bilin(i::Int64, j::Int64, L::Int64, XYZ::Vector{Float64})
        if length(XYZ) != 3
            exit()
        end
        optens = [XYZ[1] * coupling(i, j, spX, spX, L), XYZ[2] * coupling(i, j, spY, spY, L), XYZ[3] * coupling(i, j, spZ, spZ, L)]
        return optens
    end


    function kron_array(oparray::Vector{Array{ComplexF64}})
        return kron(oparray...)
    end


    function kron_array(oparray::Vector{SparseArrays.SparseMatrixCSC{ComplexF64,Int64}})
        return kron(oparray...)
    end

    "Takes the output of bilin and converts it into an operator which can act directly on hilbert space.  One can pass an anisotropy parameter-- a vector-- to get an anisotropic coupling. The components of the vector are the coefficient on the XX, YY, and ZZ couplings."
    function full_coupling(i::Int64, j::Int64, L::Int64)
        return sum(kron_array.(bilin(i, j, L)))
    end
    function full_coupling(i::Int64, j::Int64, L::Int64,XYZ::Vector{Float64})
        return sum(kron_array.(bilin(i, j, L,XYZ)))
    end



    function lattice_mod(x::Int64, L::Int64)
        if x ÷ L == 0
            return x
        elseif x ÷ L != 0
            return mod(x - 1, L) + 1
        end
    end

    function make_H(L::Int64)
        H = SparseArrays.spzeros(2^L, 2^L)
        for i in 1:L
            H += full_coupling(lattice_mod(i, L), lattice_mod(i + 1, L), L)
        end
        return H
    end
    function make_H(L::Int64,XYZ::Vector{Float64})
        H = SparseArrays.spzeros(2^L, 2^L)
        for i in 1:L
            H += full_coupling(lattice_mod(i, L), lattice_mod(i + 1, L), L,XYZ)
        end
        return H
    end


    "Only parameter is the length. Will produce a doubled copy of the hamiltonian with the original on even sites and the copy on odd sites. If you pass an additional three component anisotropy vector, Zhu Li will do the thing."
    function make_HH(L::Int64)
        H = SparseArrays.spzeros(2^(2 * L), 2^(2 * L))
        for i in 1:L
            H += (full_coupling(lattice_mod(2i, 2 * L), lattice_mod(2i + 2, 2 * L), 2L) - full_coupling(lattice_mod(2i + 1, 2 * L), lattice_mod(2i + 3, 2 * L), 2L))
        end
        return H
    end
    function make_HH(L::Int64, XYZ::Vector{Float64})
         H = SparseArrays.spzeros(2^(2 * L), 2^(2 * L))
        for i in 1:L
            H += (full_coupling(lattice_mod(2i, 2 * L), lattice_mod(2i + 2, 2 * L), 2L,XYZ) - full_coupling(lattice_mod(2i + 1, 2 * L), lattice_mod(2i + 3, 2 * L), 2L, XYZ))
        end
        return H
    end

    function disorder_H(L::Int64)
        H = SparseArrays.spzeros(2^L, 2^L)
        for i in 1:L
            H += kron_array(OS(i, disorder_mat(), L))
        end
        return H
    end

    function disorder_HH(L::Int64)
        H = SparseArrays.spzeros(2^(2L), 2^(2L))
        for i in 1:L
            dm = disorder_mat()
            H += kron_array(OS(lattice_mod(2i, 2L), dm, 2L)) + kron_array(OS(lattice_mod(2i + 1, 2L), dm, 2L))
        end
        return H
    end
end