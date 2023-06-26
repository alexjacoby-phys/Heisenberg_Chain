include("Dependencies.jl")


#Hamiltonian Maker
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
    function full_coupling(i::Int64, j::Int64, L::Int64, XYZ::Vector{Float64})
        return sum(kron_array.(bilin(i, j, L, XYZ)))
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
    function make_H(L::Int64, XYZ::Vector{Float64})
        H = SparseArrays.spzeros(2^L, 2^L)
        for i in 1:L
            H += full_coupling(lattice_mod(i, L), lattice_mod(i + 1, L), L, XYZ)
        end
        return H
    end

    function make_H(L::Int64, XYZ::Vector{Float64})
        H = SparseArrays.spzeros(2^L, 2^L)
        for i in 1:L
            H += full_coupling(lattice_mod(i, L), lattice_mod(i + 1, L), L, XYZ)
        end
        return H
    end

    function TFI(L::Int64,g::Float64; disorder::Bool = false)
        H = SparseArrays.spzeros(Float64,2^L,2^L)
        for i in 1:L
            if disorder
                ds = (2*(rand()-0.5))*g
            else
                ds = g
            end
            coupling = fill(spId, L)
            coupling[lattice_mod(i,L)] = spZ
            coupling[lattice_mod(i+1,L)] = spZ
            H+= kron_array(coupling)
            os = fill(spId, L)
            os[lattice_mod(i,L)] = spX
            H += (ds * kron_array(os))
        end
        return H
    end


end




#Core bitstrings and Qspace functionalities
begin
    function bitstrings(N::Int64)
        possibilities = collect(Iterators.product(Iterators.repeated((false, true), N)...))
        return collect.(reverse.(reshape(possibilities, (2^N))))
    end



    function basis_bitstrings(L::Int64, Q::Int64; sparse::Bool=true)
        if sparse == true
            states = SparseArrays.spzeros(Int64, binomial(L, Q), L)
        else
            states = zeros(binomial(L, Q), L)
        end
        n = 1
        for k in bitstrings(L)
            if sum(k) == L-Q
                states[n, :] = collect(Int64.(k))
                n += 1
            end
        end
        return states
    end

    function getbitindex(bitstring::Vector{Int64}, L::Int64)
        return (sum(bitstring .* (2 .^ ((L-1):-1:0))) + 1)
    end

    function getbitindex(bitstring::SparseArrays.SparseVector{Int64,Int64}, L::Int64)
        return (sum(bitstring .* (2 .^ ((L-1):-1:0)))+1)
    end

    function getbitindex(bitstring::SparseArrays.SparseVector{Int64,Bool}, L::Int64)
        return (sum(bitstring .* (2 .^ ((L-1):-1:0))) + 1)
    end

    function getbitindex(bitstring::Vector{Bool}, L::Int64)
        return (sum(bitstring .* (2 .^ ((L-1):-1:0))) + 1)
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
end


begin
    "This function is a very specific helper function. It is not recommended that you use it."
    function get_index_very_specific(X::Tuple{Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}},Tuple{SparseArrays.SparseVector{Int64,Int64},SparseArrays.SparseVector{Int64,Int64}}}; length::Int64)
        return (getbitindex(vcat(X[1][1], X[2][1]), length), getbitindex(vcat(X[1][2], X[2][2]), length))
    end

    "To transform an operator use [op[entry...] for entry in output_of_this_function]. SOMETHING GOING WRONG HERE!!!!!"
    function QtoSVD_Basis_Tuple(L::Int64, Q::Int64, A::Int64, B::Int64, QA::Int64, QAtilde::Int64)
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

    function Charge_Partitions(L::Int64,Q::Int64,A::Int64,B::Int64)
        if Q > L
            exit()
        end
        Partitions = Vector{Vector{Int64}}([])
        for QA in 0:Q
            QB = Q - QA
            if ( QA ≤ A ) && ( QB ≤ B )
                append!(Partitions, [[QA, QB]])
            end
        end
        return Partitions
    end
end


#time evolution
begin
    function Time_evolution(Evec::Vector{Float64}, S::Matrix{ComplexF64},t::Float64)
        return S*SparseArrays.spdiagm(exp.(- im * t * Evec))*S'
    end

    function Op_time_evolution(Evec::Vector{Float64}, S::Matrix{ComplexF64},t::Float64, Op::Matrix{ComplexF64})
        U = Time_evolution(Evec::Vector{Float64}, S::Matrix{ComplexF64}, t::Float64)
        return U' * Op * U
    end
    function Op_time_evolution(Evec::Vector{Float64}, S::Matrix{ComplexF64}, t::Float64, Op::SparseArrays.SparseMatrixCSC{ComplexF64, Int64})
        U = Time_evolution(Evec::Vector{Float64}, S::Matrix{ComplexF64}, t::Float64)
        return U' * Op * U
    end
end


