include("H_Maker.jl")


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



    IdVECS0 = LinearAlgebra.eigen(Id)
    XVECS0 = LinearAlgebra.eigen(X)
    YVECS0 = LinearAlgebra.eigen(Y)
    ZVECS0 = LinearAlgebra.eigen(Z)
    ##Yes the identity is a bit stupid. Just nice to treat the pauli strings consistently.

    IdVECS = (LinearAlgebra.normalize(IdVECS0.values), SparseArrays.sparse(IdVECS0.vectors), SparseArrays.sparse(conj.(IdVECS0.vectors)))
    XVECS = (LinearAlgebra.normalize(XVECS0.values), SparseArrays.sparse(XVECS0.vectors), SparseArrays.sparse(conj.(XVECS0.vectors)))
    YVECS = (LinearAlgebra.normalize(YVECS0.values), SparseArrays.sparse(YVECS0.vectors), SparseArrays.sparse(conj.(YVECS0.vectors)))
    ZVECS = (LinearAlgebra.normalize(ZVECS0.values), SparseArrays.sparse(ZVECS0.vectors), SparseArrays.sparse(conj.(ZVECS0.vectors)))
    # The first field holds indexed eigenvalues, the second field holds the eigenvectors (indexed as [:,i]), and the third field is the same as the second but with a complex conjugate (NOT HERMITIAN CONJUGATE)



    IdS = kron(IdVECS[2][:, 1], IdVECS[3][:, 1]) * IdVECS[1][1] + kron(IdVECS[2][:, 2], IdVECS[3][:, 2]) * IdVECS[1][2]

    XS = kron(XVECS[2][:, 1], XVECS[3][:, 1]) * XVECS[1][1] + kron(XVECS[2][:, 2], XVECS[3][:, 2]) * XVECS[1][2]

    YS = kron(YVECS[2][:, 1], YVECS[3][:, 1]) * YVECS[1][1] + kron(YVECS[2][:, 2], YVECS[3][:, 2]) * YVECS[1][2]

    ZS = kron(ZVECS[2][:, 1], ZVECS[3][:, 1]) * ZVECS[1][1] + kron(ZVECS[2][:, 2], ZVECS[3][:, 2]) * ZVECS[1][2]




    pstates = [IdS, XS, YS, ZS]

    states = SparseArrays.sparse.(Vector{ComplexF64}.(([1, 0], [0, 1])))

    function ps_states(PS::Vector{Int64}, L::Int64)
        if length(PS) != L
            exit()
        end
        kp = collect(Iterators.repeated(SparseArrays.spzeros(ComplexF64, 4), L))
        for n in 1:L
            kp[n] = pstates[PS[n]]
        end
        return kron(kp...)
    end

    function bitstrings(N::Int64)
        possibilities = collect(Iterators.product(Iterators.repeated((true, false), N)...))
        return reshape(possibilities, (2^N))
    end
end



begin
    function ProjLR_N(NL::Int64, N::Int64, NR::Int64)
        #returns NORMALIZED projectors with NL identities on the left and NR identities on the right (i.e. projectors on the middle N site Hilbert space). If NL or NR are 0, it will be the strings that touch an edge. Don't set N to 0 (I have no idea why you would want to do this), though setting NL or NR to 0 is supported and necessary.
        if ((NL == 0) & (NR == 0))
            #print("case1")
            return vec(collect(Iterators.product(Iterators.repeated(proj, N)...)))
        elseif ((NL == 0) & (NR != 0))
            #print("case2")
            return vec(collect(Iterators.product(Iterators.repeated(proj, N)..., Iterators.repeated([spI_N], NR)...)))
        elseif ((NL != 0) & (NR == 0))
            #print("case3")
            return vec(collect(Iterators.product(Iterators.repeated([spI_N], NL)..., Iterators.repeated(proj, N)...)))
        elseif ((NL != 0) & (NR != 0))
            #print("case4")
            return vec(collect(Iterators.product(Iterators.repeated([spI_N], NL)..., Iterators.repeated(proj, N)..., Iterators.repeated([spI_N], NR)...)))
        else
            return "Error: whatever you did broke this function"
        end
    end
    function ProjLR(NL::Int64, N::Int64, NR::Int64)
        #returns UNNORMALIZED projectors with NL identities on the left and NR identities on the right (i.e. projectors on the middle N site Hilbert space). If NL or NR are 0, it will be the strings that touch an edge. Don't set N to 0 (I have no idea why you would want to do this), though setting NL or NR to 0 is supported and necessary.
        if ((NL == 0) & (NR == 0))
            #print("case1")
            return vec(collect(Iterators.product(Iterators.repeated(proj, N)...)))
        elseif ((NL == 0) & (NR != 0))
            #print("case2")
            return vec(collect(Iterators.product(Iterators.repeated(proj, N)..., Iterators.repeated([spId], NR)...)))
        elseif ((NL != 0) & (NR == 0))
            #print("case3")
            return vec(collect(Iterators.product(Iterators.repeated([spId], NL)..., Iterators.repeated(proj, N)...)))
        elseif ((NL != 0) & (NR != 0))
            #print("case4")
            return vec(collect(Iterators.product(Iterators.repeated([spId], NL)..., Iterators.repeated(proj, N)..., Iterators.repeated([spId], NR)...)))
        else
            return "Error: whatever you did broke this function"
        end
    end
    function PJT(NL::Int64, N::Int64, NR::Int64, normalized::Bool=false)
        #removes idenity and gives the previous function a slightly more convenient name
        if !(normalized)
            Opvec = ProjLR(NL, N, NR)
        elseif normalized
            Opvec = ProjLR_N(NL, N, NR)
        end
        return Opvec
    end
end







begin
    ### Quantum Information Tools
    function Down_Projectors_LR(NL::Int64, N::Int64, NR::Int64)
        #returns UNNORMALIZED projectors with NL identities on the left and NR identities on the right (I projectors on the middle N site Hilbert space). These projectors also have the property of collapsing the relevant Hilbert space. IE this is to be used for partial traces. If NL or NR are 0, it will be the strings that touch an edge. Don't set N to 0 (I have no idea why you would want to do this), though setting NL or NR to 0 is supported and necessary.
        if ((NL == 0) & (NR == 0))
            #print("case1")
            return vec(collect(Iterators.product(Iterators.repeated(states, N)...)))
        elseif ((NL == 0) & (NR != 0))
            #print("case2")
            return vec(collect(Iterators.product(Iterators.repeated(states, N)..., Iterators.repeated([spId], NR)...)))
        elseif ((NL != 0) & (NR == 0))
            #print("case3")
            return vec(collect(Iterators.product(Iterators.repeated([spId], NL)..., Iterators.repeated(states, N)...)))
        elseif ((NL != 0) & (NR != 0))
            #print("case4")
            return vec(collect(Iterators.product(Iterators.repeated([spId], NL)..., Iterators.repeated(states, N)..., Iterators.repeated([spId], NR)...)))
        else
            return "Error: whatever you did broke this function"
        end
    end
    function Partial_Trace(rho::Matrix{ComplexF64}, NL::Int64, N::Int64, NR::Int64)
        # performs partial trace over N center sites. I would prefer not to deal with Disjoint partial traces, and they can be achieved by composing this function.
        down_converter = Down_Projectors_LR(NL, N, NR)
        newrho = zeros(ComplexF64, 2^(NL + NR), 2^(NL + NR))
        for i in 1:length(down_converter)
            newrho = newrho + adjoint(kron(down_converter[i]...)) * rho * kron(down_converter[i]...)
        end
        return newrho
    end
    "Note that this function is defined in the real space coordinate system. Not the doubled hilbert space. That is, NL+N+NR = L where L is the physical length."
    function Partial_Trace_Opstate(psi::SparseArrays.SparseVector{ComplexF64,Int64}, NL::Int64, N::Int64, NR::Int64, L::Int64)
        if L != (NR + N + NL)
            exit()
        end
        rho = SparseArrays.spzeros(4^(NL + NR), 4^(NL + NR))
        down_converter = Down_Projectors_LR(2 * NL, 2 * N, 2 * NR)
        for (k, Γ) in enumerate(down_converter)
            psi_k = kron(Γ...)' * psi
            rho += kron(psi_k, psi_k')
        end
        return rho
    end
    function SvN_Small(rho::Matrix{ComplexF64})
        Dat = LinearAlgebra.eigvals(rho)
        SvN = -1.0 * LinearAlgebra.dot(Dat, log.(Dat))
        return SvN
    end
    "Takes a state vector and makes bipartition according to A, B. The default on site dimension is set to two but can be changed freely. Furthermore, an ϵ is included to ensure convergence with a default value of 10^-14. You must have loaded linear algebra for this to work"
    function SvN(psi::Vector{ComplexF64}, A::Int64, B::Int64, dim::Int64=2, ϵ::Float64=10e-15)
        svn_mat = reshape(psi, dim^A, dim^B)
        spectrum = LinearAlgebra.svdvals(svn_mat) .^ 2
        svn_vec = spectrum .* log.(spectrum + LinearAlgebra.fill(ϵ, dim^(min(A, B))))
        return -sum(svn_vec)
    end
    function SvN(psi::SparseArrays.SparseVector{ComplexF64,Int64}, A::Int64, B::Int64, dim::Int64=2, ϵ::Float64=10e-15)
        svn_mat = reshape(psi, dim^A, dim^B)
        spectrum = LinearAlgebra.svdvals(Matrix(svn_mat)) .^ 2
        svn_vec = spectrum .* log.(spectrum + LinearAlgebra.fill(ϵ, dim^(min(A, B))))
        return -sum(svn_vec)
    end
end
#this is my trash bin of potentially useful code that is not currently in use
begin
    # this was not the right way to code this particular thing, but the code might be useful later
    # "Transforms a vector with entries 1:4 for [I,X,Y,Z] to that Pauli string's operator state. This function maps H⊗H* non-locally as H = 0:L an H*  = L+1:2L. Partial traces should be conducted with this in mind."
    # function karray_state_init(PS::Vector{Int64}, L::Int64)
    #     if length(PS) != L
    #         exit
    #     end
    #     psimat = SparseArrays.spzeros(ComplexF64,2^(2 * L), 2^L)
    #     for (k, σ) in enumerate(bitstrings(L))
    #         psi_k = collect(Iterators.repeated(SparseArrays.spzeros(ComplexF64,2), 2L))
    #         for n in 1:L
    #             psi_k[n] = ps_bitstring[PS[n]][2][:, σ[n]+1] ## YOU ARE MISSING EIGENVALUES!!!
    #             psi_k[n+L] = ps_bitstring[PS[n]][3][:,σ[n]+1]
    #         end
    #         psimat[:, k] = kron(psi_k...)
    #     end
    #     return psimat
    # end




    # function ps_state(PS::Vector{Int64},L::Int64)
    #     return sum(karray_state_init(PS,L), dims = (2))[:,1]
    # end

    # ps_bitstring = [IdVECS, XVECS, YVECS, ZVECS]

    # "Returns all bitstrings of length L as Int8 tuples (int8 used to minimize the memory allocation from this process)."
    # function bitstrings(N::Int64)
    #     possibilities = collect(Iterators.product(Iterators.repeated((true, false), N)...))
    #     return reshape(possibilities, (2^N))
    # end

    # function Partial_Trace(rho::Matrix{ComplexF64}, NL::Int64, N::Int64, NR::Int64)
    #     # performs partial trace over N center sites. I would prefer not to deal with Disjoint partial traces, and they can be achieved by composing this function.
    #     down_converter = Down_Projectors_LR(NL, N, NR)
    #     newrho = zeros(ComplexF64, 2^(NL + NR), 2^(NL + NR))
    #     for i in 1:length(down_converter)
    #         newrho = newrho + adjoint(kron(down_converter[i]...)) * rho * kron(down_converter[i]...)
    #     end
    #     return newrho
    # end
    # function Partial_Trace(rho::SparseArrays.SparseMatrixCSC{ComplexF64,Int64}, NL::Int64, N::Int64, NR::Int64)
    #     # performs partial trace over N center sites. I would prefer not to deal with Disjoint partial traces, and they can be achieved by composing this function.
    #     down_converter = Down_Projectors_LR(NL, N, NR)
    #     newrho = zeros(ComplexF64, 2^(NL + NR), 2^(NL + NR))
    #     for i in 1:length(down_converter)
    #         newrho = newrho + adjoint(kron(down_converter[i]...)) * rho * kron(down_converter[i]...)
    #     end
    #     return newrho
    # end
    # function Partial_Trace(rho::Matrix{ComplexF64}, NL::Int64, N::Int64, NR::Int64)
    #     # performs partial trace over N center sites. I would prefer not to deal with Disjoint partial traces, and they can be achieved by composing this function.
    #     down_converter = Down_Projectors_LR(NL, N, NR)
    #     newrho = zeros(ComplexF64, 2^(NL + NR), 2^(NL + NR))
    #     for i in 1:length(down_converter)
    #         newrho = newrho + adjoint(kron(down_converter[i]...)) * rho * kron(down_converter[i]...)
    #     end
    #     return newrho
    # end
end





# function choi_state(L::Int64,)
#     D =2^(2L)
#     psi = SparseArrays.spzeros(D)
#     for R in bistrings(L)
#         psi += kron
    