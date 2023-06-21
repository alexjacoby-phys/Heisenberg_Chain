
include("Functions.jl")



basis_bitstrings(10,0)

L = 14
Q = 7 
A = 7 
B = 7

T = 5.0
δ = 0.02
τ = 0.0:δ:T













eigdat = LinearAlgebra.eigen(Matrix(H_Q(L,Q)))

Evec = eigdat.values
S = eigdat.vectors

cvec = basis_converter(L, Q)
opinit = kron(coupling(1, 2, spZ, spZ, L)...)[cvec, cvec]
Op_time_evolution(Evec,S,10.,opinit)

t = 3.41
@time S * SparseArrays.spdiagm(exp.(-im * t * Evec)) * S'
@time S * LinearAlgebra.diagm(exp.(-im * t * Evec)) * S'







CTC  =  QtoSVD_Basis_Tuple(L = 10, Q = 5, A = 4, B = 6, QA = 2, QAtilde = 3)

A = LinearAlgebra.normalize(H_Q(10,5))

sdfhds = SparseArrays.sparse([ A[entry...] for entry in CTC])


H = H_Q(14,7)

eigdat = LinearAlgebra.eigen(Matrix(H))
eigdat.values

τ = 2
U = eigdat.vectors * LinearAlgebra.diagm(exp.(-im * τ * eigdat.values)) * eigdat.vectors'
LinearAlgebra.norm(U*U' - LinearAlgebra.I)

LinearAlgebra.norm(U' * H * U - H)

# This is useful if you ever forget which way around the change of basis should go.
LinearAlgebra.norm(S* LinearAlgebra.diagm(eigdat.values)* S' - H)
LinearAlgebra.norm(S' * LinearAlgebra.diagm(eigdat.values) * S - H)







function Get_OP_schmidt(Op::SparseArrays.SparseMatrixCSC{ComplexF64,Int64})



