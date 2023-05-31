include("qi_tools.jl")
import LinearAlgebra

L = 10
psi = LinearAlgebra.normalize(rand(ComplexF64, 2^L))



A = 5
B = 5
rhonew = reshape(psi, (2^A, 2^B))

info = LinearAlgebra.svd(rhonew)


SvN_Small(Matrix{ComplexF64}(LinearAlgebra.diagm(info.S .^ 2)))


rho = Partial_Trace(kron(psi, psi'), 0, 5, 5)
SvN_Small(rho)
log(2)





A .* B
import SparseArrays
SparseArrays.spzeros(10)


dim = 2
svn_mat = reshape(psi, dim^A, dim^B)
spectrum = LinearAlgebra.svdvals(svn_mat) .^ 2
svn_vec = spectrum .* log.(spectrum)
-sum(svn_vec)
ϵ = 10e-015

LinearAlgebra.fill(ϵ,10)

min(3,10)

"Takes a state vector and makes bipartition according to A, B. The default on site dimension is set to two but can be changed freely. Furthermore, an ϵ is included to ensure convergence with a default value of 10^-14. You must have loaded linear algebra for this to work"
function SvN(psi::Vector{ComplexF64}, A::Int64, B::Int64, dim::Int64=2, ϵ::Float64=10e-15)
    svn_mat = reshape(psi, dim^A, dim^B)
    spectrum = LinearAlgebra.svdvals(svn_mat) .^ 2
    svn_vec = spectrum .* log.(spectrum + LinearAlgebra.fill(ϵ,dim^(min(A,B))))
    return -sum(svn_vec)
end
function SvN(psi::SparseArrays.SparseVector{ComplexF64, Int64}, A::Int64, B::Int64, dim::Int64=2, ϵ::Float64=10e-15)
    svn_mat = reshape(psi, dim^A, dim^B)
    spectrum = LinearAlgebra.svdvals(svn_mat) .^ 2
    svn_vec = spectrum .* log.(spectrum + LinearAlgebra.fill(ϵ, dim^(min(A, B))))
    return -sum(svn_vec)
end
