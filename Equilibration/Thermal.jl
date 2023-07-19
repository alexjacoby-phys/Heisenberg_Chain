include("Functions.jl")

L=12

opinit = LinearAlgebra.normalize(kron(coupling(1, L, spZ, spZ, L)...) +kron(coupling(1, L, spY, spY, L)...)+ kron(coupling(1, L, spX, spX, L)...))



H = make_H(L)

eigdat = LinearAlgebra.eigen(Matrix(H))
Evec = eigdat.values
S = eigdat.vectors


distribution = [LinearAlgebra.norm(kron(S[:, k], S[:, k]') * opinit ) for k in 1:size(S, 1)] * (2^L)

import Plots


Plots.plot(Evec,distribution,seriestype = :scatter)


