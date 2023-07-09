import DelimitedFiles, LinearAlgebra, Plots
using LaTeXStrings

home = pwd()
L_RNG = 5:15


fnvec = [ string("L=",i,".txt") for i in L_RNG]

#fnvec = readdir("Data")


cd("Data")

# for (lind,L) in pairs(L_RNG)
#     fn = fnvec[lind]
#     dat = DelimitedFiles.readdlm(fn)
#     println(fn)
# end
begin
    plot_data = zeros(Float64, length(fnvec))
    for (lind,L) in pairs(L_RNG)
        fn = fnvec[lind]
        dat = DelimitedFiles.readdlm(fn)

        svn_vec = [-LinearAlgebra.dot(dat[i, :], log.(dat[i, :])) for i in 1:size(dat, 1)]
        svn_avg = sum(svn_vec) / length(svn_vec)
        plot_data[lind] = svn_avg
    end

    L_RNG_DIV_EVEN = ( Vector(L_RNG[2:2:11]) * 0.5) .^(-1)
    L_RNG_DIV_ODD = (Vector(L_RNG[1:2:11]))
    Plots.plot(L_RNG[2:2:11], plot_data[2:2:11] .* (L_RNG[2:2:11] .^ (-1)), label=L"{\rm Even}", yaxis=L"\mathcal{S}_{\rm \textbf{vN}}^{\rm \textbf{Op.}}, \ \textbf{Half\textrm{-}Space}", xaxis=L"{\rm \textbf{System \ Size}}", font=:times)
    Plots.plot!(L_RNG[1:2:11], plot_data[1:2:11] .* (L_RNG[1:2:11] .^ (-1)), label=L"{\rm Odd}")
end

Plots.savefig("Vol_Law_Over_L.pdf")



begin
    variance = zeros(Float64, length(fnvec))

    for (lind, L) in pairs(L_RNG)
        fn = fnvec[lind]
        dat = DelimitedFiles.readdlm(fnvec[lind])

        svn_preavg = [-LinearAlgebra.dot(dat[i, :], log.(dat[i, :])) for i in 1:size(dat, 1)]
        variance_avg = sum((svn_preavg - fill(plot_data[lind], length(svn_preavg))) .^ 2) / length(svn_preavg)
        variance[lind] = sqrt(variance_avg)
    end

    Plots.plot(L_RNG, -log.(variance .* (L_RNG .^(-1))), yaxis=L"- {\rm \textbf{Log}}\left[\sigma\left(\mathcal{S}_{\rm \textbf{vN}}^{\rm \textbf{Op.}}\right) / L \right] , \ \textbf{Half\textrm{-}Space}", xaxis=L"{\rm \textbf{System \ Size}}", font=:times, legend = :none, linewidth = 1.5, color = :red, yticks = :none)
end


Plots.savefig("Variance.pdf")




begin
    maxes = zeros(Float64, length(fnvec))

    for (lind, L) in pairs(L_RNG)
        fn = fnvec[lind]
        dat = DelimitedFiles.readdlm(fnvec[lind])

        svn_preavg = [-LinearAlgebra.dot(dat[i, :], log.(dat[i, :])) for i in 1:size(dat, 1)]
        maxes[lind] = maximum(svn_preavg)
    end

    L_RNG_DIV_EVEN = (Vector(L_RNG[2:2:11]) * 0.5) .^ (-1)
    L_RNG_DIV_ODD = (Vector(L_RNG[1:2:11]))
    Plots.plot(L_RNG[2:2:11], maxes[2:2:11] .* (L_RNG[2:2:11] .^ (-1)), label=L"{\rm Even}", yaxis=L"{ \rm \textbf{Max}} \ \ \mathcal{S}_{\rm \textbf{vN}}^{\rm \textbf{Op.}}/L, \ \textbf{Half\textrm{-}Space}", xaxis=L"{\rm \textbf{System \ Size}}", font=:times)
    Plots.plot!(L_RNG[1:2:11], maxes[1:2:11] .* (L_RNG[1:2:11] .^ (-1)), label=L"{\rm Odd}")

end


Plots.savefig("LOE_Maxes.pdf")



# size(dat,1)
# svn_vec = [ - LinearAlgebra.dot(dat[i,:], log.(dat[i,:])) for i in 1:size(dat,1)]

# sum(svn_vec)/length(svn_vec)
# for fn in fnvec
#     dat = DelimitedFiles.readdlm(fn)

# end
cd(home)