import DelimitedFiles, LinearAlgebra, Plots
using LaTeXStrings

home = pwd()
L_RNG = 4:17


fnvec = [ string("L=",i,".txt") for i in L_RNG]

#fnvec = readdir("Data")
home = "/Users/alexjacoby/Documents/Research_Code/Heisenberg_Chain/Full_Data"
fob = "/Users/alexjacoby/Documents/Research_Code/Heisenberg_Chain/Full_Data/Data"
cd(fob)

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

    # L_RNG_DIV_EVEN = ( Vector(L_RNG[2:2:12]) * 0.5) .^(-1)
    # L_RNG_DIV_ODD = (Vector(L_RNG[1:2:12]))
    Plots.plot(L_RNG[1:2:14], plot_data[1:2:14] .* (L_RNG[1:2:14] .^ (-1)), label=L"{\rm Even}", yaxis=L"\mathcal{S}_{\rm \textbf{vN}}^{\rm \textbf{Op.}}/L, \ \textbf{Half\textrm{-}Space}", xaxis=L"{\rm \textbf{System \ Size}}", font=:times)
    Plots.plot!(L_RNG[2:2:14], plot_data[2:2:14] .* (L_RNG[2:2:14] .^ (-1)), label=L"{\rm Odd}",ylims = (0,1))
end
cd(home)
Plots.savefig("Vol_Law_Over_L.pdf")
cd(fob)




begin
    plot_data = zeros(Float64, length(fnvec))
    for (lind,L) in pairs(L_RNG)
        fn = fnvec[lind]
        dat = DelimitedFiles.readdlm(fn)

        svn_vec = [-LinearAlgebra.dot(dat[i, :], log.(dat[i, :])) for i in 1:size(dat, 1)]
        svn_avg = sum(svn_vec) / length(svn_vec)
        plot_data[lind] = svn_avg
    end

    # L_RNG_DIV_EVEN = ( Vector(L_RNG[2:2:12]) * 0.5) .^(-1)
    # L_RNG_DIV_ODD = (Vector(L_RNG[1:2:12]))
    Plots.plot(L_RNG[1:2:14], plot_data[1:2:14] , label=L"{\rm Even}", yaxis=L"\mathcal{S}_{\rm \textbf{vN}}^{\rm \textbf{Op.}}, \ \textbf{Half\textrm{-}Space}", xaxis=L"{\rm \textbf{System \ Size}}", font=:times)
    Plots.plot!(L_RNG[2:2:14], plot_data[2:2:14], label=L"{\rm Odd}")
end
cd(home)
Plots.savefig("Vol_Law.pdf")
cd(fob)


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

cd(home)
Plots.savefig("Variance.pdf")
cd(fob)



begin
    maxes = zeros(Float64, length(fnvec))

    for (lind, L) in pairs(L_RNG)
        fn = fnvec[lind]
        dat = DelimitedFiles.readdlm(fnvec[lind])

        svn_preavg = [-LinearAlgebra.dot(dat[i, :], log.(dat[i, :])) for i in 1:size(dat, 1)]
        maxes[lind] = maximum(svn_preavg)
    end

    L_RNG_DIV_EVEN = (Vector(L_RNG[2:2:12]) * 0.5) .^ (-1)
    L_RNG_DIV_ODD = (Vector(L_RNG[1:2:12]))
    Plots.plot(L_RNG[1:2:14], maxes[1:2:14] .* (L_RNG[1:2:14] .^ (-1)), label=L"{\rm Even}", yaxis=L"{ \rm \textbf{Max}} \ \ \mathcal{S}_{\rm \textbf{vN}}^{\rm \textbf{Op.}}/L, \ \textbf{Half\textrm{-}Space}", xaxis=L"{\rm \textbf{System \ Size}}", font=:times)
    Plots.plot!(L_RNG[2:2:14], maxes[2:2:14] .* (L_RNG[2:2:14] .^ (-1)), label=L"{\rm Odd}")

end

cd(home)
Plots.savefig("LOE_Maxes.pdf")
cd(fob)


# size(dat,1)
# svn_vec = [ - LinearAlgebra.dot(dat[i,:], log.(dat[i,:])) for i in 1:size(dat,1)]

# sum(svn_vec)/length(svn_vec)
# for fn in fnvec
#     dat = DelimitedFiles.readdlm(fn)

# end
cd(home)










plot_data = zeros(Float64, length(fnvec))
for (lind, L) in pairs(L_RNG)
    fn = fnvec[lind]
    dat = DelimitedFiles.readdlm(fn)

    svn_vec = [-LinearAlgebra.dot(dat[i, :], log.(dat[i, :])) for i in 1:size(dat, 1)]
    svn_avg = sum(svn_vec) / length(svn_vec)
    plot_data[lind] = svn_avg
end



plot_data

plot_data[3:2:13] - plot_data[1:2:11] 
plot_data[4:2:14] - plot_data[2:2:12]


# L_RNG_DIV_EVEN = ( Vector(L_RNG[2:2:12]) * 0.5) .^(-1)
# L_RNG_DIV_ODD = (Vector(L_RNG[1:2:12]))
Plots.plot(L_RNG[1:2:14], plot_data[1:2:14] .* (L_RNG[1:2:14] .^ (-1)), label=L"{\rm Even}", yaxis=L"\mathcal{S}_{\rm \textbf{vN}}^{\rm \textbf{Op.}}/L, \ \textbf{Half\textrm{-}Space}", xaxis=L"{\rm \textbf{System \ Size}}", font=:times)
Plots.plot!(L_RNG[2:2:14], plot_data[2:2:14] .* (L_RNG[2:2:14] .^ (-1)), label=L"{\rm Odd}", ylims=(0, 1))