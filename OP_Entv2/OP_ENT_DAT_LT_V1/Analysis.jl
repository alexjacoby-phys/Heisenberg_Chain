import DelimitedFiles, LinearAlgebra, Plots

home = pwd()
L_RNG = 5:15


fnvec = [ string("L=",i,".txt") for i in L_RNG]

#fnvec = readdir("Data")
plot_data = zeros(Float64,length(fnvec))

cd("Data")

# for (lind,L) in pairs(L_RNG)
#     fn = fnvec[lind]
#     dat = DelimitedFiles.readdlm(fn)
#     println(fn)
# end
for (lind,L) in pairs(L_RNG)
    fn = fnvec[lind]
    dat = DelimitedFiles.readdlm(fn)

    size(dat, 1)
    svn_vec = [-LinearAlgebra.dot(dat[i, :], log.(dat[i, :])) for i in 1:size(dat, 1)]
    svn_avg = sum(svn_vec) / length(svn_vec)
    plot_data[lind] = svn_avg
end
using LaTeXStrings
Plots.plot(L_RNG,plot_data,legend = :none,yaxis = L"\mathcal{S}_{\rm vN}^{\rm Op.}, \ Half-Space",xaxis = L"{\rm System \ Size}", font = :times)

Plots.savefig("Vol_Law.pdf")


# size(dat,1)
# svn_vec = [ - LinearAlgebra.dot(dat[i,:], log.(dat[i,:])) for i in 1:size(dat,1)]

# sum(svn_vec)/length(svn_vec)
# for fn in fnvec
#     dat = DelimitedFiles.readdlm(fn)

# end
cd(home)