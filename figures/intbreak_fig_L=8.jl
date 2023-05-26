import Plots, DelimitedFiles


L = 8
T = 5.0
δ = 0.01
τ = 0.0:δ:T
disorder_strengths = [0.0, 0.05, 0.1, 0.25, 0.5]




dat = "/Users/alexjacoby/Documents/Research_Code/Heisenberg_Chain/Trial2/dat_Trial2.txt"



io = open(dat, "r")
data_table = DelimitedFiles.readdlm(dat)
close(io)


Labels = []
for k in disorder_strengths
    push!(Labels, string("Disorder Strength  = ", k))
end


GRAD = []
K = 5
for i in 0:K-1
    push!(GRAD, Plots.RGB(Float64(i / K), 0.0, 1- Float64(i / K)))
end

Plots.plot()
for i in 1:5
    Plots.plot!(τ, data_table[i, :], xlabel = "Time (No Units)", ylabel = "LOE", linewidth=6, color=GRAD[i], label=Labels[i], legend=:right, markersize=6, markeralpha=0.8) #=,seriestype = :scatter=#
end
Plots.plot!()
Plots.savefig("intbreak-L=8.pdf")




Plots.plot()
for i in 1:5
    Plots.plot!(τ, exp.(data_table[i, :]), xlabel = "Time (No Units)", ylabel = "LOE (Exponential)", linewidth=6, color=GRAD[i], label=Labels[i], legend=:right, markersize=6, markeralpha=0.8) #=,seriestype = :scatter=#
end
Plots.plot!()
Plots.savefig("intbreak-L=8_logtime.pdf")

