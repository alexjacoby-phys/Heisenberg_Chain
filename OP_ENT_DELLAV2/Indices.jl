# using Pkg
# Pkg.activate("jj8144/home/.julia/LOE_Code")
# Pkg.instantiate()

include("Functions.jl")


Lmin = 4
Lmax = 8
L_RNG = Lmin:1:Lmax


fob = "/Users/alexjacoby/Documents/Research_Code/Heisenberg_Chain/OP_ENT_DELLAV2/INDEX_DATA/"
cd(fob)
# L=10
# A = L ÷ 2
# B = L - A
# Q = L ÷ 2
# partitions = Charge_Partitions(L, Q, A, B)
# partition_double = reshape(collect(Iterators.product(partitions,partitions)),length(partitions)^2)



# touch("testfile.txt")
# DelimitedFiles.writedlm("testfile.txt",partitions)

f(x) = (x[1]..., x[2]...)
for L in L_RNG
    fn = string("L=", L)
    mkdir(fn)
    cd(fn)
    A = L ÷ 2
    B = L - A
    Q = L ÷ 2
    partitions = Charge_Partitions(L, Q, A, B)
    pre_partition_double = reshape(collect(Iterators.product(partitions, partitions)), length(partitions)^2)

    partition_double = f.(pre_partition_double)

    cvec = basis_converter(L, Q)

    touch("Partition_Double.txt")
    DelimitedFiles.writedlm("Partition_Double.txt", partition_double)
    touch("QBasis.txt")
    DelimitedFiles.writedlm("QBasis.txt", cvec)
    for (p_no, (QA, QB, QAtilde, QBtilde)) in pairs(partition_double)
        CTC = QtoSVD_Basis_Tuple(L, Q, A, B, QA, QAtilde)
        CTCL = [t[1] for t in CTC]
        CTCR = [t[2] for t in CTC]
        CTC_L_FILE = string("L_Partition_number=", p_no, ".txt")
        CTC_R_FILE = string("R_Partition_number=", p_no, ".txt")
        touch(CTC_L_FILE)
        touch(CTC_R_FILE)
        DelimitedFiles.writedlm(CTC_L_FILE, CTCL)
        DelimitedFiles.writedlm(CTC_R_FILE, CTCR)
    end
    cd(fob)
end


