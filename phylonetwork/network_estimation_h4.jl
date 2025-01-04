using Distributed

# how many runs for snaq
nruns = 20

# set up the max number of threads to use based on nruns
addprocs(nruns)

# load packages to parallel threads
@everywhere using Functors
@everywhere using DataFrames
@everywhere using CSV
@everywhere using PhyloNetworks

cfs = readTableCF("/home/sreinale/scratch/private/phylonetworks/snaq/CFs_sauvagesieae.csv")

# read the h=3 network
net3 = readTopology("/home/sreinale/scratch/private/phylonetworks/snaq/snaq_sauvagesieae_h3.out")

# calculate a h=4 network
net4 = snaq!(net3, cfs, hmax=4, filename="/home/sreinale/scratch/private/phylonetworks/snaq/snaq_sauvagesieae_h4", seed=2431, runs = nruns)

#exit julia
exit()
