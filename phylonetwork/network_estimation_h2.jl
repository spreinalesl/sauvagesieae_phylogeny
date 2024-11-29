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

cfs = readTableCF("/home/sreinale/scratch/private/phylonetworks/mbresults/CFs_sauvagesieae_ed.csv")

# read the h=1 network
net1 = readTopology("/home/sreinale/scratch/private/phylonetworks/pnetwork/snaq_sauvagesieae_h1.out")

# calculate a h=2 network
net2 = snaq!(net1, cfs, hmax=2, filename="/home/sreinale/scratch/private/phylonetworks/pnetwork/snaq_sauvagesieae_h2", seed=1423, runs = nruns)

#exit julia
exit()
