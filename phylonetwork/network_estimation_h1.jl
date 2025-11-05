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

cfs = readTableCF("CFs_sauvagesieae.csv")

# read the h=0 network
net0 = readTopology("snaq_sauvagesieae_h0.out")

# calculate a h=1 network
net1 = snaq!(net0, cfs, hmax=1, filename="snaq_sauvagesieae_h1", seed=1324, runs = nruns)

#exit julia
exit()
