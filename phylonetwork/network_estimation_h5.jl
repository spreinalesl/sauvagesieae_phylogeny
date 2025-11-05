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

# read the h=4 network
net4 = readTopology("snaq_sauvagesieae_h4.out")

# calculate a h=5 network
net5 = snaq!(net4, cfs, hmax=5, filename="snaq_sauvagesieae_h5", seed=5432, runs = nruns)

#exit julia
exit()
