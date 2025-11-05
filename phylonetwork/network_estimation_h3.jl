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

# read the h=2 network
net2 = readTopology("snaq_sauvagesieae_h2.out")

# calculate a h=3 network
net3 = snaq!(net2, cfs, hmax=3, filename="snaq_sauvagesieae_h3", seed=4213, runs = nruns)

#exit julia
exit()
