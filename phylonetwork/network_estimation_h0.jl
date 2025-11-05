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

# read starting tree
start_tree = readTopology("startree.tre")

# calculate a h=0 network
net0 = snaq!(start_tree, cfs, hmax=0, filename="snaq_sauvagesieae_h0", seed=1234, runs = nruns)

#exit julia
exit()
