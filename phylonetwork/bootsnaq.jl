using Distributed
using PhyloNetworks
using Functors
using DataFrames
using CSV

# use arguments for feeding tree and cfs
nruns = parse(Int, ARGS[1])
seed = parse(Int, ARGS[2])

# set up the max number of threads to use based on nruns
addprocs(nruns)

# Go to input path
cd("/home/sreinale/scratch/private/phylonetworks/iqtree")

# read the posterior trees for all loci
bt = readBootstrapTrees("/home/sreinale/scratch/private/phylonetworks/iqtree/gts_bootstrap.txt");

# read starting network
net = readTopology("/home/sreinale/scratch/private/phylonetworks/snaq/snaq_sauvagesieae_h3.out");

# calculate the bootstrap with a single replicate
bootnet = bootsnaq(net, bt, hmax=3, nrep=1, runs=nruns, filename="$seed"*"_bootsnaq_h3", seed=seed)

# save the bootstrap output just to make sure it is correct
writeMultiTopology(bootnet, /home/sreinale/scratch/private/phylonetworks/bootsnaq/"$seed"*"_bootstrapNets.tre")

#exit julia
exit()
