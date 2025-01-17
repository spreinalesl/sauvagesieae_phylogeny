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

# read the posterior trees for all loci
bt = readBootstrapTrees("~/scratch/private/phylonetworks/iqtree/lesstrees_allloci.newick");

# read starting network
net = readTopology("~/scratch/private/phylonetworks/snaq/snaq_sauvagesieae_h3.out");

# calculate the bootstrap with a single replicate
bootnet = bootsnaq(net, bt, hmax=3, nrep=1, runs=nruns, filename="$seed"*"_bootsnaq_h3", seed=seed)

# save the bootstrap output just to make sure it is correct
writeMultiTopology(bootnet, ~/scratch/private/phylonetworks/bootsnaq/"$seed"*"_bootstrapNets.tre")

#exit julia
exit()
