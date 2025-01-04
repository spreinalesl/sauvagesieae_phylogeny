## load packages
using Functors
using DataFrames
using CSV
using PhyloNetworks

## fetch the quartets without the assumption that there are the same terminales in each tree
trees = readlines("/home/sreinale/scratch/private/phylonetworks/mbresults/lesstrees_allloci.newick")

## convert from String newick to HybridNetwork with Functors.fmap
trees = Functors.fmap(readTopology, trees)

## calculate quartets
q,t = countquartetsintrees(trees)
df = writeTableCF(q,t)

CSV.write("/home/sreinale/scratch/private/phylonetworks/mbresults/CFs_sauvagesieae.csv", df)

## exit julia
exit()
