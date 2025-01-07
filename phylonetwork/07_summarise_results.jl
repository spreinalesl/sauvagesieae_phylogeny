using PhyloNetworks
using Gadfly
using PhyloPlots
using CSV
using DataFrames
using RCall

### Read snaq rus and plot the best network
net0 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h0.out")
rootatnode!(net0, "Ochna_membranacea")
PhyloPlots.plot(net0, :R);

net1 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h1.out")
rootatnode!(net1, "Ochna_membranacea")
PhyloPlots.plot(net1, :R, showGamma = true);

net2 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h2.out")
rootatnode!(net2, "Ochna_membranacea")
PhyloPlots.plot(net2, :R, showGamma = true);

net3 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h3.out")
rootatnode!(net3, "Ochna_membranacea")
PhyloPlots.plot(net3, :R, showGamma = true);

net4 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h4.out")
rootatnode!(net4, "Ochna_membranacea")
PhyloPlots.plot(net4, :R, showGamma = true);

### Select the best h and save it
#rootatnode!(net2, "panaque_cochliodon19170")
#R"pdf('net2_best.pdf', width=25, height=20)"
#PhyloPlots.plot(net2, :R, showGamma = true);
#R"dev.off()"

#net2_svg = SVG("net2.svg", 6inch, 4inch)
#draw(net2_svg, net2_plot)

### loglik scores for the best nets per iteration
scores = [net0.loglik, net1.loglik, net2.loglik, net3.loglik, net4.loglik]
scores_plot = Gadfly.plot(x=collect(0:4), y=scores, Geom.point, Geom.line,
                          Guide.xlabel("h"), Guide.ylabel("Pseudo-deviance"))

scores_svg = SVG("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/scores.svg", 6inch, 4inch)
draw(scores_svg, scores_plot)

##########################################
############ Additional plots ############

### Explore alternative networks from .networks file
## This file contains the network with the best Ploglikelihood and a set of alternative networks after some topological perturbations
## on the best network
netsh1 = readMultiTopology("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h1.networks")
rootatnode!(netsh1[3], "Ochna_membranacea")
PhyloPlots.plot(netsh1[1], :R, showGamma = true);

### Explore alternative networks from .out file
## Before plotting networks delete the Dendroscope net and the whithe space before each net
netsh1 = readMultiTopology("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h1_ed.out")
PhyloPlots.plot(netsh1[1], :R, showGamma = true);
rootatnode!(netsh1[20], "Ochna_membranacea")

# To read one topology
addnet2 = readTopology("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/attempt_iqtree/snaq/snaq_sauvagesieae_h2_ed.out")
rootatnode!(addnet2, "Ochna_membranacea")
PhyloPlots.plot(addnet2, :R, showGamma = true);


