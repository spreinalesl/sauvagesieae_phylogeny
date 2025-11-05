using PhyloNetworks
using Gadfly
using PhyloPlots
using CSV
using DataFrames
using RCall

#############################################
### Read snaq rus and plot the best network
### use absolute paths
#############################################
net0 = readSnaqNetwork("snaq_sauvagesieae_h0.out")
rootatnode!(net0, "Ochna_membranacea")
PhyloPlots.plot(net0, :R);

net1 = readSnaqNetwork("snaq_sauvagesieae_h1.out")
rootatnode!(net1, "Ochna_membranacea")
PhyloPlots.plot(net1, :R, showGamma = true);

net2 = readSnaqNetwork("snaq_sauvagesieae_h2.out")
rootatnode!(net2, "Ochna_membranacea")
PhyloPlots.plot(net2, :R, showGamma = true);

net3 = readSnaqNetwork("snaq_sauvagesieae_h3.out")
rootatnode!(net3, "Ochna_membranacea")
PhyloPlots.plot(net3, :R, showGamma = true);

net4 = readSnaqNetwork("snaq_sauvagesieae_h4.out")
rootatnode!(net4, "Ochna_membranacea")
PhyloPlots.plot(net4, :R, showGamma = true);

net5 = readSnaqNetwork("snaq_sauvagesieae_h5.out")
rootatnode!(net5, "Ochna_membranacea")
PhyloPlots.plot(net5, :R, showGamma = true);

#############################################
### Get loglik scores for the best nets
### per run/iteration
#############################################
scores = [net0.loglik, net1.loglik, net2.loglik, net3.loglik, net4.loglik, net5.loglik]

#############################################
### Plot results
#############################################
### Networks
## pdf
R"pdf('net2_best.pdf',
       width=25, height=20)"
PhyloPlots.plot(net2, :R, showGamma = true);
R"dev.off()"

## svg
net2_svg = SVG("net2.svg", 6inch, 4inch)
draw(net2_svg, net2_plot)

### Scores
scores_plot = Gadfly.plot(x=collect(0:5), y=scores, Geom.point, Geom.line,
                          Guide.xlabel("h"), Guide.ylabel("Pseudo-deviance"))

scores_svg = SVG("scores.svg", 6inch, 4inch)
draw(scores_svg, scores_plot)

##########################################
############ Additional plots ############

### Explore alternative networks from .networks file
## This file contains the network with the best Ploglikelihood and a set of alternative networks after some topological perturbations
## on the best network
netsh1 = readMultiTopology("snaq_sauvagesieae_h1.networks")
rootatnode!(netsh1[3], "Ochna_membranacea")
PhyloPlots.plot(netsh1[1], :R, showGamma = true);

### Explore alternative networks from .out file
## Before plotting networks delete the Dendroscope net and the whithe space before each net
netsh1 = readMultiTopology("snaq_sauvagesieae_h1_ed.out")
PhyloPlots.plot(netsh1[1], :R, showGamma = true);
rootatnode!(netsh1[20], "Ochna_membranacea")

# To read one topology
addnet2 = readTopology("snaq_sauvagesieae_h2_ed.out")
rootatnode!(addnet2, "Ochna_membranacea")
PhyloPlots.plot(addnet2, :R, showGamma = true);
