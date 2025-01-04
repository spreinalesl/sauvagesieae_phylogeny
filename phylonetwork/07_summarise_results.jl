using PhyloNetworks
using Gadfly
using PhyloPlots
using CSV
using DataFrames
using RCall

### Read snaq rus and plot the best network
net0 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/networks/snaq_sauvagesieae_h0.out")
rootatnode!(net0, "SRR10951426")
PhyloPlots.plot(net0, :R);

net1 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/networks/snaq_sauvagesieae_h1.out")
rootatnode!(net1, "SRR10951426")
PhyloPlots.plot(net1, :R, showGamma = true);

net2 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/networks/snaq_sauvagesieae_h2.out")
rootatnode!(net2, "SRR10951426")
PhyloPlots.plot(net2, :R, showGamma = true);

net3 = readSnaqNetwork("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/networks/snaq_sauvagesieae_h3.out")
rootatnode!(net3, "SRR10951426")
PhyloPlots.plot(net3, :R, showGamma = true);

### Select the best h and save it
#rootatnode!(net2, "panaque_cochliodon19170")
#R"pdf('net2_best.pdf', width=25, height=20)"
#PhyloPlots.plot(net2, :R, showGamma = true);
#R"dev.off()"

#net2_svg = SVG("net2.svg", 6inch, 4inch)
#draw(net2_svg, net2_plot)

### loglik scores for the best nets per iteration
scores = [net0.loglik, net1.loglik, net2.loglik, net3.loglik]
scores_plot = Gadfly.plot(x=collect(0:3), y=scores, Geom.point, Geom.line,
                          Guide.xlabel("h"), Guide.ylabel("Pseudo-deviance"))

scores_svg = SVG("scores.svg", 6inch, 4inch)
draw(scores_svg, scores_plot)


##################

### Explore alternative networks
netsh1 = readMultiTopology("/home/spreinalesl/Dropbox/Sandra/Manuscritos/NGS_phyloOchn/analyses/results/phylonetworks/networks/snaq_sauvagesieae_h1.networks")
