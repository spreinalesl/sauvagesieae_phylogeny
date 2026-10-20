# Sauvagesieae phylogenomics

In this study, we evaluated a methodological framework for analysing target-capture data obtained using the dual-hybridisation approach to integrate both the universal Angioperms353 and the Ochnaceae-specific
probe sets on the tribe Sauvagesieae (Ochnaceae). 

We also assessed the impact of multiple data-filtering strategies on phylogenetic reconstruction, including handling of fragmentary sequences and paralogs. Then, we used the best filtered
dataset to test for molecular evolutionary processes within this lineage, such as hybridisation and recent polyploidisation. 

This directory contains custom scripts used to download already-published data, clean, select and assemble DNA reads from target-capture. Scripts to detect and remove shared genes between probe sets are included as well. Read mapping and assembly were performed using the HybPiper pipeline, paralog assessments was performed using HybPiper, HybPhaser and a custom strategy. 

Phylogenetic inference was carried out using a species tree approach with Weighted ASTRAL and ASTRAL-Pro depending on the nature of data used to infer gene trees. 

Additional analyses to test reticulated evolution using SNaQ and recent polyploidisation using nQuack were also performed and the corresponding scripts are available. Scripts to produce summary tables and figures included in the paper are also provided.

This repository can be used for academic and research proposes as long as the corresponding paper is properly cited:

Reinales S, Forest F, Zuntini AR, Cardoso D, Ballen GA, Cárdenas D, Pirani JR. Integrating lineage-specific and universal genomic probes illuminates phylogenetic relationships and molecular evolution in Sauvagesieae (Ochnaceae). Molecular Phylogenetics and Evolution. 2026 Jul 14:108687.

If you use bibtex, use the following entry in your .bib file:

```
@article{Reinales2026sauv,
  title={{Integrating lineage-specific and universal genomic probes illuminates phylogenetic relationships and molecular evolution in Sauvagesieae (Ochnaceae)}},
  author={Sandra Reinales and Félix Forest and Alexandre R. Zuntini and Domingos Cardoso and Gustavo A. Ballen and Dairon Cárdenas and José R. Pirani},
  journal={Molecular Phylogenetics and Evolution},
  volume={224},
  pages={108687},
  year={2026},
  doi={10.1016/j.ympev.2026.108687},
  issn={1055-7903},
  url={https://www.sciencedirect.com/science/article/pii/S1055790326001570}
}
```
