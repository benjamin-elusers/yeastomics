# Yeastomics
-------------

![yeastomics data](https://github.com/benjamin-elusers/yeastomics/blob/main/docs/imgs/YeastOmics-logo.png "Yeastomics data"| width=100)
(Image/Logos credits to BioRender/Freepik@Flaticon)

Yeastomics gather a set of tools and functions I use to analyze yeast proteome data.

**SGD** is the reference for data related to _genome sequences_.

**UniprotKB** is the reference for _data related to proteome_.


## Sequence
-------------

- read proteomes from list of files
- read nucleic sequence 
- read protein sequence

## Alignment
-------------

make pairwise alignment

## Phylogenetic
-------------

read Rate4site (R4S) output


## Resource databases _(remote access via URL)_
-------------

- PFAM (domain assignment)
- SUPERFAMILY (domain assignment)
- UniprotKB (proteome sequences, features, subcellular localizations, families)
- SGD (genome sequences, features, gene ontology annotations)
- D2P2 (Consensus disorder predictions)
- STRING (protein-protein interactions)
- eggNog (orthology prediction and functional annotation)
- paxDB (proteome abundance)

## Published datasets _(remote access via URL)_
-------------
|Year|1st Author  |Journal             |Title                                                                                                                     |
|----|------------|--------------------|--------------------------------------------------------------------------------------------------------------------------|
|2005|Byrne       |Genome Research     |The Yeast Gene Order Browser: Combining curated homology and syntenic context reveals gene fate in polyploid species      |
|2006|Belle       |PNAS                |Quantification of protein half-lives in the budding yeast proteome                                                        |
|2008|Pu          |NAR                 |Up-to-date catalogues of yeast protein complexes                                                                          |
|2010|Costanzo    |Science             |The Genetic Landscape of a cell                                                                                           |
|2014|Lee         |Science             |Mapping the Cellular Response to Small Molecules Using Chemogenomic Fitness Signatures                                    |
|2014|Geisberg    |Cell                |Global Analysis of mRNA Isoform Half-Lives Reveals Stabilizing and Destabilizing Elements in Yeast                        |
|2016|Van Leeuwen |Science             |Exploring genetic suppression interactions on a global scale                                                              |
|2017|Villen      |Cell Systems        |Determinants and Regulation of Protein Turnover in Yeast                                                                  |
|2017|Leueunberger|Science             |Cell-wide analysis of protein thermal unfolding reveals determinants of thermostability                                   |
|2018|Ho          |Cell Systems        |Unification of Protein Abundance Datasets Yields a Quantitative Saccharomyces cerevisiae Proteome                         |
|2018|Peter       |Science             |Genome evolution across 1,011 Saccharomyces cerevisiae isolates                                                           |
|2019|Meldal      |NAR                 |Complex Portal 2018: extended content and enhanced visualization tools for macromolecular complexes                       |
|2019|Dubreuil    |JMB                 |Protein Abundance Biases the Amino Acid Composition of Disordered Regions to Minimize Non-functional Interactions         |
|2020|Van Leeuwen |Mol. Sys. Bio.      |Systematic analysis of bypass suppression of essential genes                                                              |
|2021|Dubreuil    |Frontiers Mol. Biosc|Abundance imparts evolutionary constraints of similar magnitude on the buried, surface, and disordered regions of proteins|

## Local datasets _(E. Levy inhouse data)_
-------------

- 3dcomplex *quaternary structures*
- INTACT *curated protein interactions*
- R4S data *for fungi lineage* & *for 1011 yeast strains*

## Analysis
-------------

make equal bins (1st and last bin manual)
compute network centrality
discretize variables according to quantiles

## Utilities
-------------

- String manipulations
- Shorthands
- Environment and memory
- Counting
- Filtering/Subsetting
- Precomputed-data

