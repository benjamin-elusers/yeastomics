## \# Yeastomics

![yeastomics data](https://github.com/benjamin-elusers/yeastomics/blob/main/docs/imgs/YeastOmics-logo-ldpi.png "Yeastomics data")

(Image/Logos credits to [BioRender/Freepik\@Flaticon](mailto:BioRender/Freepik@Flaticon))

**Yeastomics** gather a set of tools and functions I use to analyze yeast proteome data.

> **SGD** is the reference for data related to *genome sequences*.
>
> **UniprotKB** is the reference for *data related to proteome*.

## ## Sequence

-   read proteomes from list of files (with parallel and/or progress bar)
-   read nucleic sequence
-   read protein sequence

## ## Annotations

-   retrieve annotations from controlled vocabulary (MI)
-   get annotations from Uniprot description (localizations,complex,roles,functions...)
-   get GeneOntology annotations (localizations)

## ## Alignment

-   make pairwise alignment
-   compute summary statistics about alignment (overlap, %ID, mis/matches...)
-   convert alignment to residues table

## ## Phylogenetic

-   read Rate4site (R4S) output
-   get tabular data about YGOB yeast ohnologs
-   get eggnog data about a particular taxonomic level
-   find common ancestors between lineages written as strings

## ## Resource databases *(remote access via URL)*

-   PFAM (domain assignment)
-   SUPERFAMILY (domain assignment)
-   UniprotKB (proteome sequences, features, subcellular localizations, families)
-   SGD (genome sequences, features, gene ontology annotations)
-   D2P2 (Consensus disorder predictions)
-   STRING (protein-protein interactions)
-   eggNog (orthology prediction and functional annotation)
-   paxDB (proteome abundance)

## ## Published datasets *(remote access via URL)*

\* proteome data not available (computed locally)

| Year   | 1st Author     | Journal                 | Title                                                                                                                           |
|--------|----------------|-------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| 2005   | Byrne          | Genome Research         | The Yeast Gene Order Browser: Combining curated homology and syntenic context reveals gene fate in polyploid species            |
| 2006   | Belle          | PNAS                    | Quantification of protein half-lives in the budding yeast proteome                                                              |
| 2008   | Pu             | NAR                     | Up-to-date catalogues of yeast protein complexes                                                                                |
| 2010   | Costanzo       | Science                 | The Genetic Landscape of a cell                                                                                                 |
| 2010   | Barton         | PLOS One                | Evolutionary Systems Biology of Amino Acid Biosynthetic Cost in Yeast                                                           |
| 2012   | Marguerat      | Cell                    | Quantitative Analysis of Fission Yeast Transcriptomes and Proteomes in Proliferating and Quiescent Cells                        |
| 2014   | Christiano     | Cell Reports            | Global Proteome Turnover Analyses of the Yeasts *S. cerevisiae* and *S. pombe*                                                  |
| 2014   | Lee            | Science                 | Mapping the Cellular Response to Small Molecules Using Chemogenomic Fitness Signatures                                          |
| 2014   | Geisberg       | Cell                    | Global Analysis of mRNA Isoform Half-Lives Reveals Stabilizing and Destabilizing Elements in Yeast                              |
| 2014   | Dana           | G3                      | Mean of the Typical Decoding Rates: A New Translation Efficiency Index Based on the Analysis of Ribosome Profiling Data         |
| 2014   | Lancaster      | Bioinformatics          | PLAAC: a web and command-line application to identify proteins with Prion-Like Amino Acid Composition                           |
| 2015   | Filleton       | Epigenetics & Chromatin | The complex pattern of epigenomic variation between natural yeast strains at single-nucleosome resolution                       |
| 2016   | Van Leeuwen    | Science                 | Exploring genetic suppression interactions on a global scale                                                                    |
| 2016\* | Bolognesi      | Cell Reports            | A concentration-dependent liquid phase separation can cause toxicity upon increased protein expression                          |
| 2017   | Lahtvee        | Cell Systems            | Absolute Quantification of Protein and mRNA Abundances Demonstrate Variability in Gene-Specific Translation Efficiency in Yeast |
| 2017   | Villen         | Cell Systems            | Determinants and Regulation of Protein Turnover in Yeast                                                                        |
| 2017   | Mittal         | Nature Communications   | The Gcn4 transcription factor reduces protein synthesis capacity and extends yeast lifespan                                     |
| 2017   | Leueunberger   | Science                 | Cell-wide analysis of protein thermal unfolding reveals determinants of thermostability                                         |
| 2018   | Ho             | Cell Systems            | Unification of Protein Abundance Datasets Yields a Quantitative Saccharomyces cerevisiae Proteome                               |
| 2018   | Peter          | Science                 | Genome evolution across 1,011 *Saccharomyces cerevisiae* isolates                                                               |
| 2018\* | Vernon         | eLife                   | Pi-Pi contacts are an overlooked protein feature relevant to phase separation                                                   |
| 2019   | Meldal         | NAR                     | Complex Portal 2018: extended content and enhanced visualization tools for macromolecular complexes                             |
| 2019   | Hausser        | Nature Communications   | Central dogma rates and the trade-off between precision and economy in gene expression                                          |
| 2019   | Dubreuil       | JMB                     | Protein Abundance Biases the Amino Acid Composition of Disordered Regions to Minimize Non-functional Interactions               |
| 2020   | Szavits-Nossan | Nucleic Acids Research  | Inferring efficiency of translation initiation and elongation from ribosome profiling                                           |
| 2020   | Van Leeuwen    | Mol. Sys. Bio.          | Systematic analysis of bypass suppression of essential genes                                                                    |
| 2021   | Dubreuil       | Frontiers Mol. Biosc    | Abundance imparts evolutionary constraints of similar magnitude on the buried, surface, and disordered regions of proteins      |

## ## Local datasets *(E. Levy inhouse data)*

-   3dcomplex *quaternary structures*
-   INTACT *curated protein interactions*
-   R4S data *for fungi lineage* & *for 1011 yeast strains*

## ## Analysis

-   make equally sized/distributed bins (1st and last can be set manually)
-   get values at decile x (`d*()` where \* = 10,40,60,90)
-   get values at quantile y (`q*()` where \* = 25,75)
-   filter data based on deciles/quantile values at x (with/out ties)
-   get values in interquantile range (`get.iqr()`)
-   slice tibble data for specific quantile or interquantile range (`slice_d*()` and `slice_q*()`)
-   compute network centrality measures (`network.centrality()`)
-   discretize variables into 5 custom bins (`fivebins()`)
-   get extreme (outliers) values/index/boundaries from a vector (`get_extremes()` ,`get_outliers`,`get_outliers_index`,`get_outliers_boundary()`)
-   compute protein amino acid scores from residues count (`AACOUNT2SCORE()`)
-   detect rows with NAs (`find_na_rows()`)
-   coalesce-join between dataframe (NA are replaced in shared columns for which one of the value is defined `coalesce_join()`)
-   find closest value in vector (`nearest()`)
-   Toolbox for regression and unidimensional analysis...

## ## Utilities

-   String manipulations
-   Shorthands and wrappers to useful functions
-   Environment and memory
-   Counting elements
-   Filtering/Subsetting objects
-   Testing/Checking values
-   I/O
-   Values formatting
-   Precomputed-data (Sequences/Amino-acids/AA scores)
