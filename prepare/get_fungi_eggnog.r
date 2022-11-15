source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)

fuNOG = get_eggnog_node(node = 4751)
fu_species = get_eggnog_species(node = 4751)

fu_yeast = eggnog_annotations_species(node = 4751, species = c(4932,4896))
fu_taxons = count_taxons_eggnog_node(4751)

fu_ali = get_eggnog_alignment(node=4751, use_trimmed = F)
Biostrings::AAStringSet(fu_ali$`3NTVE`)

eggnog_fasta = Biostrings::readAAStringSet("http://eggnog5.embl.de/download/latest/e5.proteomes.faa")
