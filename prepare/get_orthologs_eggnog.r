source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
URL_FASTA_EGGNOG = "http://eggnogapi5.embl.de/nog_data/text/fasta"

# Fungi = 4751

fuNOG = get_eggnog_node(node = 4751) %>% mutate(url_fasta = sprintf("%s/%s",URL_FASTA_EGGNOG,OG) )
fu_species = get_eggnog_species(node = 4751)

fu_yeast = eggnog_annotations_species(node = 4751, species = c(4932,4896))
fu_taxons = count_taxons_eggnog_node(4751)

fu_dir = here::here('data','eggnog','4751-fuNOG')
fu_urls =
fu_fastafiles = sprintf('%s/%s.fasta',fu_dir,fuNOG$OG)
safely(download_file(fu_urls, fu_fastafiles)
fu_fasta=load_eggnog_fasta()

#fu_ali = get_eggnog_alignment(node=4751, use_trimmed = F)


#eggnog_fasta = Biostrings::readAAStringSet("http://eggnog5.embl.de/download/latest/e5.proteomes.faa")

# Mammalia = 40674

maNOG = get_eggnog_node(40674)
ma_species = get_eggnog_species(node = 40674)

ma_human = eggnog_annotations_species(node = 40674, species = c(9606))
ma_taxons = count_taxons_eggnog_node(40674)

#ma_ali = get_eggnog_alignment(node=40674, use_trimmed = F)

ma_fasta = load_eggnog_fasta(maNOG$OG)
