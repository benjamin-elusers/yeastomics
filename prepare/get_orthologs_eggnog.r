source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)

# Fungi = 4751
fu_dir = here::here('data','eggnog','4751-fuNOG')

fuNOG = get_eggnog_node(node = 4751)
fu_species = get_eggnog_species(node = 4751)

fu_yeast = eggnog_annotations_species(node = 4751, species = c(4932,4896))
fu_taxons = count_taxons_eggnog_node(4751)

fu_fastafiles = sprintf('%s/%s.fasta',fu_dir,fuNOG$OG)
download.file(fuNOG$url_fasta, fu_fastafiles, quiet=T,mode = 'wb')



walk2(fuNOG$url_fasta, fu_fastafiles, safe_download)

pbmcapply::pbmclapply(fuNOG$url_fasta, safe_download, mc.cores=14, path=fu_dir)

#fu_fasta=load_eggnog_fasta()
#fu_ali = get_eggnog_alignment(node=4751, use_trimmed = F)


#eggnog_fasta = Biostrings::readAAStringSet("http://eggnog5.embl.de/download/latest/e5.proteomes.faa")

# Mammalia = 40674

maNOG = get_eggnog_node(40674)
ma_species = get_eggnog_species(node = 40674)

ma_human = eggnog_annotations_species(node = 40674, species = c(9606))
ma_taxons = count_taxons_eggnog_node(40674)

#ma_ali = get_eggnog_alignment(node=40674, use_trimmed = F)

ma_fasta = load_eggnog_fasta(maNOG$OG)
