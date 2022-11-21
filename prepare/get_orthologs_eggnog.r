source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)

#eggnog_fasta = Biostrings::readAAStringSet("http://eggnog5.embl.de/download/latest/e5.proteomes.faa")
fungi_clades = c('4890'='ascomycota', '5204'='basidiomycota')
ascomycota_clades = c('451866'='taphrimycotina','4891'='saccharomycetes','147541'='dothideomycetes',
                      '147545'='eurotiomycetes','147550'='sordariomycetes','147548'='leotiomycetes')
metazoa_clades = c('7742'='vertebrata','40674'='mammalia','50557'='insecta','6231'='nematoda')
clades = lst('33154'=
          lst('4751'=lst('5204','4890'=ascomycota_clades)),
          lst('33208'=metazoa_clades)
         )

# Fungi = 4751
fu_dir = here::here('data','eggnog','4751-fungi')
fuNOG      = get_eggnog_node(node = 4751)
fu_tax     = get_eggnog_taxonomy(4751)
fu_species = get_eggnog_species(node = 4751)
fu_yeast   = eggnog_annotations_species(node = 4751, species = c(4932,4896))
fu_taxons  = count_taxons_eggnog_node(4751, subnode=c(4890,5204,451866,4891,147541,147545,147550,147548)  )


## Make the fungi species tree
### Select orthogroups with 179 species and at most 200 proteins ###


# Filter orthogroups fasta to keep at most 179 species (ortholog closest to yeast)



fu_fastafiles = sprintf('%s/%s.fasta',fu_dir,fuNOG$OG)
#pbmcapply::pbmclapply(fuNOG$url_fasta, safe_download, mc.cores=14, path=fu_dir, ext='.fasta')
#fu_ali = get_eggnog_alignment(node=4751, use_trimmed = F)

fu_clades = map_dfr(names(fungi_clades), ~get_eggnog_taxonomy(.x))
fu_4890 = count_eggnog_node(4751, 4890)
fu_5204 = count_eggnog_node(4751, 5204)


# Ascomycetes = 4890
am_dir = here::here('data','eggnog','4890-ascomycetes')
amNOG = get_eggnog_node(node = 4890)
am_tax=get_eggnog_taxonomy(4890)
am_species = get_eggnog_species(node = 4890)

am_yeast = eggnog_annotations_species(node = 4890, species = c(4932,4896))
am_taxons = count_taxons_eggnog_node(4890)

am_fastafiles = sprintf('%s/%s.fasta',am_dir,amNOG$OG)
#pbmcapply::pbmclapply(amNOG$url_fasta, safe_download, mc.cores=14, path=am_dir, ext='.fasta')
#am_ali = get_eggnog_alignment(node=4890, use_trimmed = F)

am_clades = map_dfr(names(ascomycota_clades), ~get_eggnog_taxonomy(.x))

am_451866 = count_eggnog_node(4890,451866)
am_4891 = count_eggnog_node(4890,4891)
am_147541 = count_eggnog_node(4890,147541)
am_147545 = count_eggnog_node(4890,147545)
am_147550 = count_eggnog_node(4890,147550)
am_147548 = count_eggnog_node(4890,147548)

# Mammalia = 40674
ma_dir = here::here('data','eggnog','40674-mammalia')

maNOG = get_eggnog_node(40674)
ma_tax=get_eggnog_taxonomy(40674)
ma_species = get_eggnog_species(node = 40674)

ma_human = eggnog_annotations_species(node = 40674, species = c(9606))
ma_taxons = count_taxons_eggnog_node(40674)

pbmcapply::pbmclapply(maNOG$url_fasta, safe_download, mc.cores=14, path=ma_dir, ext='.fasta')
#ma_ali = get_eggnog_alignment(node=40674, use_trimmed = F)

# Metazoa = 33208
mz_clades = str_split_fixed(metazoa_clades,pattern = '-',n=2) %>% as_tibble() %>% pull(2,1)
mz_dir = here::here('data','eggnog','33208-metazoa')

mzNOG = get_eggnog_node(33208)
mz_tax=get_eggnog_taxonomy(33208)
mz_species = get_eggnog_species(node = 33208)

mz_human = eggnog_annotations_species(node = 33208, species = c(9606))
mz_taxons = count_taxons_eggnog_node(33208)

pbmcapply::pbmclapply(mzNOG$url_fasta, safe_download, mc.cores=14, path=mz_dir, ext='.fasta')
#mz_ali = get_eggnog_alignment(node=33208, use_trimmed = F)
mz_clades = map(metazoa_clades, ~get_eggnog_taxonomy(.x))

get_eggnog_species(7742)
mz_7742 = count_eggnog_node(33208, 7742)
mz_40674 = count_eggnog_node(33208, 40674)
mz_50557 = count_eggnog_node(33208, 50557)
mz_6231 = count_eggnog_node(33208, 6231)
