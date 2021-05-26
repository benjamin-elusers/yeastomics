library(rotl)
library(taxize)
#remotes::install_github("lindeloev/job")
library(job)
library(tidyverse)

my_git_repo="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main"
source(file.path(my_git_repo,"src/utils.r"))
source(file.path(my_git_repo,"src/function_annotation.r"))
source(file.path(my_git_repo,"src/function_sequence.r"))
source(file.path(my_git_repo,"src/function_phylogenetic.r"))
source(file.path(my_git_repo,"src/function_analysis.r"))
source(file.path(my_git_repo,"src/function_datalocal.r"))
source(file.path(my_git_repo,"src/function_datapub.r"))
# turn off annoying messages from dplyr::summarise
options(dplyr.summarise.inform = FALSE,dplyr.width=Inf)

# load fungi taxons data from gtrnadb
gtrna_fungal = rio::import("data/GtRNAdb-fungal-genomes.txt") %>%
separate(tax_ID, sep=":", c("taxid","lineage"))

# load fungi taxons data from eggnog
enog.fungi =load.eggnog.node(node="4751")
n_distinct(enog.fungi$taxid)
enog.tax =sort(unique(enog.fungi$taxid))

# compare fungi taxons from gtrnadb and eggnog
sum(enog.tax %in% gtrna_fungal$taxid)

# look for children of eggnog taxons
job::job(enog.child = {taxize::children(enog.tax, "ncbi")})
job::job(enog.down = {taxize::downstream(enog.tax,'ncbi',"no rank",intermediate=T)}, title='retrieve downstream eggnox taxons')


tmp = bind_rows( enog.child)
class(enog.child)
sapply(enog.down,dim)
getCDS(
  db = "refseq",
  organism,
  reference = FALSE,
  release = NULL,
  gunzip = FALSE,
  path = file.path("_ncbi_downloads", "CDS")
)

#taxize::id2name(id = names(trna.count),db='ncbi')
#fungis = stringr::str_split(unique(df.ortho$members_tax_ids)," ",simplify = T)[1,]
#ORTHO=load.eggnog.node('4890')
#n_distinct(ORTHO$taxid)
#trna.sp[!trna.sp %in% fungis]

#test= eggnog_taxons %>% filter(taxid %in% unique(ORTHO$taxid) ) %>% dplyr::select(taxid,sci_name)
#readr::write_tsv(test, file='~/Desktop/4890_Ascomycota.tax.name.tsv')
#head(ORTHO)
