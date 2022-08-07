library(tidyverse)
library(tictoc)
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("analysis","function_evorate_fitting.R"))
pkg.net=c('igraph' ,'network','CINNA', 'centiserve', 'linkcomm','sna','netrankr','tidygraph')
xfun::pkg_load2(pkg.net)

#### PICKLE3 ####
pickle3 = read_delim("~/Downloads/UniProtNormalizedTabular-default.txt")
head(pickle3)

hs_ppi = pickle3 %>% dplyr::select(InteractorA,InteractorB) %>% arrange(InteractorA)
pickle3.cent = network.centrality()
saveRDS(pickle3.cent,"prepare/human_pickle3_ppi_centrality.rds")


mean(pickle3.cent$cent_deg)
median(pickle3.cent$cent_deg)
quantile(pickle3.cent$cent_deg,seq(0,1,len=21))
n_distinct(pickle3.cent$ids)


median(table(hs_ppi$InteractorA))
mean(table(hs_ppi$InteractorA))
sum(hs_ppi$InteractorA==hs_ppi$InteractorB)


#### IntAct ####
library(PItools)

# load all interactions from IntAct FTP storage (https://www.ebi.ac.uk/intact/downloads)
intact_ppi = fullInteractome(taxid = 9606, database = "IntActFTP",
                        format = "tab27", within_species=T,
                        clean = TRUE, # parse into usable format (takes 5-10 minutes)
                        protein_only = F, # filter protein interactions
                        directory = "data/intact/") # keep data files inside R library,
# or specify your directory