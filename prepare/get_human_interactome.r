library(tidyverse)
library(tictoc)
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("analysis","function_evorate_fitting.R"))
pkg.net=c('igraph' ,'network','CINNA', 'centiserve', 'linkcomm','sna','netrankr','tidygraph')
xfun::pkg_load2(pkg.net)

#### PICKLE3 ####
URL_HS_PICKLE3="http://www.pickle.gr/Data/3.3/PICKLE_HUMAN_3_3_UniProtNormalizedTabular-default.zip"
pickle3 = rio::import(URL_HS_PICKLE3)
hs_pickle3 = pickle3 %>% dplyr::select(InteractorA,InteractorB) %>% arrange(InteractorA)
pickle3.cent = preload(here::here("prepare","human_pickle3_ppi_centrality.rds"),
                       network.centrality(hs_pickle3),
                       'compute centrality of human interactome from PICKLE3...')
# statistics human interactome
n_distinct(pickle3.cent$ids) # Interacting proteins
#hs_ppi %>% group_by(InteractorA) %>% summarize(nself = sum(InteractorA==InteractorB)) %>% ungroup() %>% pull(nself) %>% sum
sum(hs_pickle3$InteractorA==hs_pickle3$InteractorB) # Self-interactions

mean(pickle3.cent$cent_deg) # Average interactors (degree)
median(pickle3.cent$cent_deg) # Median interactors (degree)
quantile(pickle3.cent$cent_deg,seq(0,1,len=21)) # quantile degree

#### IntAct ####
#BiocManager::install('vitkl/PItools')
library(PItools)
URL_INTACT="https://ftp.ebi.ac.uk/pub/databases/intact/"
#options(timeout = max(600, getOption("timeout")))
#intact.zip = here::here("data", "intact",  paste0("IntAct-",lastIntActRelease(),".zip"))
#download(paste0(URL_INTACT,"current/psimitab/intact.zip"), destfile = intact.zip)
intact.txt=here::here('data','intact',paste0("IntActRelease_",lastIntActRelease()),"intact.txt")
intact = read_delim(intact.txt)
result = list(data = intact,
             metadata = sprintf("This object contains the data from %scurrent/psimitab/intact.txt and contains all molecular interaction data from the following databases:
                                \"IntAct\", \"MINT\", \"DIP\", \"bhf-ucl\", \"MPIDB\", \"MatrixDB\", \"HPIDb\",\"I2D-IMEx\",\"InnateDB-IMEx\", \"MolCon\", \"UniProt\", \"MBInfo\"",URL_INTACT))
result$data$taxid_A =  str_extract(result$data$`Taxid interactor A`,"(?<=^taxid:)[0-9]+")
result$data$taxid_B =  str_extract(result$data$`Taxid interactor B`,"(?<=taxid:)[0-9]+(?=\\(.+\\)$)")
result$data$within_species = result$data$taxid_A==result$data$taxid_B
class(result) = "RAW_MItab27"

# library(PSICQUIC)
# psicquic <- PSICQUIC()
# providers(psicquic)
# hs_ppi <- interactions(psicquic, id=c("TP53", "MYC"), species="9606")
hs_intact = result$data %>%
         filter(within_species & taxid_A == 9606 & taxid_B == 9606) %>%
         dplyr::rename( id_interactor_A = "#ID(s) interactor A", id_interactor_B = "ID(s) interactor B"  ) %>%
         mutate( id_A = str_remove(id_interactor_A,"^(uniprotkb:|.+:)"),
                 id_B = str_remove(id_interactor_B,'^(uniprotkb:|.+:)')) %>%
         dplyr::select(id_A,id_B) %>% ungroup %>% distinct %>% arrange(id_A)
# intact_ppi = fullInteractome(taxid = 9606,database = 'IntActFTP',
#                              directory =  here::here('data','intact'),
#                              format = "tab27", within_species=T,
#                              clean = TRUE, # parse into usable format (takes 5-10 minutes)
#                              protein_only = F) # filter protein interactions
intact.cent = preload(here::here("prepare","human_intact_ppi_centrality.rds"),
        network.centrality(hs_intact),
        'compute centrality of human interactome from PICKLE3...')

# statistics human interactome
n_distinct(intact.cent$ids) # Interacting proteins
#hs_ppi %>% group_by(InteractorA) %>% summarize(nself = sum(InteractorA==InteractorB)) %>% ungroup() %>% pull(nself) %>% sum
sum(hs_intact$id_A==hs_intact$id_B) # Self-interactions

mean(intact.cent$cent_deg) # Average interactors (degree)
median(intact.cent$cent_deg) # Median interactors (degree)
quantile(intact.cent$cent_deg,seq(0,1,len=21)) # quantile degree

