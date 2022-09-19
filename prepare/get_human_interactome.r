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



### homo-oligomers
library(RMySQL)
#To get the refset:
mydb = dbConnect(MySQL(), user='elevy', password='Mysql1!', dbname='3dcomplexV0', host='els141g.weizmann.ac.il')
full_refset = dbSendQuery(mydb, paste0("select * from proteome where refset = '1'"))
data_full_refset = fetch(full_refset, n=-1)
#Code above will give you all the refset proteins from the four model organisms.
count_uniref= data_full_refset %>% ungroup() %>% group_by(org) %>% summarize(nprot=n_distinct(code_sub))
#The column with symmetry info is sym_grom_best.
# please take note that if you want to have all the homomers (and not only the refset) you need to do the following query:
prot_homo = dbSendQuery(mydb, paste0("select * from proteome where refset='1' or
                                    (PDBsym_homo is not null and PDBsym_homo != 'NPS')"))
data_prot_homo = fetch(prot_homo, n=-1) %>% as_tibble() %>% type_convert()


homomers = data_prot_homo %>%
           dplyr::select(OX,OS,org,org_code,
                         uniprot = code_sub, gene, GN, AF=code, pdb=PDBcode,
                         symAF, sym_geom_best, qs_nsub, L=length_full, dom_arch_pfam, ab,abQ,abMax) %>%
          mutate( homomers = factor(qs_nsub,ordered = T, labels = unique(na.omit(map_chr(sort(data_prot_homo$qs_nsub),to_oligomer))) ),
                  multimer = fct_other(homomers,keep=c('monomer','dimer'),other_level = 'multimer'),
                  multi = fct_other(sym_geom_best,keep=c('NPS','C2'),other_level = 'multimer') )

count_homomers = homomers %>%
                  group_by(org) %>%
                  count(multi) %>%
                  add_tally(n,name = 'total_org') %>%
                  ungroup() %>%
                  mutate( fr = 100*n/total_org) %>%
                  mutate( fr = 100*n/total_org)


ggplot(count_homomers, aes(x=multi,y=fr,color=org,fill=org)) +
  geom_col(position = position_dodge(0.8),width = 0.5) +
  geom_text(aes(y=0,label=sprintf("%s",n)),color='white',angle=90, hjust=-0.1,position=position_dodge(width = 0.8),size=3.5) +
  scale_fill_metro() + scale_color_metro()

tab_count_homomers = homomers %>%
  janitor::tabyl(multi,org) %>%
  janitor::adorn_percentages('col') %>%
  janitor::adorn_pct_formatting(digits = 2)

#  janitor::adorn_totals()

#This is because some known homomers are not detected with Alphafold

#You can also access a code that generate the distribution of symmetries per organism + the fraction of homomers in the proteomes of the four organisms there:
#  /media/elusers/users/hugo/15_alphafold/06_proteome_analysis/plot_symmetries_refset.R



##
#/home/benjamin/Desktop/GitHub/yeastomics/data/cavidb
cavspace = here::here('data','cavityspace','AF-all_cavities.idx') %>% read_delim('\t') %>% janitor::clean_names()
bindingdb = here::here('data','bindingdb','BindingDB_All.tsv') %>% read_delim('\t') %>% janitor::clean_names()
dnaprodb =  here::here('data','dnaprodb','collection.txt') %>% rjson::fromJSON(file=.) %>% janitor::clean_names()

dim(cavspace)
head(cavspace)

head(bindingdb)

enumerate_seq = function(BS,k=5){
  w=nchar(BS)
  words = lapply(1:(w-k), function(p0){ subseq(BS, start = p0, width = k) %>% as.character }) %>% unlist
  return(words)
}

enumerate = function(proteome, K=3:12){
  library(pbmcapply)
  N = length(proteome)
  words_by_size=list()
  for( k in K ){
    message(sprintf("Enumerate %s-mer in %s sequences...",k,N))
    words_by_size[[paste0("k",k)]] = pbmclapply(sc_aa, enumerate_seq,k=k,mc.cores=14)
  }
  return(words_by_size)
}

sc_aa=load.sgd.proteome()
sc_kmer = enumerate(sc_aa, K=3:15)

k3 = sc_kmer$k3 %>% unlist %>% table
k4 = sc_kmer$k4 %>% unlist %>% table
length(k4) / min(20^unique(nchar(sc_kmer$k4%>%unlist)), length(sc_kmer$k4 %>% unlist))

nmer  = sapply(sc_kmer, function(x){ unlist(x) %>% length })
kspace = sapply(sc_kmer, function(x){ mers = unlist(x); 20^(unique(nchar(mers))) })

get_repeats = function(words_by_seq,nmin=3){ map_dfrsapply(words_by_seq, function(w){ table(w) %>% as_tibble }) }

sc_words = map_dfr(sc_kmer,stack) %>%
            set_names('word','orf') %>%
            filter( !str_detect(word,"\\*") ) %>%
            mutate(k=nchar(word)) %>%
            group_by(word,orf,k) %>%
            summarize( nr = n()) %>%
            group_by(word,k) %>%
            mutate( NR = sum(nr), NP = n() )

dim(sc_words)
sc_repeats = sc_words %>% filter(nr>1 & NP > 100)
head(sc_repeats)
