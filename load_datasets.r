set.seed(01042021)
script.to.env = function(src,nm){ # load script object into separate environment
  if( !(nm %in% search()) ){
    env = attach(what = NULL,name = nm);
  }else{
    env = as.environment(nm)
  }
  sys.source(file=src,env)
}

# LOAD FUNCTIONS
library(AnnotationDbi)
library(org.Sc.sgd.db)
library(GO.db)
library(Biostrings)
library(tidyverse)
library(stringr)
library(tictoc)
library(openxlsx)

script.to.env("src/utils.r",'utils')
script.to.env("src/function_sequence.r",'seq')
script.to.env("src/function_annotation.r",'annot')
script.to.env("src/function_datalocal.r",'localdata')
script.to.env("src/function_datapub.r",'remotedata')

setup_backup = function(backup_dir ){
  ddmmyyyy=timestamp(stamp = format(Sys.time(), "%d_%m_%Y"),prefix = "", suffix = "",quiet = T)
  next.backup = normalizePath(file.path(backup_dir,ddmmyyyy))
  dir.create(next.backup,showWarnings = F)
  if(!dir.exists(next.backup)){ stop("Failed to create a new backup directory!") }

  message("Data may be backed up at: ",next.backup)
  return(next.backup)
}


# REFERENCE DATA
where = file.path(".downloaded/")
#last.backup = list.dirs(backup)
setup_backup(where)

# REFERENCE DATA URL
## MAIN PROTEIN FEATURES
SGD  = load.sgd.features()
UNIPROT = load.uniprot.features()
#saveRDS(UNIPROT, "data/uniprot-features.rds")

# s. cerevisiae ohnologs (wgd)
OHNO = load.read.csv("http://ygob.ucd.ie/browser/ohnologs.csv",stringsAsFactors = F)[, -1]


# PRIMARY SEQUENCES
S288C  = load.sgd.proteome()
UNI.SC = load.uniprot.proteome('yeast')
#UNI.HS = load.uniprot.proteome('human')
# QUITE LONG! 2hours for 7000 ids (6600 retrieved)
# d2p2.4932 = load.d2p2(ids = names(S288C),saved="data/d2p2-yeast-orf.rds")
# d2p2.4932 = load.d2p2(ids = names(UNI.SC),saved="data/d2p2-yeast-uniprotKB.rds")
# VERY LONG! 5hours for 20000 ids (18000 retrieved)
# d2p2.9606 = load.d2p2(ids = names(UNI.HS),saveto="data/d2p2-human-uniprotKB.rds")

## LOAD DATASETS ---------------------------------------------------------------
#==============================================================================#
# LOCAL DATA
EDL = load.emmanuel.data()                     # Yeast data from emmanuel levy
TAI = load.codon.usage(inputseq=load.sgd.CDS())
#
WAP = load.wapinsky2007.data()       # Evolutionary rate from fungi lineage
K11.R4S = load.rate4site_1011.data() # Evolutionary rate from 1011 strains
#K11 = load.1011.strains()           # 1011 strains proteomes sequences
head(K11.R4S)

R4S = load.aligned.data(data.path="data/" ) # Aligned evolutionary rate

# REMOTE DATA
DUB  = load.dubreuil2019.data(1)        # Stickiness disorder
LEU  = load.leunberger2017.data()       # Limited proteolysis (stability)
VAN  = load.vanleeuwen2020.data()       # Gene dispensability (essentiality)
HO   = load.ho2018.data()               # Unified protein abundance (MS,TAP,GFP)
VIL  = load.villen2017.data()           # Protein Turnover (half-life)
BEL  = load.belle2006.data()            # Protein Half-lives (NOT CORRELATED WITH THE REST)
GEI  = load.geisberg2014.data(nodesc=T) # mRNA half-lives

# load.byrne2005.data()                 # ohnologs are on different columns (~550 rows)
WGD  = get.sc.ohno()                    # Whole-genome duplication data (ohnologs are a single column ~ 1000 rows)
JAC  = load.peter2018.data()          # 1011 strains project details


#saveRDS(SGD,file='/media/elusers/users/benjamin/A-PROJECTS/02_Scripts/R/forMeta/sgd-features.rds')
#saveRDS(DUB,file='/media/elusers/users/benjamin/A-PROJECTS/02_Scripts/R/forMeta/dubreuil-2019.rds')
#saveRDS(UNIPROT, "/media/elusers/users/benjamin/A-PROJECTS/02_Scripts/R/forMeta/uniprot-features.rds")
#saveRDS(VAN,"/media/elusers/users/benjamin/A-PROJECTS/02_Scripts/R/forMeta/essential-dispensable.rds")

GO   = get.uniprot.go(names(UNI.SC))
LOC  = get.uniprot.localization(annot = UNIPROT,loc_to_columns = T)
PMID = get.uniprot.pmid(names(UNI.SC))
UNISGD  = get.uniprot.sgd(names(UNI.SC))

head(PMID)
head(LOC)
head(GO)


get.uniprot.sgd()
load.uniprot.features(tax = 559292)
