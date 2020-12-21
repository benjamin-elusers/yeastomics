set.seed(09122020)
script.to.env = function(src,nm){ # load script object into separate environment
  if( !(nm %in% search()) ){
    env = attach(what = NULL,name = nm);
  }else{
    env = as.environment(nm)
  }
  sys.source(file=src,env)
}
d2p2.orf=readRDS("data/d2p2-yeast-orf.rds")
d2p2.uni=readRDS("data/d2p2-yeast-uniprotKB.rds")

# REFERENCE DATA
library(AnnotationDbi)
library(org.Sc.sgd.db)
library(GO.db)
library(Biostrings)
library(tidyverse)
library(stringr)
library(tictoc)
script.to.env("src/utils.r",'utils')
script.to.env("src/function_sequence.r",'seq')

SGD  = load.sgd.features()
S288C  = load.sgd.proteome()

## LOAD DATASETS ---------------------------------------------------------------
#==============================================================================#
# LOCAL DATA
script.to.env("src/function_datalocal.r",'localdata')
EDL = load.emmanuel.data()                     # Yeast data from emmanuel levy
TAI = load.codon.usage(inputseq=load.sgd.CDS())
#DP2 = load.d2p2(ids = names(S288C)) # QUITE LONG ~ 2hours for 7000ids -> load as RDS object
WAP = load.wapinsky2007.data()
R4S = load.aligned.data(data.path="data/" )    # Aligned evolutionary rate

# REMOTE DATA
script.to.env("src/function_datapub.r",'remotedata')
DUB  = load.dubreuil2019.data()         # Stickiness disorder
LEU  = load.leunberger2017.data()       # Limited proteolysis (stability)
VAN  = load.vanleeuwen2020.data()       # Gene dispensability (essentiality)
HO   = load.ho2018.data()               # Unified protein abundance (MS,TAP,GFP)
VIL  = load.villen2017.data()           # Protein Turnover (half-life)
BEL  = load.belle2006.data()            # Protein Half-lives (NOT CORRELATED WITH THE REST)
GEI  = load.geisberg2014.data(nodesc=T) # mRNA half-lives
# load.byrne2005.data()                 # ohnologs are on different columns (~550 rows)
WGD  = get.sc.ohno()                    # Whole-genome duplication data (ohnologs are a single column ~ 1000 rows)
JAC  = load.jackson2018.data()          # 1011 strains project details

G1011 = load.proteome