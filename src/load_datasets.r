set.seed(09122020)
## LOAD DATASETS ---------------------------------------------------------------
#==============================================================================#
# LOCAL DATA
source("function_dataloc.r")
EDL = load.emmanuel.data()          # Yeast data from emmanuel levy
R4S = load.aligned.data()           # Aligned evolutionary rate
# REMOTE DATA
source("function_datapub.r")
LEU  = load.leunenberger2017.data()     # Limited proteolysis (stability)
VAN  = load.vanleeuwen2020.data()       # Gene dispensability (essentiality)
HO   = load.ho2018.data()               # Unified protein abundance (MS,TAP,GFP)
VIL  = load.villen2017.data()           # Protein Turnover (half-life)
BEL  = load.belle2006.data()            # Protein Half-lives
GEI  = load.geisberg2014.data(nodesc=T) # mRNA half-lives
WGD  = load.byrne2005.data()            # Whole-genome duplication data (ohnologs)
