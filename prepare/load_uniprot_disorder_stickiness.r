library(tidyverse)
library(here)
# requires Yeastomics from my github
# can be loaded directly from url
url_yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/"
source(paste0(url_yeastomics,"__setup_yeastomics__.r"))

library(Biostrings)
library(dplyr)
library(furrr)
library(progressr)
library(pbmcapply)
NCPUS=parallel::detectCores()-2 # on my workstation (16-2) = 14 CPUS

UNIPROT_IDR_ATAR = c( "H0WFA5", "Q7T226", "O00444", "Q15468", "Q6UVJ0", "O95613",
                      "Q14676", "P43351", "Q12888", "P23771", "P03372", "P45973",
                      "P15992", "P31539", "Q12329", "P27695", "Q03690", "P19097",
                      "P05453", "Q06787", "Q04637", "Q9U2F5", "P35637", "Q9NQI0",
                      "Q7Z3E1", "P54253", "Q12988", "Q16082", "P09651", "P17133",
                      "P02545", "P20591", "Q8WUM0", "Q02630", "P14907", "P52948",
                      "O43791", "P10636", "Q9UER7", "O60885", "P78545", "Q13618",
                      "O60563", "Q14781", "Q01130", "O43670", "Q15648", "Q9NWB1",
                      "P06748", "Q92547", "Q9NY12", "P22087", "Q13148", "P34689",
                      "G5EBV6", "Q95XR4", "Q9N303")


PHASESEP_DB_LLPS = readxl::read_excel("/home/benjamin/Downloads/Compressed/phasepdb/phasepdbv2_llps.xlsx")
#PHASESEP_DB_MLO  = readxl::read_excel("/home/benjamin/Downloads/Compressed/phasepdb/phasepdbv2_mlolt_mloht.xlsx")
#colnames(PHASESEP_DB_LLPS)


# 1. Find the taxon id matching your keywords (case-insensitive) ===============
# For example, to find the bacteria E. coli
#ecoli = find.uniprot_refprot(c('bacteria','Escherichia','coli','K12'))
#print(ecoli$tax_id) # Taxon id = 83333

# 2. Get the reference proteome sequences for the selected taxon ===============
# ec_aa = get.uniprot.proteome(taxid = 83333, DNA = F)

IDR_prot = Biostrings::readAAStringSet("/home/benjamin/Desktop/GitHub/yeastomics/prepare/uniprot_idr.fasta")
names(IDR_prot) = str_extract(names(IDR_prot),UNIPROT.nomenclature())

# 3. Get the reference uniprot accession for the selected taxon ================
# ec_uniref = names(ec_aa)

# 4. Get all disorder predictions (from taxon id or from list of ids) ==========
# Took 130 sec. for 4402 identifiers from e. coli
# ------> With uniprot accession <------
MOBIDB = fetch.mobidb(UNIPROT_IDR_ATAR)
        #left_join(PHASESEP_DB_LLPS, by=c('acc'='uniprot_entry'))
PHASESEP_DB_LLPS[ PHASESEP_DB_LLPS$uniprot_entry %in% UNIPROT_IDR_ATAR,]

# PHASE SEPARATION FEATURES
phasepro = MOBIDB %>% filter(feature == 'phase_separation') %>%  mutate(PS_acc_pos = paste0(acc,"_",S,"-",E))
phasepro_wide = phasepro %>%
                pivot_wider(id_cols=c('acc','PS_acc_pos'), names_from='source', values_from = c('S','E'))

MOBIDB_PS = MOBIDB %>% ungroup() %>%
  filter( feature != 'phase_separation') %>%
  left_join(phasepro_wide, by = 'acc') %>%
  mutate(PS_overlap_cter = (S_phasepro <= E & S < S_phasepro & E_phasepro > E),
         PS_overlap_nter = (S_phasepro <= S & (E_phasepro > S & E_phasepro < E)),
         PS_overlap_nested = (S <= S_phasepro & E_phasepro <= E),
         PS_overlap_inside = (S >= S_phasepro & E <= E_phasepro),
         PS_overlap = PS_overlap_cter | PS_overlap_nter | PS_overlap_inside | PS_overlap_nested
  ) %>%
  group_by(acc) %>% mutate(has_PS_features = !is.na(PS_overlap))


# ALL DISORDER FEATURES
idr_hiconf=c("mobidb_lite",'priority','th_50','merge','disprot')
MOBIDB_DISO = MOBIDB_PS %>%
              filter( feature == 'disorder' & source %in% idr_hiconf) %>%
              group_by(acc,length,source,content_fraction,content_count) %>%
              mutate(is_longest_feature = feature_len == max_(feature_len)) %>%
              group_by(acc,length) %>% mutate(is_longest = feature_len == max_(feature_len))

# LONGEST DISORDER FEATURES
MOBIDB_LONGDISO = MOBIDB_DISO %>%
        filter( (is_longest & feature_len > 35) | PS_overlap ) %>%
        arrange(acc,has_PS_features,evidence,feature,source,S,E,is_longest_feature,is_longest) %>%
        relocate(acc,has_PS_features,evidence,feature,source,S,E,feature_len,is_longest_feature,is_longest)
table(MOBIDB_LONGDISO$acc)


# ------> With ncbi taxon id <------
# ec_mobidb = load.mobidb(83333) %>% # Took 20 sec. for the protein identifiers of the taxon
#             dplyr::filter(acc %in% ec_uniref) %>% # Only use uniprot in reference proteome
#             group_by(acc) %>%
#             mutate(feature_maxlen = max_(content_count),
#                    longest_idr = feature == 'disorder' & content_count == feature_maxlen )
# n_distinct(ec_mobidb$acc)

# 5. Get the sequence for disorder stretches ===================================
## Might be slow on single-cpu depending on no. of (features+proteins) to process
## Filtering out unnecessary features/proteins would make it faster

tic('retrieve regions sequence...') # Time the task of retrieving feature sequences
#message('get sequences of mobidb features/regions...')
irows=1:nrow(MOBIDB_LONGDISO)

# Preparing data for accessing feature sequences
ACC = MOBIDB_LONGDISO$acc
START = MOBIDB_LONGDISO$S %>% as.integer()
END = MOBIDB_LONGDISO$E %>% as.integer()

# pbmclapply to loop across rows of mobidb features in parallel on (n-2) cpus
diso_seq <- pbmcapply::pbmclapply(
                X=irows,
                # extract subsequence of mobidb feature from uniprot protein using positions
                FUN = function(x){
                        feature_seq = subseq(x=IDR_prot[[ACC[x]]], start=START[x],end=END[x]) %>% as.character
                }, mc.cores = NCPUS)

MOBIDB_LONGDISO$feature_seq = unlist(diso_seq)
toc()
### TOOK ~20 SECONDS with pbmcapply() (parallel on (n-2) cpus)

# 6. Compute residue propensities on sequence of mobidb features ==================
# !! LONG COMPUTATION (>5mn on (n-2) cpus) !!

# Residue propensities used are the same as in Fig3. of Dubreuil et al. 2019 (JMB):
# - hydrophobicity (wimleywhite, kytedoolittle, roseman, camsol)
# - amyloidogenicity (foldamyloid, pawar, aggrescan)
# - interaction propensity (stickiness, voronoi_sickiness)

tic('compute aa scores of mobidb features ...')
#message('get amino acid scores of mobidb features/regions...')
iseq = 1:length(diso_seq)
diso_score <- pbmcapply::pbmclapply(
                  X=iseq,
                  # Using subsequence of mobidb features, compute sum of residue propensities
                  FUN = function(x){ get_aa_score(string = diso_seq[[x]]) },
                  mc.cores = NCPUS)
toc(log=T)
### TOOK ~340 SECONDS with pbmcapply() (parallel on (n-2) cpus)

# 7. Get average residue propensity for all proteins ===========================
# Average propensity = sum aa score / feature length
names(diso_seq) = iseq
# AA count
aacount = lapply(diso_seq, function(x){ alphabetFrequency(Biostrings::AAString(x),as.prob = F)[AA_STANDARD] }) %>%
          bind_rows %>%
          dplyr::rename(setNames(get.AA1(),get.AA3()))

aacount %>%
  rowwise() %>%
  mutate( AA1 = colnames(.)[which.max(c_across(1:20))],
          AA2 = colnames(.)[which.max(c_across(-AA1))],
          AA3 = colnames(.)[which.max(c_across(-c(AA1,AA2)))],
          AA4 = colnames(.)[which.max(c_across(-c(AA1,AA2,AA3)))]
          )


mobidb_scores = MOBIDB_LONGDISO %>% dplyr::select(-S_merge,-E_merge,-is_longest_feature) %>%
                  bind_cols( aacount ) %>%
                  bind_cols( bind_rows(diso_score) ) %>%
                  group_by(acc,length,evidence,feature,source,is_longest,feature_len,S,E,
                             has_PS_features,PS_overlap,S_phasepro,E_phasepro,
                             PS_overlap_cter,PS_overlap_nter,PS_overlap_nested,PS_overlap_inside) %>%
                  distinct() %>%
                  mutate( across(aggrescan:wimleywhite, ~ sum(.x)/feature_len))

# 8. YOUR OWN FILTER ===========================================================
# ... You can add your own code here to filter the data ...

# EXAMPLE 1: Selecting longest disordered region in each protein
mobidb_scores %>% filter( is_longest_idr )
# EXAMPLE 2: get the mobidb sub-regions enriched for various characteristics such as:
# polyampholyte
# positive/negative electrolyte
# G/C/P rich
# polar
# ...

table(mobidb_scores$source)
mobidb_scores %>% filter(source == 'mobidb_lite_sub')
mobidb_scores %>% filter(source == 'phasepro')

# Check description of features from MobiDB at: https://mobidb.bio.unipd.it/about/mobidb


# RULE1 : TAKE ALL FEATURES THAT OVERLAP PHASESEPDB or PHASEPRO

mobidb_scores %>% filter(is_longest_idr)
RULE1 = mobidb_scores %>%
  ungroup() %>%
  filter( source != 'phasepro') %>%
  left_join(phasepro_wide, by = 'acc') %>%
  mutate(overlap_phasepro_cter = (S_phasepro <= E & S < S_phasepro),
         overlap_phasepro_nter = (S_phasepro <= S & E_phasepro < E),
         overlap_phasepro_nested = (S <= S_phasepro & E_phasepro <= E),
         overlap_phasepro_inside = (S >= S_phasepro & E <= E_phasepro),
         overlap_phasepro = overlap_phasepro_cter | overlap_phasepro_nter | overlap_phasepro_inside | overlap_phasepro_nested
  ) %>% filter(overlap_phasepro)

n_distinct(RULE1$acc)

# RULE2 : FILTER FEATURE BELOW 35 RESIDUES
RULE2 = RULE1 %>% filter( feature_len > 35 )
n_distinct(RULE2$acc)
# RULE3 : CHECK IF FEATURE OVERLAP INTERACTION FROM ring
#ring = mobidb_scores %>% filter(source == 'ring') %>%
#       mutate(acc_pos = paste0(acc,"_",S,"-",E))


pivot_wider(id_cols=c('acc','acc_pos'), names_from='source', values_from = c('S','E'))
RULE3 = RULE2 %>% group_by(acc) %>% mutate(

# RULE4 : LONGEST FEATURE OR MORE ABUNDANT STRETCHES
# RULE5 : CHECK IF IN TOP/BOTTOM 5/10% OF AVERAGE STICKINESS
# RULE6 : CHECK THE ENRICHMENT FOR MobiDB-lite sub regions




  # 1. longest region, boundaries
  # 2. protein size
  # 2. region size,
  # 3-6. Four most frequent AA and frequency
  # 7-10. Four most over-represented AA and fold increase relative to mean
  # 11. number of + charges
  # 12. number of - charges
  # 13. number of his
  # 14. stickiness
  # 15. MOBI info (phase separation, zzz)
  # 16. sequence
