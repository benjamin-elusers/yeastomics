source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
library(Biostrings)
library(dplyr)
library(furrr)
library(progressr)
library(pbmcapply)
NCPUS=parallel::detectCores()-2 # In my workstation (16-2) = 14 CPUS

# 1. Find the taxon id matching your keywords (case-insensitive) ===============
find.uniprot_refprot(c('human','homo','sapiens','eukaryota'))
# For example, to find the bacteria E. coli
ecoli = find.uniprot_refprot(c('bacteria','Escherichia','coli','K12'))
print(ecoli$tax_id) # Taxon id = 83333

# 2. Get the reference proteome sequences for the selected taxon ===============
ec_aa = get.uniprot.proteome(taxid = 83333, DNA = F)

# 3. Get the reference uniprot accession for the selected taxon ================
ec_uniref = names(ec_aa)

# 4. Get all disorder predictions (from taxon id or from list of ids) ==========
# ------> With uniprot accession <------
# EC_MOBIDB = fetch.mobidb(ec_uniref) # Took 130 sec. for 4402 identifiers from e. coli
# n_distinct(EC_MOBIDB$acc)

# ------> With ncbi taxon id <------
ec_mobidb = load.mobidb(83333) %>% # Took 20 sec. for the protein identifiers of the taxon
            dplyr::filter(acc %in% ec_uniref) # Only use uniprot in reference proteome
n_distinct(ec_mobidb$acc)

# Preparing
ACC = ec_mobidb$acc
START = ec_mobidb$S %>% as.integer()
END = ec_mobidb$E %>% as.integer()

# 5. Get the sequence for disorder stretches ===================================
## Might be slow on single-cpu depending on no. of (features+proteins) to process
## Filtering out unnecessary features/proteins would make it faster


tic('retrieve regions sequence...') # Time the task of retrieving feature sequences
#message('get sequences of mobidb features/regions...')
irows=1:nrow(ec_mobidb)
# pbmclapply to loop across rows of mobidb features in parallel on (n-2) cpus
diso_seq <- pbmcapply::pbmclapply(
                X=irows,
                # extract subsequence of mobidb feature from uniprot protein using positions
                FUN = function(x){
                        feature_seq = subseq(x=ec_aa[[ACC[x]]], start=START[x],end=END[x])
                        return(as.character(feature_seq))
                      },
                mc.cores=NCPUS)
# })
toc(log=T)
### TOOK ~20 SECONDS with pbmcapply() (parallel on (n-2) cpus)

# 6. Compute residue propensities on sequence of mobidb features ==================
# !! LONG COMPUTATION (>5mn on (n-2) cpus) !!

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
# Propsenties used are:
# - hydrophobicity (wimleywhite, kytedoolittle, roseman, camsol)
# - amyloidogenicity (foldamyloid, pawar, aggrescan)
# - interaction propensity (stickiness, voronoi_sickiness)

ec_mobidb_scores = ec_mobidb %>%
                    bind_cols( bind_rows(diso_score) ) %>%
                    group_by(acc,evidence,feature,source,length) %>%
                    summarize( feat_count= unique(content_count),
                               feat_frac = unique(content_fraction),
                               across(aggrescan:wimleywhite, ~ sum(.x)/length)) %>%
                    distinct()


# 8. FILTER CONSENSUS PREDICTION OR INTERESTING FEATURES =======================

# Check description of features from MobiDB at: https://mobidb.bio.unipd.it/about/mobidb
