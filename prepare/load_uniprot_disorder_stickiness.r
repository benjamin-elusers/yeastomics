source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
library(Biostrings)
library(furrr)
library(progressr)
library(pbmcapply)
#handlers(global = TRUE)
#handlers("progress")
library(dplyr)
plan(multisession, workers = availableCores()-2)

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
            dplyr::filter(acc %in% ec_uniref) %>%  # Only use uniprot reference proteome
            rowwise()
n_distinct(ec_mobidb$acc)

# 5. Get the sequence for disorder stretches and compute amino-acid scores =====
# !!!!! This is longest step (might help to filter disorder regions first) !!!!!

### PARALLEL
tic('retrieve regions sequence...')
# with_progress({
#   p <- progressor(steps = nrow(ec_sticky))
message('get sequences of mobidb features/regions...')
diso_seq <- pbmcapply::pbmclapply(X=1:nrow(ec_mobidb),
                                  mc.cores = parallel::detectCores()-2,
                                  FUN = function(x){
                                        ACC = ec_mobidb$acc[x]
                                        START = ec_mobidb$S[x] %>% as.integer()
                                        END = ec_mobidb$E[x] %>% as.integer()
                                        #p(message = sprintf('get regions from %s [%s/%s]',ACC,x,nrow(ec_mobidb)))
                                        feature_seq = subseq(x=ec_aa[[ACC]], start=START,end=END)
                                        return(as.character(feature_seq))
                                  })
# })
toc(log=T)
### TOOK 15 SECONDS with pbmcapply() (parallel on 14 cpus)

tic('compute aa scores of mobidb features ...')
diso_score <- pbmcapply::pbmclapply(X=1:length(diso_seq),
                                  FUN = function(x){
                                    SEQ = diso_seq[[x]]
                                    cat(SEQ)
                                    get_aa_score(string = SEQ)
                                  },mc.cores = parallel::detectCores()-2)
toc(log=T)
### TOOK  SECONDS with pbmcapply() (parallel on 14 cpus)

# with_progress({
#   p <- progressor(steps = length(diso_seq))
#
#   diso_score <- furrr::future_map(seq_along(diso_seq), ~{
#     p(message = 'get aa scores...')
#   })
# })
