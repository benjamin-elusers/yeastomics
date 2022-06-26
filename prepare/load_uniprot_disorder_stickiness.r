source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
library(Biostrings)



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
# ------> Using a vector of uniprot ids <------
# EC_MOBIDB = fetch.mobidb(ec_uniref) # Took 130 sec. for 4402 identifiers from e. coli
# n_distinct(EC_MOBIDB$acc)

# Using a taxon id
ec_mobidb = load.mobidb(83333) %>% # Took 20 sec. for the protein identifiers of the taxon
            dplyr::filter(acc %in% ec_uniref) # Only use uniprot reference proteome
n_distinct(ec_mobidb$acc)

# 5. Get the sequence for disorder stretches and compute amino-acid scores =====
ec_sticky = ec_mobidb %>%
            rowwise()


library(furrr)
library(progressr)
library(dplyr)

tic('retrieve regions sequence...')
plan(multisession, workers = 14)
with_progress({
  p <- progressor(steps = nrow(ec_sticky), message = 'extract sequence from region start-end positions...')

  diso_seq <- furrr::future_map(1:nrow(ec_sticky), ~{
    p()
    subseq(x=ec_aa[[ec_sticky$acc[.x]]], start=as.integer(ec_sticky$S[.x]), end=as.integer(ec_sticky$E[.x]))
  })
})
toc()

tic('compute aa score...')
with_progress({
  p <- progressor(steps = length(diso_seq), message = 'calculate amino acid scores for sequences...')

  diso_score <- furrr::future_map(seq_along(diso_seq), ~{
    p()
    get_aa_score(diso_seq[[.x]])
  })
})

