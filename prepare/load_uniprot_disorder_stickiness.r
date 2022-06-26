source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
library(Biostrings)

fetch.mobidb = function(ids,to.df=T){
  # Check mobidb data description at the following url:
  # https://mobidb.bio.unipd.it/about/mobidb
  library(pbmcapply)
  library(parallel)
  library(tictoc)
  cpus=detectCores()-2
  tic('get disorder from MobiDB')
  message(sprintf('retrieve disorder predictions from MobiDB for %s identifiers (uniprot) ...',n_distinct(ids)))
  mobidb = pbmcapply::pbmclapply(ids, get.mobidb.id,mc.cores=cpus)
  names(mobidb) = ids
  toc()

  if(to.df){
    return( purrr::compact(mobidb) %>% bind_rows() %>% )
  }

  return(mobidb)
}

load.mobidb = function(taxon){
  # Check mobidb data description at the following url:
  # https://mobidb.bio.unipd.it/about/mobidb
  tic('get mobidb for taxon from url....')
  mobidb = readr::read_delim(sprintf("https://mobidb.bio.unipd.it/api/download?ncbi_taxon_id=%s&format=tsv",taxon)) %>%
           separate(col=feature, sep='-',into=c('evidence','feature','source')) %>%
           separate_rows('start..end',sep=',') %>%
           separate('start..end',into = c('S','E'),convert = T) %>%
           mutate( feature_len = E-S+1 )

  toc()
  return(mobidb)
}


# 1. Find the taxon id matching your keywords (case-insensitive)
find.uniprot_refprot(c('human','homo','sapiens','eukaryota'))
# For example, to find the bacteria E. coli
ecoli = find.uniprot_refprot(c('bacteria','Escherichia','coli','K12'))
print(ecoli$tax_id) # Taxon id = 83333

# 2. Get the reference proteome sequences for the selected taxon
ec_aa = get.uniprot.proteome(taxid = 83333, DNA = F)

# 3. Get the reference uniprot accession for the selected taxon
ec_uniref = names(ecoli_prot)

# 4. Get all disorder predictions (from taxon id or from list of ids)
# using a vector of uniprot ids
EC_MOBIDB = fetch.mobidb(ec_uniref) # Took 130 sec. for 4402 identifiers from e. coli
# using a taxon id
ec_mobidb = load.mobidb(83333) # much faster and return potentially more results

n_distinct(ec_diso$acc)
# 5. Get the sequence for disorder stretches and compute amino-acid scores


ec_sticky =
            rowwise() %>%

                    feature_aa  = as.character(subseq(ec_aa[[acc]],start=as.integer(S), end=as.integer(E))) )

head(ec_sticky)
table( ec_sticky$feature_len )

Biostrings::extract_character_from_XString_by_ranges(ec_aa[[1]],start = 10L,width=10L)
