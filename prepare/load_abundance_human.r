source(here::here("src","__setup_yeastomics__.r"))

### LOAD PAXDB DATASETS
hs.ppm = get.paxdb(tax = 9606, abundance=c('integrated','median'))
# write_rds(hs.ppm,"data/paxdb-human-integrated.rds")

### MAPPING OF ENSEMBL TO UNIPROT IDENTIFIERS
hs.pax2uni = get.hs.ens2uni(transcript = F, ids.ens = unique(hs.ppm$protid)) %>%
   mutate(uniprot = str_extract(UNIPROTID,UNIPROT.nomenclature())) %>% # keep only the UNIPROT Accession without the version number
   dplyr::select(-UNIPROTID) %>% distinct() # keep unique Uniprot identifiers
# write_rds(hs.pax2uni,"data/ensembl-human-uniprot.rds")

### RETRIEVAL OF UNIPROT FEATURES (for all human proteins or only proteins in paxDB)
# !!! TAKES VERY LONG !!! MORE THAN 190K protein identifiers
# uni.feat = load.uniprot.features(tax=9606)
# write_rds(uni.feat,"data/uniprot-human-features.rds")

# pax.feat = uni.feat %>% filter(UNIPROTKB %in% hs.pax2uni$uniprot)
# write_rds(pax.feat,"data/uniprot-human_paxdb-features.rds")

### MERGE UNIPROT/ENSEMBL/PAXDB DATA
# selected.features = c("UNIPROTKB","REVIEWED","EXISTENCE","SCORE","FAMILIES","PATHWAY","PNAME","L")
hs.ppm.uni = left_join(hs.ppm, hs.pax2uni, by=c('protid'='PROTEINID')) #%>%
             #left_join(pax.feat %>% dplyr::select(all_of(selected.features)), by=c("uniprot"='UNIPROTKB'))
# write_rds(hs.ppm.uni,"data/paxdb_integrated-human_uniprot_features.rds")

