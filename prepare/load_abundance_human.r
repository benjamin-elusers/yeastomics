source(here::here("src","__setup_yeastomics__.r"))
### LOAD PAXDB DATASETS
hs_uniref=get.uniprot.proteome(taxid=9606) %>% names
hs_abundance = load.paxdb(tax=9606)
#hs_integrated = get.paxdb(tax = 9606, abundance=c('integrated'))

# CHECK FOR TOTAL SUM OF PPM IN INTEGRATED DATASET (should be close to 1e6 ppm)
# hs_abundance %>%
#   dplyr::filter(is_integrated) %>%
#   group_by(organ) %>%
#   summarize( total_ppm = sum(ppm) ) %>%
#   arrange( desc(total_ppm) ) %>%
#   ggplot(.) + geom_col(aes(y=organ,x=total_ppm),orientation='y') +
#   scale_x_log10()

hs_ppm = summarise_paxdb_abundance(9606)

#hs.ppm = get.paxdb(tax = 9606, abundance=c('integrated','median','mean'))
# write_rds(hs.ppm,"data/paxdb-human-integrated.rds")
### RETRIEVAL OF UNIPROT FEATURES (for all human proteins or only proteins in paxDB)
# !!! TAKES VERY LONG !!! MORE THAN 190K protein identifiers
# uni.feat = load.uniprot.features(tax=9606)
# write_rds(uni.feat,"data/uniprot-human-features.rds")

# pax.feat = uni.feat %>% filter(UNIPROTKB %in% hs.pax2uni$uniprot)
# write_rds(pax.feat,"data/uniprot-human_paxdb-features.rds")

### MERGE UNIPROT/ENSEMBL/PAXDB DATA
# selected.features = c("UNIPROTKB","REVIEWED","EXISTENCE","SCORE","FAMILIES","PATHWAY","PNAME","L")
#hs.ppm.uni = left_join(hs.ppm, hs.pax2uni, by=c('protid'='PROTEINID')) #%>%
             #left_join(pax.feat %>% dplyr::select(all_of(selected.features)), by=c("uniprot"='UNIPROTKB'))
# write_rds(hs.ppm.uni,"data/paxdb_integrated-human_uniprot_features.rds")

