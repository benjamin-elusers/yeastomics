source(here::here("src","__setup_yeastomics__.r"))

### LOAD PAXDB DATASETS
hs_abundance = load.paxdb(tax=9606)
hs_ppm_int = get.paxdb(tax = 9606, abundance=c('integrated'))

hs_ppm = hs_ppm_int %>%
  group_by(taxid,protid) %>%
  summarise(
          ppm_whole_org = ppm_int[ organ == 'WHOLE_ORGANISM' ],
          ppm_max = max_(ppm_int),
          ppm_min = min_(ppm_int),
          ppm_sd = sd_(ppm_int),
          ppm_md = median_(ppm_int),
          n_organ = n_distinct(organ),
          #rk_organ = rank(-(ppm_int),ties.method = 'average'),
          max_organ = paste0(organ[ppm_int == ppm_max],collapse='/'),
          min_organ = paste0(organ[ppm_int == ppm_min],collapse='/')) %>%
          distinct() %>%
  group_by(taxid) %>%
  mutate( pc_whole_org = 100*percent_rank(-ppm_whole_org),
          pc_max = 100*percent_rank(-ppm_max),
          pc_min = 100*percent_rank(-ppm_min),
          pc_md = 100*percent_rank(-ppm_md),
          q_whole_org =  cume_dist(ppm_whole_org),
          q_max = cume_dist(ppm_max),
          q_min = cume_dist(ppm_min),
          q_md =  cume_dist(ppm_md)
        )


### MAPPING OF ENSEMBL TO UNIPROT IDENTIFIERS
hs_pax2uni = get.hs.ens2uni(transcript = F, ids.ens = unique(hs_ppm$protid)) %>%
  mutate(uniprot = str_extract(UNIPROTID,UNIPROT.nomenclature())) %>% # keep only the UNIPROT Accession without the version number
  dplyr::select(-UNIPROTID) %>% distinct() # keep unique Uniprot identifiers
# write_rds(hs.pax2uni,"data/ensembl-human-uniprot.rds")

hs_uni_feat = read_rds( here::here("data","uniprot-human-features.rds") )
hs_ppm_uni = left_join(hs_ppm, hs_pax2uni, by=c('protid'='PROTEINID')) %>%

write_rds(hs_ppm_uni,"data/paxdb_integrated_human.rds")

read_rds("data/paxdb_integrated_human.rds")
#hs.ppm = get.paxdb(tax = 9606, abundance=c('median','mean'))

test = hs_abundance %>%
  group_by(organ,protid) %>%
  mutate(ppm_max_organ = max(ppm), rk_organ = dense_rank(ppm)) %>%
  group_by(protid) %>%
  mutate(nd = n_distinct(id), breadth = n_distinct(organ), ppm_max = max(ppm), rk = dense_rank(ppm)) %>%
  arrange(protid,organ)

#%>%
#  group_by(is_integrated==T,protid) %>% mutate(ppm_int = mean_(ppm), rk_int= dense_rank(ppm))

hs.ppm = get.paxdb(tax = 9606, abundance=c('integrated','median','mean'))
# write_rds(hs.ppm,"data/paxdb-human-integrated.rds")


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

