# LOAD DATASETS ################################################################
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source(here::here("prepare","setup_idr_analysis.r"),echo = F,chdir = T)
library(Biostrings)
library(hablar)

# BUILD HUMAN IDR DATASET ######################################################
# Filter consensus disordered regions (at least 50% agreement of predictors)
hs_diso =  dplyr::filter(hs_mobidb, is_uniref & feature == 'disorder' & source %in% 'th_50')

cluster <- new_cluster(14)
cluster_library(cluster, "dplyr")
hs_idr = add_idr_sequence(hs_diso, cl=cluster)
hs_peptides = get_peptstats(hs_idr,cl=cluster)
rm(cluster)

# Compute amino acid counts in human IDRs
hs_aacount = get_aa_count(hs_idr,'feature_seq')
# Compute molecular features based on human IDR amino acid counts
hs_aafreq = get_aa_freq(hs_aacount)
hs_charge  = get_aa_charge(hs_aacount)
hs_scores  = get_aa_scores(hs_aacount)
hs_topfreq = get_aa_topfreq(hs_aacount)
hs_foldchange = get_aa_foldchange(hs_aacount)
hs_topfc = get_aa_topfc(hs_foldchange)

# Combine molecular features of human IDRs
hs_features = bind_cols(hs_aacount, hs_aafreq, hs_charge,
                        hs_peptides, hs_foldchange,hs_topfreq, hs_topfc)

# Combine human phase-separation regions
PS_raw = bind_rows(phasepro_wide,phasep_wide)

# Phase-separation proteins with uniprot infos
hs_PS = full_join(PS_raw,UNI_PS, by=c('acc'='AC'))
n0_ps = n_distinct(hs_PS$PS_id)
n0_psprot = n_distinct(hs_PS$acc)
.info$log(sprintf("Phase-separating regions: %s (in %s proteins)",n0_ps,n0_psprot))

# Merging the regions that phase-separates
PS_merged = merge_ps_region(hs_PS)
n1_ps = n_distinct(PS_merged$PS_id)
n1_psprot = n_distinct(PS_merged$acc)
.info$log(sprintf("Phase-separating MERGED regions: %s (in %s proteins)",n1_ps,n1_psprot))


# Dataset for human IDRs with molecular features + phase-separating regions
hs_idr_ps = left_join(hs_idr,hs_features,by=c("IDR_id")) |>
            ungroup() |>
            left_join(PS_merged, by = 'acc') |>
            rowwise() |>
            mutate(has_ps = has_overlap(START,END,PS_START,PS_END))
                   #ps_overlap = which_overlap(START,END,PS_START,PS_END))
# Subset human IDR based on length (at least 35 residues)
# Rename and reorganize the columns
hs_diso = inner_join(hs_uni,hs_diso,by=c('AC'='acc')) %>%
  left_join( HS_AA_feat, by=c('IDR_id','feature_len'='IDR_len') ) %>%
  filter(feature_len>35) %>%
  dplyr::select(-DB,-id_cdna,-OX,-length) %>%
  dplyr::rename(IDR_len=feature_len,IDR_seq=feature_seq,
                IDR_frac=content_fraction,IDR_count=content_count) %>%
  relocate(OS,OS,AC,ID,GN,ensp,NAME,PE,SV,
           PROT_len,IDR_frac,IDR_count,
           IDR_id,IDR_len,S,E,source,feature,IDR_seq,
           positive,negative,netcharge,
           starts_with('pep_'),starts_with("fr_"),starts_with("fc_"), aggrescan:wimleywhite) %>%
  ungroup() %>% distinct()

summary(hs_diso)
write_tsv(hs_diso,file=here::here('prepare','HUMAN_MOBIDB_FEATURES.tsv'))

# BUILD ATAR IDR DATASET #######################################################
# LAST UPDATE OF INPUT DATA *FEB 2023*
IDR_ATAR = rio::import(here::here('prepare',"ATAR_candidates_table.xlsx"))
df_atar = left_join(IDR_ATAR, df_hs_seq, by=c("UNIPROT"='uniprot_id')) |>
          dplyr::select(PROTEIN,UNIPROT,uniprot_seq,atar_sequence=SEQUENCE,
                        START="Start Position",END="End Position", LEN = SEQ_LENGTH) |>
          mutate( IDR_id  = paste0(UNIPROT,"_",START,"..",END))

cluster <- new_cluster(2)
atar_idr = add_idr_sequence(df_atar,cluster)
atar_peptides = get_peptstats(atar_idr,cl=cluster)
rm(cluster)

# compute amino acid counts in idr
atar_aacount = get_aa_count(atar_idr,'atar_sequence')
# compute molecular features based on amino acid counts
atar_aafreq = get_aa_freq(atar_aacount,col_aa = get.AAA())
atar_charge  = get_aa_charge(atar_aacount)
atar_scores  = get_aa_scores(atar_aacount)
atar_topfreq = get_aa_topfreq(atar_aacount)
atar_foldchange = get_aa_foldchange(atar_aacount,col_aa = get.AAA())
atar_topfc = get_aa_topfc(atar_foldchange,col_aa = paste0("fc_",get.AAA()))
#relocate(acc:IDR_id, feature_len, sum_aa)


# All combined molecular features of IDRs
atar_features = bind_cols(atar_aacount, atar_aafreq, atar_charge,
                        atar_peptides, atar_foldchange, atar_topfreq, atar_topfc)

# ALL DISORDER FEATURES
#idr_hiconf=c('th_50')#c("mobidb_lite",'priority','th_50','merge','disprot')
atar_diso = left_join(atar_idr,atar_features,by=c("IDR_id"="pep_IDR_id")) %>% ungroup() %>%
  group_by(UNIPROT,LEN) %>%
  left_join(PS_merged, by = c('UNIPROT'='acc')) |>
  distinct() |>
  mutate(from_atar=TRUE,
         has_ps = has_overlap(START,END,PS_S,PS_E),
         ps_overlap = which_overlap(START,END,PS_S,PS_E),
         has_PS_features = !is.na(PS_overlap)) |>
  distinct()


mobidb_scores = MOBIDB_LONGDISO %>% dplyr::select(-is_longest_feature) %>%
  bind_cols( bind_rows(diso_score) ) %>%
  left_join( AA_feat, by='IDR_id' ) %>% distinct() %>%
  group_by(acc,feature,source,is_longest,feature_len,S,E,
           has_PS_features,) %>%
  distinct() %>%
  mutate( across(aggrescan:wimleywhite, ~ sum(.x)/feature_len)) %>%
  dplyr::select(-c(ORG,feature_len,length,ID,reviewed,PNAME,GNAME,starts_with('PS_ol.'))) %>%
  arrange(PROTEIN,acc,S,E) %>%
 %>%
  relocate(PROTEIN,acc,PROT_len,IDR_frac,IDR_count,has_PS_features,
           IDR_id,S,E,IDR_len,source, feature, is_longest,
           atar_overlap, S_atar,E_atar,
           PS_overlap,PS_id,PS_db,PS_full,PS_len,PS_n,PS_S,PS_E,
           AA3,all_of(aa.class),PLUS,MINUS,NETCHARGE,
           ends_with("_fr"),ends_with("_fc"),
           starts_with('pep'), aggrescan:wimleywhite) %>%
  ungroup() %>% distinct()

write_tsv(mobidb_scores, file = here::here('prepare','ATAR-IDR-FEATURES.tsv'))

#save.image(file = here::here('prepare', 'IDR-features-data.rdata'))
# IDR UMAP #####################################################################

# IDR FEATURES TO USE
# AA Frequencies + by chemical group
# AA Foldchanges
# AA scores (stickiness, roseman aggrescan)
# Peptide stats (average MW, Netcharge, PI, IDR_frac)

umap_data =

col_to_log10 =  c('PROT_len','PEP_len','PEP_mw','IDR_len','IDR_count')

df_atar = mobidb_scores %>% dplyr::select(-starts_with('PS'),-region) %>%
  mutate( across(.cols =col_to_log10,.fns = log10,.names = "{.col}_log10")) %>%
  mutate( across(c(AA3,aa.class,IDR_frac), ~ .x*100) ) %>%
  dplyr::select(-col_to_log10) %>%
  distinct() %>%
  dplyr::relocate(PROTEIN,acc,IDR_frac,has_PS_features,IDR_id,S,E,S_atar,E_atar)

HS_DISO = hs_diso %>% filter( !duplicated(IDR_id) )
hs_mobi_data = HS_DISO %>% ungroup() %>%
  dplyr::select(-c(OS,ID,PE,SV,NAME,region,ensp,source,feature,IDR_seq, starts_with('PS_ol.'))) %>%
  mutate( across(.cols =col_to_log10,.fns = log10,.names = "{.col}_log10")) %>%
  mutate( across(c(AA3,aa.class,IDR_frac), ~ .x*100) ) %>%
  dplyr::select(-col_to_log10) %>%
  distinct() %>% relocate(IDR_id)

### Select features for umap
# all numeric
hs_mobi_num = hs_mobi_data %>% dplyr::select(where(~ is.numeric(.x))) %>%
  dplyr::select(-c(S,E,PS_S,PS_E,PS_len,PS_n))
# Check correlogram of numeric features to remove redundancy
# (high absoluter correlation == redundant features)
ggcorrplot::ggcorrplot(cor(hs_mobi_num),outline.color = 'transparent', tl.cex = 8, tl.srt = 70)

col_AA_fc = paste0(AA3,"_fc")
# non-redundant features
hs_mobi_selected = hs_mobi_data %>%
  dplyr::select(
    all_of(AA3), #aa.class,
    all_of(col_AA_fc),
    voronoi_sickiness, roseman, aggrescan,
    PEP_netcharge, PEP_PI, PEP_len_log10) # PEP_avg_mw  # IDR_frac

# select the same features for atar's IDR
atar_num     = df_atar %>% dplyr::select(colnames(hs_mobi_selected))

# Get non-features columns (protein id, IDR boundaries...)
hs_mobi_info = hs_mobi_data %>% dplyr::select( -colnames(hs_mobi_selected) )
atar_info = df_atar %>% dplyr::select(-colnames(atar_num))

# features for umap distance
df_num = bind_rows(hs_mobi_selected,atar_num)
colnames(df_num)
# non-features for umap identification
df_info = bind_rows(hs_mobi_info,atar_info)


library(umap)
set.seed(142)
umap.config = umap.defaults
umap.config$n_neighbors = 15
# Compute umap based on scaled features
mobi_map <- df_num %>% t() %>% scale() %>% t() %>% umap::umap(seed = 142, config=umap.config)

# mark outliers on umap coordinates (in 1/99% percentiles or outside the interval defined below for X1/X2)
mobi_umap = mobi_map$layout %>% set_colnames(c('X1','X2')) %>%
  bind_cols(mobi_map$data %>% as_tibble() %>% dplyr::rename_with(.fn=xxS, s=".", sx='scaled')) %>%
  bind_cols(df_num) %>%
  bind_cols(df_info) %>%
  mutate(atar_proteins = AC %in% acc[atar_overlap] & !is.na(atar_overlap),
         outliers_x1 = !between(percent_rank(X1),0.01,0.99),
         outliers_x2 = !between(percent_rank(X2),0.01,0.99),
         outliers = !between(X1,-8,8) | !between(X2,-6,6))
summary(mobi_umap[,c('X1','X2')])

#library(ggtrace)
library(ggalt)
library(ggiraph)
library(ggforce)
umap_data = mobi_umap %>% filter(!outliers) %>%
  mutate(idr_lab=paste0(PROTEIN," ",S,"-",E))
atar_idr = subset(umap_data,atar_proteins & atar_overlap)
atar_neighbor = subset(umap_data,atar_proteins & !atar_overlap)

summary(umap_data[,c('X1','X2')])
colnames(df_num)
# Plot the umap
UMAP = ggplot(data=umap_data ,aes(x = X1,y = X2)) +
  # all idrs
  geom_point(size=1.5,shape=16,color='gray',alpha=1) +
  # circle PS idr
  #geom_point(data=subset(umap_data,PS_overlap ), size=3, shape=21, color='black',stroke=0.5) +

  # highlight atar idr
  geom_point(data=atar_idr, aes(color=PROTEIN), shape=16, size=4,alpha=0.9) +
  # geom_text_repel(data=atar_idr,aes(label=idr_lab,color=PROTEIN),
  #                 size=4,fontface='bold', force=5,  max.overlaps=15,force_pull = -1,
  #                 seed=291222) +

  # highlight idr in atar's protein (not overlapping)
  #geom_point(data=atar_neighbor, aes(color=PROTEIN), fill='white',shape=21, stroke=1, size=4, alpha=0.7) +
  #geom_text_repel(data=atar_neighbor,aes(label=idr_lab,color=PROTEIN),
  #               size=3,fontface='italic', force=5, max.overlaps=15, force_pull = -1,
  #              seed=291222) +

  # graphical parameters
labs(x = "UMAP1", y = "UMAP2", subtitle="")+
  scale_color_metro(palette = 'rainbow', discrete = T) +
  theme(legend.position="bottom") + #theme_blackboard()+
  ggeasy::easy_text_size(c("axis.title","axis.text.x", "axis.text.y"), size = 20)+
  ggeasy::easy_remove_legend() +
  ggeasy::easy_remove_x_axis('text') + ggeasy::easy_remove_y_axis('text') +
  theme(aspect.ratio = 1)

# Make interactive points
# ggiraph::geom_point_interactive(size=0.5, color='white',alpha=0, mapping=aes(tooltip=IDR_id)) +
plot(UMAP)
#ggiraph::girafe( code = print(UMAP) )

# save umap in PNG/PDF
ggsave(UMAP, filename = here::here('prepare','umap-atar-idr-human-labelled.png'), height=12,width=12, bg = 'black')
ggsave(UMAP, filename = here::here('prepare','umap-atar-idr-human-labelled.pdf'), height=12,width=12, bg = 'black')

#subset(umap_data,atar_proteins) %>% arrange(IDR_id) %>%
#  dplyr::relocate(PROTEIN,acc,IDR_id,S,E,S_atar,E_atar,atar_proteins,atar_overlap,PS_overlap) %>% View()

