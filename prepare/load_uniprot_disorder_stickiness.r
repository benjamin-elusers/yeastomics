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

## AMINO ACID CLASSES
aa.prop = seqinr::SEQINR.UTIL$AA.PROPERTY %>%
  append( list('Alcohol'=c('S','T'),
               'Turnlike'=c('A','C','D','E','G','H','K','N','Q','R','S','T'))
  )
sum.aa.fr = function(BS,alphabet,to100=F){
  alphabet.freq = rowSums(letterFrequency(BS,as.prob = F,letters = alphabet))
  prot.enrich = tibble(id=names(BS),fr=alphabet.freq)
  if(to100){ prot.enrich$fr = 100*prot.enrich$fr }
  return(prot.enrich)
}
aa.set = map_chr(aa.prop,str_c,collapse='')
aa.class = paste0(tolower(sub(x=names(aa.set),"\\.","")),"_",aa.set)


# UNIPROT_IDR_ATAR = c( "H0WFA5", "Q7T226", "O00444", "Q15468", "Q6UVJ0", "O95613",
#                       "Q14676", "P43351", "Q12888", "P23771", "P03372", "P45973",
#                       "P15992", "P31539", "Q12329", "P27695", "Q03690", "P19097",
#                       "P05453", "Q06787", "Q04637", "Q9U2F5", "P35637", "Q9NQI0",
#                       "Q7Z3E1", "P54253", "Q12988", "Q16082", "P09651", "P17133",
#                       "P02545", "P20591", "Q8WUM0", "Q02630", "P14907", "P52948",
#                       "O43791", "P10636", "Q9UER7", "O60885", "P78545", "Q13618",
#                       "O60563", "Q14781", "Q01130", "O43670", "Q15648", "Q9NWB1",
#                       "P06748", "Q92547", "Q9NY12", "P22087", "Q13148", "P34689",
#                       "G5EBV6", "Q95XR4", "Q9N303")
IDR_ATAR = rio::import(here::here('prepare',"IDR_Candidate_new.xlsx"))

#PHASESEP_DB_LLPS = readxl::read_excel("/home/benjamin/Downloads/Compressed/phasepdb/phasepdbv2_llps.xlsx")
#PHASESEP_DB_MLO  = readxl::read_excel("/home/benjamin/Downloads/Compressed/phasepdb/phasepdbv2_mlolt_mloht.xlsx")
#colnames(PHASESEP_DB_LLPS)


# 1. Find the taxon id matching your keywords (case-insensitive) ===============
# For example, to find the bacteria E. coli
#ecoli = find.uniprot_refprot(c('bacteria','Escherichia','coli','K12'))
#print(ecoli$tax_id) # Taxon id = 83333

human = find.uniprot_refprot(c('9606','HUMAN','homo sapiens'))  %>%
        arrange(desc(keyword_matched))
# 2. Get the reference proteome sequences for the selected taxon ===============
# ec_aa = get.uniprot.proteome(taxid = 83333, DNA = F)
hs_aa = get.uniprot.proteome(taxid = human$tax_id, DNA = F)

IDR_prot = Biostrings::readAAStringSet(here::here("prepare","IDR_Candidate_new.fasta"))
#IDR_prot = Biostrings::readAAStringSet("/home/benjamin/Desktop/GitHub/yeastomics/prepare/uniprot_idr.fasta")
names(IDR_prot) = str_extract(names(IDR_prot),UNIPROT.nomenclature())

# 3. Get the reference uniprot accession for the selected taxon ================
# ec_uniref = names(ec_aa)
hs_uniref = names(hs_aa)

# 4. Get all disorder predictions (from taxon id or from list of ids) ==========
# Took 130 sec. for 4402 identifiers from e. coli
# ------> With uniprot accession <------
MOBIDB = fetch.mobidb(IDR_ATAR$UNIPROT)
UNI = pbmcapply::pbmclapply(IDR_ATAR$UNIPROT,get_uniprot_id,mc.cores = 14) %>%
      purrr::compact() %>%
      bind_rows() %>%
      dplyr::rename(AC=Entry,ID='Entry Name',reviewed=Reviewed,
                    PNAME='Protein names', GNAME = 'Gene Names',
                    ORG = Organism, prot_len = Length)

MOBIDB = left_join(MOBIDB,UNI,by=c('acc'='AC'))
#left_join(PHASESEP_DB_LLPS, by=c('acc'='uniprot_entry'))
#PHASESEP_DB_LLPS[ PHASESEP_DB_LLPS$uniprot_entry %in% IDR_ATAR$UNIPROT,]

# PHASE SEPARATION FEATURES
phasepro = MOBIDB %>% filter(feature == 'phase_separation') %>%  mutate(PS_acc_pos = paste0(acc,"_",S,"..",E))
phasepro_wide = phasepro %>% filter(source=='phasepro') %>%
                pivot_wider(id_cols=c('acc','PS_acc_pos'), names_from='source', values_from = c('S','E'), values_fn=unique)

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
        filter( (is_longest & feature_len > 35) | (PS_overlap & feature_len > 35)) %>%
        mutate(ID_region = paste0(acc,"_",S,"..",E)) %>%
        arrange(acc,has_PS_features,evidence,ID_region,feature,source,S,E,is_longest_feature,is_longest) %>%
        relocate(acc,has_PS_features,evidence,ID_region,feature,source,S,E,feature_len,is_longest_feature,is_longest)

#table(MOBIDB_LONGDISO$acc)

# ------> With ncbi taxon id <------
# ec_mobidb = load.mobidb(83333) %>% # Took 20 sec. for the protein identifiers of the taxon
#             dplyr::filter(acc %in% ec_uniref) %>% # Only use uniprot in reference proteome
#             group_by(acc) %>%
#             mutate(feature_maxlen = max_(content_count),
#                    longest_idr = feature == 'disorder' & content_count == feature_maxlen )
# n_distinct(ec_mobidb$acc)


# 5. Get all human disorder predictions ========================================
hs_mobidb = preload(here::here('data','mobidb-human-features.rds'),
                    load.mobidb(human$tax_id),
                    'retrieve human mobidb features...')

df_hs_diso = hs_mobidb %>%
  dplyr::filter(acc %in% hs_uniref & feature == 'disorder' & source %in% 'th_50') %>%
  mutate(ID_region = paste0(acc,"_",S,"..",E))

irows=1:nrow(df_hs_diso)
hs_diso_seq_list <- pbmcapply::pbmclapply(
  X=irows,
  # extract subsequence of mobidb feature from uniprot protein using positions
  FUN = function(x){
    feature_seq = tryCatch(
      subseq(x=hs_aa[[df_hs_diso$acc[x]]], start=df_hs_diso$S[x],end=df_hs_diso$E[x]) %>% as.character,
      error = function(e){ finally = print(paste0(df_hs_diso$acc[x],df_hs_diso$S[x],df_hs_diso$E[x])); return(NA) })
  }, mc.cores = NCPUS)

df_hs_diso$feature_seq = unlist(hs_diso_seq_list)

hs_diso_seq = AAStringSet(df_hs_diso$feature_seq[!is.na(df_hs_diso$feature_seq)])
names(hs_diso_seq) = unique(df_hs_diso$ID_region[!is.na(df_hs_diso$feature_seq)])

###### HS AA DISO FEATURES #######
AA1=get.AA1() %>% as.vector()
AA3=get.AA3() %>% as.vector()

HS_AA_COUNT = (letterFrequency(hs_diso_seq,as.prob = F,letters = get.AA1())) %>%
  bind_cols( ID_region=names(hs_diso_seq), region_len = widths(hs_diso_seq)) %>%
  dplyr::rename(setNames(AA1,AA3))
HS_AA_FR = HS_AA_COUNT %>% mutate( across(AA3, ~ . / region_len) )

hs_mobidb_aa = left_join(df_hs_diso , HS_AA_COUNT, by='ID_region')
tot_aa_hs=sum(widths(hs_aa[unique(df_hs_diso$acc)]))
HS_DISO_AACOUNT = hs_mobidb_aa  %>% group_by(source) %>%
  summarize( across(.cols=AA3, .fns = sum_ ),
             total_count = sum_(c_across(AA3)),
             total_freq =total_count/ tot_aa_hs)
HS_DISO_AAFREQ = (HS_DISO_AACOUNT[,AA3] / HS_DISO_AACOUNT$total_count) %>% as.double()

# COUNT CLASS OF AA
HS_AACLASS_COUNT = map(aa.prop, sum.aa.fr, BS=hs_diso_seq) %>% bind_cols %>%
  rename_with(~aa.class, starts_with('fr')) %>%
  dplyr::select(-starts_with('id')) %>%
  mutate(ID_region=names(hs_diso_seq), region_len = widths(hs_diso_seq))

HS_AACLASS_FR = HS_AACLASS_COUNT %>% mutate( across(aa.class, ~ . / region_len) )

HS_TOP4 = HS_AA_FR %>% mutate( across(AA3, ~ . *100)) %>%
  pivot_longer(cols = 1:20) %>%
  group_by(ID_region,region_len) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,2))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=c(ID_region,region_len), names_from='rk', names_glue = 'AA{rk}_fr',
              values_from = c('name.val'))

HS_AA_FC = sweep( HS_AA_FR[,AA3],2,HS_DISO_AAFREQ,"/")
HS_TOP4_FC = HS_AA_FC %>%
  bind_cols(ID_region= HS_AA_FR$ID_region, region_len =HS_AA_FR$region_len) %>%
  pivot_longer(cols = 1:20) %>%
  group_by(ID_region) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,2))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=c(ID_region,region_len), names_from='rk', names_glue = 'AA{rk}_fc',
              values_from = c('name.val'))

HS_AA_CHARGE = HS_AA_FR %>% group_by(ID_region,region_len) %>%
  summarize(PLUS=sum(LYS+ARG+HIS),MINUS=sum(ASP+GLU),
            NETCHARGE=PLUS-MINUS,
            pep_len = Peptides::lengthpep(hs_diso_seq),
            pep_mw  = Peptides::mw(hs_diso_seq),
            pep_avg_mw = pep_mw / pep_len,
            pep_netcharge = Peptides::charge(hs_diso_seq),
            pep_PI = Peptides::pI(hs_diso_seq))

HS_AA_feat = left_join(HS_AA_FR,HS_AACLASS_FR) %>%
  left_join(HS_AA_CHARGE,by=c('ID_region','region_len')) %>%
  left_join(HS_TOP4,by=c('ID_region','region_len')) %>%
  left_join(HS_TOP4_FC,by=c('ID_region','region_len'))


tic('compute aa scores of mobidb features ...')
#message('get amino acid scores of mobidb features/regions...')
iseq = 1:length(hs_diso_seq_list)
hs_diso_score <- pbmcapply::pbmclapply(
  X=iseq,
  # Using subsequence of mobidb features, compute sum of residue propensities
  FUN = function(x){ get_aa_score(string = hs_diso_seq_list[[x]]) },
  mc.cores = NCPUS)
toc(log=T)
### TOOK ~340 SECONDS with pbmcapply() (parallel on (n-2) cpus)

hs_uni = get_uniprot_reference(9606)

hs_mobidb_scores = df_hs_diso %>%
  left_join( HS_AA_feat, by=c('ID_region','feature_len'='region_len') ) %>%
  bind_cols( bind_rows(hs_diso_score) ) %>%
  group_by(acc,feature,source,feature_len,S,E) %>%
  distinct() %>%
  mutate( across(aggrescan:wimleywhite, ~ sum(.x)/feature_len)) %>%
  dplyr::select(-evidence) %>%
  relocate(acc, ID_region,S,E,feature_len, source, feature_seq,
           get.AA3() %>% as.vector,ends_with("_fr")
  ) %>% filter(feature_len > 35)

hs_diso$atar_selection = hs_diso$ID_region %in% MOBIDB_LONGDISO$ID_region

hs_diso = inner_join(hs_uni,hs_mobidb_scores,by=c('AC'='acc')) %>%
           #left_join(get.width(hs_aa),by=c('AC'='orf')) %>%
  filter(feature_len>35) %>%
  mutate(atar_selection = ID_region %in% MOBIDB_LONGDISO$ID_region) %>%
  dplyr::select(-DB,-id_cdna,-OX) %>%
  dplyr::rename(prot_len=length,idr_len=feature_len,idr_seq=feature_seq,
                idr_frac=content_fraction,idr_count=content_count) %>%
  relocate(OS,OS,AC,ID,GN,ensp,NAME,PE,SV,
           prot_len,idr_frac,idr_count,
           atar_selection,ID_region,idr_len,S,E,source,feature,idr_seq,
           AA3,ends_with("_fr"),all_of(aa.class),PLUS,MINUS,
           ends_with("_fc"), aggrescan:wimleywhite)


write_tsv(hs_diso,file=here::here('prepare','HUMAN_MOBIDB_FEATURES.tsv'))
#View(hs_diso)

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
diso_seq_list <- pbmcapply::pbmclapply(
                X=irows,
                # extract subsequence of mobidb feature from uniprot protein using positions
                FUN = function(x){
                        feature_seq = subseq(x=IDR_prot[[ACC[x]]], start=START[x],end=END[x]) %>% as.character
                }, mc.cores = NCPUS)

diso_seq = AAStringSet(unique(unlist(diso_seq_list)))
names(diso_seq) = unique(paste0(ACC,"_",START,"..",END))

MOBIDB_LONGDISO$feature_seq = unlist(diso_seq_list)
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
iseq = 1:length(diso_seq_list)
diso_score <- pbmcapply::pbmclapply(
                  X=iseq,
                  # Using subsequence of mobidb features, compute sum of residue propensities
                  FUN = function(x){ get_aa_score(string = diso_seq_list[[x]]) },
                  mc.cores = NCPUS)
toc(log=T)
### TOOK ~340 SECONDS with pbmcapply() (parallel on (n-2) cpus)

# 7. Get average residue propensity for all proteins ===========================
# Average propensity = sum aa score / feature length
names(diso_seq_list) = iseq

###### ATAR DISO AA FEATURES ######
# COUNT AA
AA.FR = (letterFrequency(diso_seq,as.prob = T,letters = get.AA1()) * 100) %>%
        bind_cols( ID_region=names(diso_seq)) %>%
          dplyr::rename(setNames(get.AA1(),get.AA3()))

hs_aa_fr_diso = matrix(HS_DISO_AAFREQ , byrow=T,ncol=20,nrow=nrow(AA.FR))

IDR_AAFR =  AA.FR[,get.AA3()]
AA.FC = sweep(IDR_AAFR,2,100*HS_DISO_AAFREQ,"/")

TOP4_FC = AA.FC %>% add_column(ID_region= AA.FR$ID_region) %>%
  pivot_longer(cols = 1:20) %>%
  group_by(ID_region) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,1))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=ID_region, names_from='rk', names_glue = 'AA{rk}_fc',
              values_from = c('name.val'))

AACLASS.FR = map(aa.prop, sum.aa.fr, BS=diso_seq) %>% bind_cols %>%
             rename_with(~aa.class, starts_with('fr')) %>%
             dplyr::select(-starts_with('id')) %>%
             mutate(ID_region=names(diso_seq))

TOP4 = AA.FR %>%
  pivot_longer(cols = 1:20) %>%
  group_by(ID_region) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,1))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=ID_region, names_from='rk', names_glue = 'AA{rk}_fr',
              values_from = c('name.val'))

AA.CHARGE = AA.FR %>% group_by(ID_region) %>%
  summarize(PLUS=sum(LYS+ARG+HIS),MINUS=sum(ASP+GLU),
            NETCHARGE=PLUS-MINUS,
            pep_len = Peptides::lengthpep(diso_seq),
            pep_mw  = Peptides::mw(diso_seq),
            pep_avg_mw = pep_mw / pep_len,
            pep_netcharge = Peptides::charge(diso_seq),
            pep_PI = Peptides::pI(diso_seq))

AA_feat = left_join(AA.FR,AACLASS.FR) %>%
  left_join(AA.CHARGE,by="ID_region") %>%
  left_join(TOP4,by='ID_region') %>%
  left_join(TOP4_FC,by='ID_region')


mobidb_scores = MOBIDB_LONGDISO %>% dplyr::select(-is_longest_feature) %>%
                  left_join( AA_feat, by='ID_region' ) %>%
                  bind_cols( bind_rows(diso_score) ) %>%
                  group_by(acc,feature,source,is_longest,feature_len,S,E,
                             has_PS_features,PS_overlap,S_phasepro,E_phasepro,
                             PS_overlap_cter,PS_overlap_nter,PS_overlap_nested,PS_overlap_inside) %>%
                  distinct() %>%
                  mutate( across(aggrescan:wimleywhite, ~ sum(.x)/feature_len)) %>%
                  dplyr::select(-evidence,-length) %>%
                 relocate(ORG,acc,prot_len,ID,reviewed,PNAME,GNAME,has_PS_features,
                          PS_acc_pos,S_phasepro,E_phasepro, starts_with("PS_overlap"),
                          ID_region,S,E,feature_len, is_longest, source, feature_seq,
                          get.AA3() %>% as.vector,ends_with("_fr"),ends_with("_fc")
                          )
colnames(mobidb_scores)


write_tsv(mobidb_scores, file = here::here('prepare','ATAR-IDR-FEATURES.tsv'))

save.image(file = here::here('prepare', 'IDR-features-data.rdata'))

# 8. CLUSTERING OF DISORDERED REGIONS ==========================================
### UMAP + Clustering
load(here::here('prepare', 'IDR-features-data.rdata'))

#hs_diso = hs_mobidb_scores %>% filter(feature_len>35)
#hs_diso$atar_selection = hs_diso$ID_region %in% mobidb_scores$ID_region

hs_mobi_data = hs_diso %>% ungroup() %>% drop_na() %>% dplyr::select(where(is.numeric),-S,-E,-PE,-SV)
mobi_id =  hs_diso %>% ungroup() %>% drop_na() %>% pull(ID_region)
atar_selection = hs_diso %>% ungroup() %>% drop_na() %>% pull(atar_selection)
set.seed(142)
library(umap)

umap.config = umap.defaults
umap.config$n_neighbors = 7
library(M3C)

hs_scaled_mobi = t(hs_mobi_data) %>%  scale() %>% t() %>% as_tibble()
rownames(hs_scaled_mobi) = mobi_id
mobi_map <- hs_scaled_mobi %>% umap::umap(seed = 142, config=umap.config)
mobi_umap = mobi_map$layout %>% set_colnames(c('X1','X2')) %>%
            bind_cols(mobi_map$data,id=mobi_id) %>%
            mutate(uni=str_extract(id,UNIPROT.nomenclature()))

library(ggiraph)
UMAP = ggplot(data=mobi_umap,aes(x = X1,y = X2, tooltip=id)) +
  geom_point(size=1, alpha=0.2, aes(color=NET)) +
  # atar selection
  geom_point(data=subset(mobi_umap,atar_selection), size=2, color='red') +
  geom_text_repel(data=subset(mobi_umap,atar_selection),aes(label=id),max.overlaps = 50,size=3.5,color='red') +
  ggiraph::geom_point_interactive(size=0.5, mapping=aes(tooltip=id)) +
  labs(x = "X1", y = "X2", subtitle="UMAP plot")+
  theme(legend.position="bottom")


plot(UMAP)
ggiraph::girafe( code = print(UMAP) )

#K <- M3C::M3C(mobi_data_t, method=2,maxK = 30, seed = 142)
#ATAR=subset(mobi_map$data,atar_selection)
#LOWX2 = mobi_map$data %>% group_by(AC)  %>% filter(rk_X2 < 20)
ggsave(UMAP, filename = here::here('prepare','umap-atar-idr-human.png'), height=12,width=12, bg = 'white')
