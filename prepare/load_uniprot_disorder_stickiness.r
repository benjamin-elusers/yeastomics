# requires Yeastomics from my github; can be loaded directly from url
url_yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/"
source(paste0(url_yeastomics,"__setup_yeastomics__.r"))
pkg=c('here','tidyverse','dplyr','furrr','progressr','pbmcapply','Biostrings')
xfun::pkg_load2(pkg)

NCPUS=parallel::detectCores()-2 # on my workstation (16-2) = 14 CPUS

## AMINO ACID CLASSES
AA1=get.AA1() %>% as.vector()
AA3=get.AA3() %>% as.vector()

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

# 0. Find taxon id matching keywords ===========================================
human = find.uniprot_refprot(c('9606','HUMAN','homo sapiens'))  %>%
        arrange(desc(keyword_matched))

hs_uni = get_uniprot_reference(human$tax_id)
## Reference proteome sequences ================================================
hs_aa = get.uniprot.proteome(taxid = human$tax_id, DNA = F)
## Reference uniprot accession =================================================
hs_uniref = names(hs_aa)

#### ...................................................................... ####
# A.) HUMAN --------------------------------------------------------------------
# 1. Get all human disorder predictions ========================================
hs_mobidb = preload(here::here('data','mobidb-human-features.rds'),
                    load.mobidb(human$tax_id),
                    'retrieve human mobidb features...')

df_hs_diso = hs_mobidb %>%
  dplyr::filter(acc %in% hs_uniref & feature == 'disorder' & source %in% 'th_50') %>%
  mutate(IDR_id = paste0(acc,"_",S,"..",E))

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

# 2. Propensity scores for human disorder ======================================
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

df_hs_mobidb = df_hs_diso %>%
  bind_cols( bind_rows(hs_diso_score) ) %>%
  group_by(acc,feature,source,feature_len,S,E) %>%
  distinct() %>%
  mutate( across(aggrescan:wimleywhite, ~ sum(.x)/feature_len)) %>%
  dplyr::select(-evidence) %>%
  filter( !is.na(feature_seq) )

hs_diso_seq = AAStringSet(df_hs_mobidb$feature_seq)
names(hs_diso_seq) = unique(df_hs_mobidb$IDR_id)

# 3. Molecular features for human disorder =====================================
HS_AA_COUNT = letterFrequency(hs_diso_seq,as.prob = F,letters = AA1) %>%
  bind_cols( IDR_id=names(hs_diso_seq), IDR_len = widths(hs_diso_seq)) %>%
  dplyr::rename(setNames(AA1,AA3))
HS_AA_FR = HS_AA_COUNT %>% mutate( across(AA3, ~ . / IDR_len) )

hs_mobidb_aa = left_join(df_hs_diso , HS_AA_COUNT, by='IDR_id')
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
  mutate(IDR_id=names(hs_diso_seq), IDR_len = widths(hs_diso_seq))

HS_AACLASS_FR = HS_AACLASS_COUNT %>% mutate( across(aa.class, ~ . / IDR_len) )

HS_TOP4 = HS_AA_FR %>% mutate( across(AA3, ~ . *100)) %>%
  pivot_longer(cols = 1:20) %>%
  group_by(IDR_id,IDR_len) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,2))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=c(IDR_id,IDR_len), names_from='rk', names_glue = 'AA{rk}_fr',
              values_from = c('name.val'))

HS_AA_FC = sweep( HS_AA_FR[,AA3],2,HS_DISO_AAFREQ,"/")
HS_TOP4_FC = HS_AA_FC %>%
  bind_cols(IDR_id= HS_AA_FR$IDR_id, IDR_len =HS_AA_FR$IDR_len) %>%
  pivot_longer(cols = 1:20) %>%
  group_by(IDR_id) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,2))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=c(IDR_id,IDR_len), names_from='rk', names_glue = 'AA{rk}_fc',
              values_from = c('name.val'))

HS_AA_CHARGE = HS_AA_COUNT %>% group_by(IDR_id,IDR_len) %>%
  summarize(PLUS=sum(LYS+ARG+HIS),MINUS=sum(ASP+GLU),
            NETCHARGE=(PLUS-MINUS), PLUS=PLUS/IDR_len, MINUS=MINUS/IDR_len)

HS_PEP = df_hs_diso %>% ungroup %>%
  dplyr::select(IDR_id, IDR_len = feature_len,feature_seq) %>%
  rowwise %>%
  mutate(
    PEP_len = Peptides::lengthpep(feature_seq),
    PEP_mw  = Peptides::mw(feature_seq),
    PEP_avg_mw = PEP_mw / PEP_len,
    PEP_netcharge = Peptides::charge(feature_seq),
    PEP_PI = Peptides::pI(feature_seq)) %>%
  dplyr::select(-feature_seq)

HS_AA_feat = left_join(HS_AA_FR,HS_AACLASS_FR) %>%
  left_join(HS_AA_CHARGE,by=c('IDR_id','IDR_len')) %>%
  left_join(HS_PEP,by=c('IDR_id','IDR_len')) %>%
  left_join(HS_TOP4,by=c('IDR_id','IDR_len')) %>%
  left_join(HS_TOP4_FC,by=c('IDR_id','IDR_len'))

# 4. Phase-separation features for human =======================================
# Phasepro
# phasepro = RJSONIO::fromJSON("https://phasepro.elte.hu/download_full.json") %>%
#   bind_rows %>% type_convert()
phasepro = hs_mobidb %>%
  filter(feature == 'phase_separation') %>%
  mutate(PS_id = paste0(acc,"_",S,"..",E))
phasepro_wide = phasepro %>%
  filter(source=='phasepro') %>%
  pivot_wider(id_cols=c('acc','PS_id'),
              names_from='source',
              values_from = c('S','E'), values_fn=unique) %>%
  dplyr::rename(PS_S=S_phasepro,PS_E=E_phasepro) %>% mutate(PS_db='phasepro')
# Phasepdb
phasesep_db_llps = rio::import("http://db.phasep.pro/static/db/database/phaseodbv2_1_llps.xlsx", na = c("","_")) %>%
  mutate(region = str_replace(region,"–","-") )%>%
  filter(str_detect(organism,'Homo sapiens')) %>%
  separate_rows(region,sep=',') %>%
  separate(region,remove = F, into=c('PS_S','PS_E'),sep='-') %>%
  rowwise() %>%
  mutate( PS_S = ifelse(is_number(PS_S),parse_integer(PS_S),NA),
          PS_E = ifelse(is_number(PS_E),parse_integer(PS_E),NA)) %>%
  type_convert() %>%
  mutate(PS_db='phasepdb',
         PS_id = ifelse( is_number(PS_S) & is_number(PS_E),
                         paste0(uniprot_entry,"_",PS_S,"..",PS_E),
                         paste0(uniprot_entry,"_",str_replace_all(region," ","-")))
  ) %>%
  filter( !is.na(PS_S) & !is.na(PS_E) | !is.na(region))  %>%
  arrange(PS_id)

phasesep_wide = phasesep_db_llps %>%
  dplyr::select(acc=uniprot_entry,PS_db,PS_id,PS_S,PS_E,region)

PS_raw = bind_rows(phasepro_wide,phasesep_wide)

uni_ps = pbmcapply::pbmclapply(PS_raw$acc,get_uniprot_id,mc.cores = 14) %>%
         purrr::compact() %>%
         bind_rows() %>%
         dplyr::select(AC=Entry,PROT_len = Length)

PS_all  = PS_raw %>%
          left_join(uni_ps,by=c('acc'='AC')) %>%
          mutate(PS_len = PS_E-PS_S+1, PS_full = PROT_len == PS_len) %>%
          distinct() %>%
          add_count(acc,name='PS_n') %>%
          ungroup() %>%
          mutate(row=row_number())

library(GenomicRanges)
g0 = PS_all %>% filter(PS_full | PS_n==1 )
g1 = PS_all %>% filter(!(PS_full | PS_n==1 ))
# Remove PS regions with no clear boundaries
#PS_all %>% filter( !(row %in% c(g0$row,g1$row) ) )

g1.1 = g1 %>%
       dplyr::select(db=PS_db,seqnames=acc,start=PS_S,end=PS_E) %>%
       as("GRanges") %>%
       reduce() %>%
       as_tibble() %>%
       dplyr::select(-strand) %>%
       dplyr::rename(acc=seqnames, PS_S=start,PS_E=end,PS_len=width) %>%
       mutate(PS_db ='merged',PS_id = paste0(acc,"_",PS_S,"..",PS_E)) %>%
       mutate(region = paste0(PS_S,"-",PS_E))

PS_merged = bind_rows(g1.1,g0) %>%
  dplyr::select(-c(PROT_len,PS_full,PS_n,row) ) %>%
  left_join(uni_ps,by=c('acc'='AC')) %>%
  mutate(PS_full = PROT_len == PS_len) %>%
  ungroup() %>% distinct() %>%
  add_count(acc,name='PS_n') %>%
  arrange(PS_id) %>%
  dplyr::select(-PROT_len) %>%
  mutate(PS_n = PS_n - ((PS_n>1)*PS_full) )

df_hs_diso = df_hs_mobidb %>% ungroup() %>%
  filter( feature != 'phase_separation') %>%
  left_join(PS_merged, by = 'acc') %>%
  mutate(PS_ol.cter = (PS_S <= E & S < PS_S & PS_E > E),
         PS_ol.nter = (PS_S <= S & (PS_E > S & PS_E < E)),
         PS_ol.nest = (S <= PS_S & PS_E <= E),
         PS_ol.in = (S >= PS_S & E <= PS_E),
         PS_overlap = PS_ol.cter | PS_ol.nter | PS_ol.in | PS_ol.nest
  ) %>%
  group_by(acc) %>% mutate(has_PS_features = !is.na(PS_overlap))

hs_diso = inner_join(hs_uni,df_hs_diso,by=c('AC'='acc')) %>%
  left_join( HS_AA_feat, by=c('IDR_id','feature_len'='IDR_len') ) %>%
  filter(feature_len>35) %>%
  dplyr::select(-DB,-id_cdna,-OX) %>%
  dplyr::rename(PROT_len=length,IDR_len=feature_len,IDR_seq=feature_seq,
                IDR_frac=content_fraction,IDR_count=content_count) %>%
  relocate(OS,OS,AC,ID,GN,ensp,NAME,PE,SV,
           PROT_len,IDR_frac,IDR_count,
           IDR_id,IDR_len,S,E,source,feature,IDR_seq,
           all_of(AA3),ends_with("_fr"),all_of(aa.class),
           PLUS,MINUS,NETCHARGE, starts_with('PEP'),
           ends_with("_fc"), aggrescan:wimleywhite) %>%
  ungroup() %>% distinct()

summary(hs_diso)
write_tsv(hs_diso,file=here::here('prepare','HUMAN_MOBIDB_FEATURES.tsv'))

#### ...................................................................... ####

# B.) ATAR'S PROTEIN  ----------------------------------------------------------
###
# 5. Get disorder for list of uniprot ==========================================
# ------> With uniprot accession <------
#IDR_ATAR = rio::import(here::here('prepare',"IDR_Candidate_new.xlsx"))
IDR_ATAR = rio::import(here::here('prepare',"ATAR_candidates_human.xlsx"))
idr_atar = IDR_ATAR %>%
           dplyr::select(PROTEIN,UNIPROT,S_atar="Start Position",E_atar="End Position")

MOBIDB_ATAR = fetch.mobidb(IDR_ATAR$UNIPROT)
UNI_ATAR = pbmcapply::pbmclapply(IDR_ATAR$UNIPROT,get_uniprot_id,mc.cores = 14) %>%
  purrr::compact() %>%
  bind_rows() %>%
  dplyr::rename(AC=Entry,ID='Entry Name',reviewed=Reviewed,
                PNAME='Protein names', GNAME = 'Gene Names',
                ORG = Organism, PROT_len = Length)

MOBIDB_ATAR = left_join(MOBIDB_ATAR,UNI_ATAR,by=c('acc'='AC'))

IDR_prot = hs_aa[IDR_ATAR$UNIPROT]
# Biostrings::readAAStringSet(here::here("prepare","IDR_Candidate_new.fasta"))
# IDR_prot = Biostrings::readAAStringSet("/home/benjamin/Desktop/GitHub/yeastomics/prepare/uniprot_idr.fasta")

# ALL DISORDER FEATURES
idr_hiconf=c('th_50')#c("mobidb_lite",'priority','th_50','merge','disprot')
MOBIDB_DISO = left_join(MOBIDB_ATAR,PS_merged,by=c('acc')) %>%
              filter( feature == 'disorder' & source %in% idr_hiconf) %>%
              group_by(acc,length,source,content_fraction,content_count) %>%
              mutate(is_longest_feature = feature_len == max_(feature_len)) %>%
              group_by(acc,length) %>% mutate(is_longest = feature_len == max_(feature_len)) %>%
              distinct() %>%
              mutate(PS_ol.cter = (PS_S <= E & S < PS_S & PS_E > E),
                     PS_ol.nter = (PS_S <= S & (PS_E > S & PS_E < E)),
                     PS_ol.nest = (S <= PS_S & PS_E <= E),
                     PS_ol.in = (S >= PS_S & E <= PS_E),
                     PS_overlap = PS_ol.cter | PS_ol.nter | PS_ol.in | PS_ol.nest
              ) %>%
              group_by(acc) %>%
              mutate(has_PS_features = !is.na(PS_overlap)) %>%
              left_join(idr_atar,by=c('acc'='UNIPROT')) %>%
              ungroup() %>% rowwise() %>%
              mutate(atar_ol.cter = (S_atar <= E & S < S_atar & E_atar > E),
                     atar_ol.nter = (S_atar <= S & (E_atar > S & E_atar < E)),
                     atar_ol.nest = (S <= S_atar & E_atar <= E),
                     atar_ol.in = (S >= S_atar & E <= E_atar),
                     atar_overlap = atar_ol.cter | atar_ol.nter | atar_ol.nest | atar_ol.in
              )


# LONGEST DISORDER FEATURES
MOBIDB_LONGDISO = MOBIDB_DISO %>%
        filter( (is_longest & feature_len > 35) | (PS_overlap & feature_len > 35)) %>%
        mutate(IDR_id = paste0(acc,"_",S,"..",E)) %>%
        dplyr::select(-evidence) %>%
        arrange(PROTEIN,PROT_len,acc,S,E,IDR_id,atar_overlap,has_PS_features,PS_overlap) %>%
        relocate(PROTEIN,PROT_len,acc,IDR_id,feature_len,atar_overlap,has_PS_features,PS_overlap,is_longest_feature,is_longest,feature,source,S,E)
#View(MOBIDB_LONGDISO)

# 6. Get the sequence for disorder stretches ===================================
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

# 7. Compute residue propensities on sequence of mobidb features ==================
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

# 8. Get average residue propensity for all proteins ===========================
# Average propensity = sum aa score / feature length
names(diso_seq_list) = iseq

###### ATAR DISO AA FEATURES ######
# COUNT AA
AA.COUNT = letterFrequency(diso_seq,as.prob = F,letters = get.AA1()) %>%
  bind_cols( IDR_id=names(diso_seq), IDR_len = widths(diso_seq)) %>%
  dplyr::rename(setNames(get.AA1(),get.AA3()))

AA.FR = AA.COUNT %>% mutate(across(get.AA3() %>% as.vector, ~ ./IDR_len ))

hs_aa_fr_diso = matrix(HS_DISO_AAFREQ , byrow=T,ncol=20,nrow=nrow(AA.FR))

IDR_AAFR =  AA.FR[,get.AA3()]
AA.FC = sweep(IDR_AAFR,2,HS_DISO_AAFREQ,"/")

TOP4_FC = AA.FC %>% add_column(IDR_id= AA.FR$IDR_id) %>%
  pivot_longer(cols = 1:20) %>%
  group_by(IDR_id) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,1))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=IDR_id, names_from='rk', names_glue = 'AA{rk}_fc',
              values_from = c('name.val'))


AACLASS_COUNT = map(aa.prop, sum.aa.fr, BS=diso_seq) %>% bind_cols %>%
             rename_with(~aa.class, starts_with('fr')) %>%
             dplyr::select(-starts_with('id')) %>%
             mutate(IDR_id=names(diso_seq), IDR_len=width(diso_seq))
AACLASS_FR = AACLASS_COUNT %>% mutate( across(aa.class, ~ . / IDR_len) )

TOP4 = AA.FR %>% mutate(across(AA3,~ .*100)) %>%
  pivot_longer(cols = 1:20) %>%
  group_by(IDR_id) %>%
  slice_max(order_by = value,n = 4,with_ties = F) %>%
  mutate(rk = rank(-value,ties.method = 'first'),
         name.val = paste0(name,":",round(value,2))) %>%
  dplyr::select(-name,-value) %>%
  pivot_wider(id_cols=IDR_id, names_from='rk', names_glue = 'AA{rk}_fr',
              values_from = c('name.val'))


AA.CHARGE = AA.COUNT %>% group_by(IDR_id,IDR_len) %>%
  summarize(PLUS=sum(LYS+ARG+HIS),MINUS=sum(ASP+GLU),
            NETCHARGE=(PLUS-MINUS), PLUS=PLUS/IDR_len, MINUS=MINUS/IDR_len)

AA.PEP = MOBIDB_LONGDISO %>% ungroup() %>%
         dplyr::select(IDR_id, feature_seq) %>%
         rowwise() %>%
         mutate(PEP_len = Peptides::lengthpep(feature_seq),
                PEP_mw  = Peptides::mw(feature_seq),
                PEP_avg_mw = PEP_mw / PEP_len,
                PEP_netcharge = Peptides::charge(feature_seq),
                PEP_PI = Peptides::pI(feature_seq)) %>%
        dplyr::select(-feature_seq)

AA_feat = left_join(AA.FR,AACLASS_FR,by=c('IDR_id','IDR_len')) %>%
  left_join(AA.CHARGE,by=c('IDR_id','IDR_len')) %>%
  left_join(AA.PEP,by=c('IDR_id')) %>%
  left_join(TOP4,by=c('IDR_id')) %>%
  left_join(TOP4_FC,by=c('IDR_id'))


mobidb_scores = MOBIDB_LONGDISO %>% dplyr::select(-is_longest_feature) %>%
                  bind_cols( bind_rows(diso_score) ) %>%
                  left_join( AA_feat, by='IDR_id' ) %>% distinct() %>%
                  group_by(acc,feature,source,is_longest,feature_len,S,E,
                             has_PS_features,) %>%
                  distinct() %>%
                  mutate( across(aggrescan:wimleywhite, ~ sum(.x)/feature_len)) %>%
                  dplyr::select(-c(ORG,feature_len,length,ID,reviewed,PNAME,GNAME,starts_with('PS_ol.'))) %>%
                  arrange(PROTEIN,acc,S,E) %>%
                  dplyr::rename(IDR_seq=feature_seq,
                                IDR_frac=content_fraction,
                                IDR_count=content_count) %>%
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
# 9. CLUSTERING OF DISORDERED REGIONS ==========================================
### UMAP + Clustering
#load(here::here('prepare', 'IDR-features-data.rdata'))
#hs_diso = hs_mobidb_scores %>% filter(feature_len>35)
#hs_diso$atar_selection = hs_diso$IDR_id %in% mobidb_scores$IDR_id
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

hs_mobi_num = hs_mobi_data %>% dplyr::select(where(~ is.numeric(.x))) %>%
              dplyr::select(-c(S,E,PS_S,PS_E,PS_len,PS_n))
atar_num     = df_atar %>% dplyr::select(colnames(hs_mobi_num))

hs_mobi_info = hs_mobi_data %>% dplyr::select( -colnames(hs_mobi_num) )
atar_info = df_atar %>% dplyr::select(-colnames(atar_num)) %>% add_column(is_atar_idr = T)

df_num = bind_rows(hs_mobi_num,atar_num)
colnames(df_num)
df_info = bind_rows(hs_mobi_info,atar_info)

#hs_scaled_mobi = t( ) %>% drop_na) %>% scale() %>% t() %>% as_tibble()
#mobi_id =  hs_mobi_data$IDR_id
#rownames(hs_scaled_mobi) = mobi_id

set.seed(142)
library(umap)

umap.config = umap.defaults
umap.config$n_neighbors = 15
#library(M3C)
#hs_diso %>% filter( IDR_id %in% get.dup(IDR_id) )

mobi_map <- df_num %>% t() %>% scale() %>% t() %>% umap::umap(seed = 142, config=umap.config)

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
library(chameleon)
library(ggforce)
umap_data = mobi_umap %>% filter(!outliers)
summary(umap_data[,c('X1','X2')])



UMAP = ggplot(data=umap_data ,aes(x = X1,y = X2)) +
  labs(x = "X1", y = "X2", subtitle="")+
  #ggiraph::geom_point_interactive(size=0.5, color='white',alpha=0, mapping=aes(tooltip=IDR_id)) +
  geom_point(size=1.5,shape=16,color='white',alpha=0.3) +
  geom_point(data=subset(umap_data,PS_overlap ), size=3, shape=21, color='white',stroke=0.5) +

  # atar selection
  geom_point(data=subset(umap_data,atar_proteins & atar_overlap), aes(color=PROTEIN), shape=16, size=4,alpha=0.9) +
  geom_text_repel(data=subset(umap_data,atar_proteins & atar_overlap),aes(label=paste0(PROTEIN,"\n",S,"-",E),color=PROTEIN),max.overlaps = 50,size=4,fontface='bold') +
  geom_point(data=subset(umap_data,atar_proteins & !atar_overlap), aes(color=PROTEIN), fill='white',shape=21, stroke=1, size=4, alpha=0.7) +
  geom_text_repel(data=subset(umap_data,atar_proteins & !atar_overlap),aes(label=paste0(PROTEIN,"\n",S,"-",E),color=PROTEIN),max.overlaps = 50,size=3,fontface='italic') +


  # geom_point(data=subset(umap_data,is_atar_idr), size=5, alpha=1, color='green',shape=16) +
  # geom_text_repel(data=subset(umap_data,is_atar_idr),
  #                 aes(label=paste0(PROTEIN,"_",S_atar,"-",E_atar)),
  #                 max.overlaps = 50,size=3,color='green',force=20) +
  #geom_point(data=subset(umap_data,atar_overlap), size=2, shape=16, alpha=1, color='white') +
  #geom_text_repel(data=subset(umap_data,atar_overlap & atar_proteins),
  #                aes(label=paste0(PROTEIN,"_",S,"-",E)),
  #                max.overlaps = 50,size=3,color='white',force=20) +
  # geom_text_repel(data=subset(umap_data,atar_proteins & PS_overlap),
  #                 aes(label=PS_id),
  #                 max.overlaps = 50,size=3,color='purple',force=20)
  scale_color_metro(palette = 'rainbow', discrete = T) +
  theme(legend.position="bottom") + theme_blackboard()+
  ggeasy::easy_text_size(c("axis.title","axis.text.x", "axis.text.y"), size = 20)+
  ggeasy::easy_remove_legend()

plot(UMAP)
#ggiraph::girafe( code = print(UMAP) )

#K <- M3C::M3C(mobi_data_t, method=2,maxK = 30, seed = 142)
#ATAR=subset(mobi_map$data,atar_selection)
#LOWX2 = mobi_map$data %>% group_by(AC)  %>% filter(rk_X2 < 20)
ggsave(UMAP, filename = here::here('prepare','umap-atar-idr-human.png'), height=12,width=12, bg = 'black')
ggsave(UMAP, filename = here::here('prepare','umap-atar-idr-human.pdf'), height=12,width=12, bg = 'black')

#subset(umap_data,atar_overlap) %>% arrange(IDR_id) %>% print(n=30)
#subset(umap_data,atar_proteins) %>% arrange(IDR_id) %>%
#  dplyr::relocate(PROTEIN,acc,IDR_id,S,E,S_atar,E_atar,atar_proteins,atar_overlap,PS_overlap) %>% View()

idr_atar
