load(here::here('output','hs_datasets.rdata'))
load(here::here('output','hs_integrated_datasets.rdata'))
load(here::here('output','hs_predictors_proteome.rdata'))

source(here::here("src","__setup_yeastomics__.r"))
source(here::here("analysis","function_evorate_fitting.R"))

col_ens = setNames(nm =c('uniprot','ensp','ensg'),
                   object=c('uniprotswissprot','ensembl_peptide_id','ensembl_gene_id'))

# Reference Identifiers --------------------------------------------------------
hs_uniprot = get_uniprot_reference()
hs_map_uni = get.uniprot.mapping(9606)
hs_uniref = hs_uniprot$AC
hs_enspref = str_subset(hs_uniprot$ensp, ENSEMBL.nomenclature())
hs_hgnc = load.hgnc(with_protein = F, all_fields = F) %>% filter(uni %in% hs_uniref)
hs_ensgref = drop_na(hs_hgnc,ensg) %>% pull(ensg)

# Sequences --------------------------------------------------------------------
hs_prot = get.uniprot.proteome(9606,DNA = F)
hs_cdna = get.uniprot.proteome(9606,DNA = T)

hs_len = get.width(hs_prot) %>% rename(uniprot=orf,prot_len=len) %>%
  left_join(get.width(hs_cdna), by=c('uniprot'='orf')) %>% rename(cdna_len=len)

hs_aa_freq = (letterFrequency(hs_prot,as.prob = T,letters = get.AA1()) * 100) %>% bind_cols( id=names(hs_prot))
hs_aa_count = letterFrequency(hs_prot,as.prob = F,letters = get.AA1()) %>% bind_cols( id=names(hs_prot))
hs_cdna_gc = (100*rowSums(letterFrequency(hs_cdna, letters="CG",as.prob = T))) %>% round(digits = 2)
hs_codon_freq = Biostrings::trinucleotideFrequency(hs_cdna,step = 3,as.prob = T) *100

# Proteome of reference  (Ensembl and Uniprot) ---------------------------------
hs_ref = left_join(hs_uniprot,hs_hgnc, by=c('AC'='uni')) %>%
         arrange(AC,ensg,ensp,GN) %>%
         mutate(
           uniprot=AC,
           is_uniref = AC %in% hs_uniref,
           has_ensp = ensp %in% hs_enspref,
           has_ensg = ensg %in% hs_ensgref,
           has_cdna = !is.na(id_cdna),
           gene_group = fct_explicit_na(gene_group,na_level = 'unannotated'),
           locus_group = fct_explicit_na(locus_group,na_level = 'unannotated'),
           locus_type = fct_explicit_na(locus_type,na_level = 'unannotated'),
           name = if_na(name,NAME)
           ) %>%
        left_join(hs_len, by='uniprot') %>%
        left_join(tibble(uniprot=names(hs_cdna),UP.GC_cdna  =hs_cdna_gc ),by='uniprot')


# Genomics (%GC, and chromosome number) ----------------------------------------
hs_gc =  get_ensembl_gc(hs_ref$ensg) %>%
         rename(ENS.GC_gene=percentage_gene_gc_content)

hs_chr = get_ensembl_chr(remove_patches = F,ENSG=hs_ref$ensg)

# Transcriptomics --------------------------------------------------------------
hs_transcript = get_ensembl_tx(ENSG=hs_ref$ensg,ENSP=hs_ref$ensp)

id_cols = c("OS","uniprot","AC","ID","GN","NAME","SV","ensg","enst","ensp","ensp_canonical")
uni_col  = c("DB","OX","PE","has_cdna","id_cdna","is_uniref","has_ensp","has_ensg","cdna_len","prot_len")
hgnc_col = c("symbol","name","location","gene_group","locus_group","locus_type")
ens_col  = c("canonical","is_canonical","is_enspref","has_ensp_canonical",
             "has_introns",'n_proteins','n_transcripts','n_exons','n_exons_mini',
             paste0('chr_',c(1:22,'X','Y','MT','other')),
             "cds_len","transcript_len","gene_len")

HS_CODING = left_join(hs_ref,hs_transcript, by='ensg', suffix=c('','_canonical')) %>%
            left_join(hs_gc) %>%
            left_join(hs_chr) %>%
            relocate(all_of(id_cols),all_of(uni_col),all_of(hgnc_col),all_of(ens_col)) %>%
            dplyr::rename_with(all_of(uni_col),.fn = Pxx, 'UP') %>%
            dplyr::rename_with(all_of(hgnc_col),.fn = Pxx, 'HGNC') %>%
            dplyr::rename_with(all_of(ens_col),.fn = Pxx, 'ENS')

# Codons -----------------------------------------------------------------------
#hs_codons = read_delim("/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/codonR/CODON-COUNTS/9606_hs-uniprot.ffn")
library(coRdon)
codon_table = get_codon_table()
hs_codon = bind_cols(uniprot=names(hs_cdna),hs_codon_freq) %>%
  rename(all_of(set_names(codon_table$CODON,codon_table$codon_aa))) %>%
  rename_with(-uniprot, .fn=Pxx, 'UP_cdna')

# Add amino acid with its associated codons
human_codon_usage.rds = here::here('output','hs-codon-usage.rds')
hs_CU = preload(human_codon_usage.rds, load.codon.usage(cds=hs_cdna,with.counts=F,sp = 'hsa') ) %>%
        dplyr::rename_with(starts_with('CU_'),.fn=Pxx,'elek2022')

# Single amino-acid frequencies ------------------------------------------------
hs_aa = hs_aa_freq %>% rename_with(.fn = Pxx, px="UP_prot.f",s='_',.cols=-id)
# Grouped amino-acid frequencies -----------------------------------------------
aa.prop = seqinr::SEQINR.UTIL$AA.PROPERTY %>%
  append( list('Alcohol'=c('S','T'),
               'Turnlike'=c('A','C','D','E','G','H','K','N','Q','R','S','T'))
  )
sum.aa.fr = function(BS,alphabet){
  alphabet.freq = rowSums(letterFrequency(BS,as.prob = T,letters = alphabet))
  prot.enrich = tibble(id=names(BS),fr=alphabet.freq)
}
aa.set = map_chr(aa.prop,str_c,collapse='')
aa.class = paste0(tolower(sub(x=names(aa.set),"\\.","")),"_",aa.set)
AACLASS.FR = map(aa.prop, sum.aa.fr, BS=hs_prot )
hs_aa_class = AACLASS.FR %>% purrr::reduce(full_join,by='id') %>%
              rename_with(~aa.class, starts_with('fr')) %>%
              rename_with(.fn = Pxx, px="UP_prot.f",s='_',.cols=-id)

HS_COUNT = left_join(hs_aa,hs_aa_class) %>% left_join(hs_codon,by=c('id'='uniprot')) %>%
            dplyr::rename(uniprot=id) %>% relocate(uniprot)
# Conservation -----------------------------------------------------------------
#hs_r4s = load.evorate(resdir = "/data/benjamin/Evolution/HUMAN",ref = NULL,ext.r4s = '.r4s')
hs_r4s = read_rds(here("data","RESIDUE-EVORATE-HUMAN.rds")) %>% group_by(id) %>%
             summarize(lmsa=max(msa_pos), lref=max(ref_pos),
                       pid=mean(matched==total-1), pdiv=mean(mismatched>0),pgap=mean(indel>0),
                       nsp = mean(total),
                       rate = mean_(r4s_rate),
                       rate_norm = abs(2*min_(r4s_rate)) + mean_(r4s_rate)) %>%
             dplyr::filter(!is.na(rate)) %>%
             dplyr::rename_with(-id, .fn = Pxx, 'r4s_mammals')


#plot(density(x = hs_r4s_prot$r4s_mammals))

# Disorder ---------------------------------------------------------------------
hs_d2p2 =  read_rds(here("data","d2p2-human-uniprotKB.rds")) %>%
        get.d2p2.diso(.,as.df = T) %>%
        mutate(d2p2.seg = find.consecutive(d2p2.diso>=7, TRUE, min=3),
               d2p2.gap = find.consecutive(d2p2.diso>=7, FALSE, min=1)) %>%
        group_by(d2p2.seg) %>% mutate( d2p2.seglen = sum_(d2p2.seg!=0)) %>%
        group_by(d2p2.gap) %>% mutate( d2p2.gaplen = sum_(d2p2.gap!=0)) %>%
        dplyr::filter(has.d2p2) %>%
        dplyr::select(-c(has.d2p2,d2p2.size)) %>%
        group_by(d2p2.id) %>%
        summarise(D2P2.diso_len = sum_(d2p2.diso>=7),
                  D2P2.diso_frac = mean_(d2p2.diso>=7),
                  D2P2.diso_seg.count = n_distinct(d2p2.seg),
                  D2P2.diso_segmax.len = max(d2p2.seglen)) %>%
        mutate(uniprot = str_extract(d2p2.id,UNIPROT.nomenclature())) %>%
        dplyr::select(-d2p2.id)

hs_pepstats = load.dubreuil2019.data(8) %>%
              dplyr::select(UNIPROT,PEPSTATS.netcharge=netcharge,PEPSTATS.mw=MW,PEPSTATS.pI=pI,UP.prot_len=prot.size) %>%
              group_by(UNIPROT) %>% mutate(PEPSTATS.mean_MW = PEPSTATS.mw/UP.prot_len) %>%
              dplyr::filter(UP.prot_len > 50) %>%
              mutate(PEPSTATS.AA_costly = PEPSTATS.mean_MW > 118, PEPSTATS.AA_cheap=PEPSTATS.mean_MW <= 105)

### Amino-acid interactions propensities
# fullsti = tibble( uni=AA.COUNT$id,
#                   uniprot.sti_full = AACOUNT2SCORE(COUNT=AA.COUNT[,1:20],SCORE=get.stickiness()) )
# df.aascales = DUB %>%
#   mutate(has_iup = !is.na(L.IUP20) | !is.na(L.IUP30) | !is.na(L.IUP40),
#          has_dom = !is.na(L.domain) ) %>%
#   dplyr::filter(has_iup | has_dom) %>%
#   dplyr::select(UNIPROT,ends_with(c('.dom','.iup20'),ignore.case = F)) %>%
#   rename_with(.fn=str_replace, pattern="(.+)\\.(.+)",replacement = 'dubreuil2019.\\1_\\2') %>%
#   left_join(fullsti,by=c('UNIPROT'='uni'))

# Domains ----------------------------------------------------------------------
hs_pfam=load.pfam(tax = 9606) %>%
             filter(seq_id %in% hs_uniref) %>%
             group_by(clan) %>% mutate(pfam.clansize = n_distinct(seq_id)) %>%
              group_by(seq_id) %>% mutate(pfam.ndom=n_distinct(hmm_acc)) %>%
              add_count(hmm_acc,name="pfam.repeat")

#dup_dom = hs_pfam %>% filter(duplicated(seq_id,hmm_acc)) %>% pull(seq_id)
# pfam_dup = hs_pfam %>% filter(is.dup(seq_id)) %>% arrange(seq_id,alignment_start,alignment_end)
# pf_ol_list = pbmcapply::pbmclapply(X= unique(pfam_dup$seq_id),
#                               FUN = calculate_pfam_overlap,
#                              pf=pfam_dup, verbose=F,
#                               mc.cores =  14)
#
# pf_ol = pf_ol_list %>% bind_rows()
#summary(pf_ol)
#pf_ol %>% filter( seq_id %in% pf_ol$seq_id[pf_ol$overlap & pf_ol$no_overlap])

hs_pfam_count = hs_pfam %>% group_by(seq_id) %>%
             summarize( pfam.HMM_none=pfam.ndom==0,
                        pfam.HMM_single=pfam.ndom==1,
                        pfam.HMM_pair=pfam.ndom==2,
                        pfam.HMM_multi=pfam.ndom>=3)%>% distinct()
hs_pfam_dom = pivot_wider(hs_pfam %>% mutate(pfam_val=T), id_cols=seq_id,
                          names_from = 'hmm_name',names_prefix = 'pfam.dom_',
                          values_from = 'pfam_val', values_fill = F, values_fn = sum ) %>%
               mutate(across(where(is.integer), as.logical)) %>%
              dplyr::select(seq_id,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()

hs_pfam_clan = pivot_wider(hs_pfam %>% mutate(pfam_val=T), id_cols=seq_id,
                          names_from = 'clan_name',names_prefix = 'pfam.clan_',
                          values_from = 'pfam_val', values_fill = F, values_fn = sum ) %>%
              mutate(across(where(is.integer), as.logical)) %>%
              dplyr::select(seq_id,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()

HS_PFAM = left_join(hs_pfam_count,hs_pfam_dom) %>% left_join(hs_pfam_clan) %>% distinct()

hs_supfam=load.superfamily(tax = 'hs') %>%
          dplyr::rename(seqid = "sequence_id") %>%
          janitor::clean_names() %>%
          dplyr::filter(seqid %in% hs_ref$ensp) %>%
          group_by(seqid) %>% mutate(superfamily.ndom=n_distinct(superfamily_id)) %>%
          add_count(superfamily_id,name="pfam.repeat")

hs_supfam_count = hs_supfam %>% group_by(seqid) %>%
                 summarize(superfamily.supfam_none=superfamily.ndom==0,
                 superfamily.supfam_single=superfamily.ndom==1,
                 superfamily.supfam_pair=superfamily.ndom==2,
                 superfamily.supfam_multi =superfamily.ndom>=3)
hs_superfamilies = pivot_wider(hs_supfam %>% mutate(supfam_val=T), id_cols=seqid,
                          names_from = 'superfamily_description',names_prefix = 'superfamily.SF_',
                          values_from = 'supfam_val', values_fill = F, values_fn = sum ) %>%
                  mutate(across(where(is.integer), as.logical)) %>%
                  dplyr::select(seqid,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()

hs_families = pivot_wider(hs_supfam %>% mutate(supfam_val=T), id_cols=seqid,
                               names_from = 'family_description',names_prefix = 'superfamily.F_',
                               values_from = 'supfam_val', values_fill = F, values_fn = sum ) %>%
  mutate(across(where(is.integer), as.logical)) %>%
  dplyr::select(seqid,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()


HS_SUPFAM = left_join(hs_supfam_count,hs_superfamilies) %>% left_join(hs_families) %>% distinct()


# Folding energy and stability -------------------------------------------------
hs_stab = load.leuenberger2017.data("Human HeLa Cells",rawdata = F) %>%
          add_count(protein_id,name='npep') %>%
          distinct() %>%
          mutate( nres = round(0.01*protein_coverage*length),
                  Tm_stable = protinfo=="Stable",
                  Tm_medium = protinfo=="Medium",
                  Tm_unstable = protinfo=="Unstable") %>%
          rename_with(.fn=Pxx, px='leuenberger2017.LIP',s='_',.cols=-protein_id) %>%
          dplyr::select(-c(leuenberger2017.LIP_protinfo,
                        leuenberger2017.LIP_protein_coverage,
                        leuenberger2017.LIP_length,
                        leuenberger2017.LIP_measured_domains,
                        leuenberger2017.LIP_theoretical_number_of_domain,
                        leuenberger2017.LIP_nres))


hs_tm = load.jarzab2020.data(org = "H.sapiens") %>%
        dplyr::select(UNIPROT,GENENAME,Tm_celsius,Tm_type,AUC) %>%
         group_by(UNIPROT,GENENAME) %>%
         mutate(Tm_mean=mean_(Tm_celsius),Tm_max=max_(Tm_celsius))

tm_low_celsius = max_(hs_tm$Tm_celsius[hs_tm$Tm_type == 'Low-Tm'])
tm_med_celsius = max_(hs_tm$Tm_celsius[hs_tm$Tm_type == 'Medium-Tm'])
tm_nomelt_celsius =max_(hs_tm$Tm_celsius[hs_tm$Tm_type =='Non-melter' ])

hs_tm = hs_tm %>% mutate(Tm_low = Tm_max < tm_low_celsius,
                         Tm_med = between(Tm_max,tm_low_celsius,tm_med_celsius),
                         Tm_high = Tm_max > tm_med_celsius,
                         Tm_nomelt = is.na(Tm_max)) %>%
        dplyr::select(-c(Tm_type,AUC,Tm_celsius)) %>% distinct() %>%
        rename_with(.fn=Pxx, px='jarzab2020.TPP',s='_',.cols=-c(UNIPROT,GENENAME))

HS_FOLD = left_join(hs_tm,hs_stab, by=c('UNIPROT'='protein_id'))
# Complexes --------------------------------------------------------------------
nmers= paste0(c('mono','di','tri','tetra','penta','hexa','septa','octa','nona','deca'),"mer")
hs_complex = load.meldal.2019.data(species = 'human') %>%
  filter(is_uniprot) %>% # BASED ON UNIPROT
  mutate(oligomers = cut(n_members, breaks = c(1:10,20,81),
                         labels = paste0("meldal2019.CPX_",c(nmers[2:10],'high_oligomer','molecular_machine'))),
         CPLX_ASSEMBLY = gsub("(.+)(\\.$)","\\1",CPLX_ASSEMBLY),
         CPLX_ASSEMBLY = str_replace_all(CPLX_ASSEMBLY," ","_"),
         CPLX_NAME = str_replace_all(CPLX_NAME," ","_")) %>%
  mutate(oligomers_val=T)

hs_oligomers = pivot_wider(hs_complex, id_cols = members, names_from = 'oligomers',
                           names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum,values_fill = F) %>%
  mutate(across(where(is.integer), as.logical)) %>%
  dplyr::select(uniprot=members,where(~ is.logical(.x) && sum(.x) > 9 )) %>%
  dplyr::rename(meldal2019.unknown_oligomer = meldal2019.NA)

hs_assembly = pivot_wider(hs_complex, id_cols = members,
                          names_from = 'CPLX_ASSEMBLY', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum,values_fill = F) %>%
  mutate(across(where(is.integer), as.logical)) %>%
  dplyr::select(uniprot=members,where(~ is.logical(.x) && sum(.x) > 9 )) %>%
  dplyr::rename(meldal2019.unknown_assembly = meldal2019.NA)

hs_complexes = pivot_wider(hs_complex, id_cols = members,
                           names_from = 'CPLX_NAME', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum, values_fill = F) %>%
                mutate(across(where(is.integer), as.logical)) %>%
                dplyr::select(uniprot=members,where(~ is.logical(.x) && sum(.x) > 4 ))
HS_COMPLEX = left_join(hs_oligomers,hs_assembly) %>% left_join(hs_complexes)

# Functional interactions ------------------------------------------------------

hs_string = load.string(tax="9606",phy=F, ful=T, min.score = 900) %>%
  mutate(ens1 = str_extract(protein1,ENSEMBL.nomenclature()),
         ens2 = str_extract(protein2,ENSEMBL.nomenclature())
  ) %>% relocate(ens1,ens2) %>% dplyr::select(-c(protein1,protein2))

human_centrality.rds = here::here('output','hs-string-centralities.rds')
hs_string_centralities = preload(human_centrality.rds, network.centrality(hs_string %>% dplyr::select(ens1,ens2)))


hs_string_centralities$megahub_func = hs_string_centralities$cent_deg >= 300
hs_string_centralities$superhub_func = between(hs_string_centralities$cent_deg,100,300)
hs_string_centralities = hs_string_centralities %>% type_convert() %>%
  dplyr::rename_with(.cols = -ids, .fn = str_replace_all, pattern='string.',replacement="") %>%
  dplyr::rename_with(-ids,.fn=Pxx, 'STRING')
#Biological pathways -----------------------------------------------------------

human_pathways.rds = here('output','hs-kegg-pathways.rds')
hs_pathways=preload(human_pathways.rds,
                    get.KEGG(sp='hsa',type='pathway',as.df=T,to_uniprot = T)) %>%
            mutate(desc=str_replace_all(tolower(desc),' ','_'))


kegg_pathways = hs_pathways %>% mutate(pathway_val=T) %>%
                pivot_wider(id_cols=id,names_from = 'desc',names_prefix = 'KEGG.path_',
                     values_from = 'pathway_val', values_fill = F, values_fn = sum ) %>%
                mutate(across(where(is.integer), as.logical)) %>%
                dplyr::select(id,where(~ is.logical(.x) && sum(.x) >9 ))

human_modules.rds = here('output','hs-kegg-modules.rds')
hs_modules=preload(human_modules.rds,
                   get.KEGG(sp='hsa',type='module',as.df=T,to_uniprot = T)) %>%
  mutate(desc=str_replace_all(tolower(desc),' ','_'))

kegg_modules = hs_modules %>% mutate(module_val=T) %>%
  pivot_wider(id_cols=id,names_from = 'desc',names_prefix = 'KEGG.mod_',
              values_from = 'module_val', values_fill = F, values_fn = sum ) %>%
          mutate(across(where(is.integer), as.logical)) %>%
          dplyr::select(id,where(~ is.logical(.x) && sum(.x) >9 ))

HS_KEGG = left_join(kegg_pathways,kegg_modules)

# Integrate all datasets -------------------------------------------------------
load(here::here('output','hs_datasets.rdata'))
load(here::here('output','hs_integrated_datasets.rdata'))

# Based on uniprot accession
head(hs_r4s) # id = uniprot AC
hs_orthologs = unique(hs_r4s$id)

#head(hs_codon) # ID = uniprot AC
#head(hs_aa) # id = uniprot AC
#head(hs_aa_class) # id = uniprot AC
#head(hs_CU) # ID = uniprot AC
#head(hs_d2p2) # id = uniprot AC

# head(hs_pfam) # seq_id = uniprot AC
# head(hs_pfam_dom) # seq_id = uniprot AC
# head(hs_pfam_clan) # seq_id = uniprot AC

#head(hs_stab) # protein_id = uniprot AC
#head(hs_tm) # UNIPROT = uniprot AC; jarzab2020.TPP_GENENAME = gene symbol

#head(hs_oligomers) # members = uniprot AC
#head(hs_assembly) # members = uniprot AC
#head(hs_complexes) # members = uniprot AC


#head(hs_pathways) # uniprot = uniprot AC
#head(hs_modules)  # uniprot = uniprot AC
# Based on Ensembl peptide identifiers
# head(hs_supfam) # seqid = ensembl peptide
# head(hs_superfamilies)  # seqid = ensembl peptide
# head(hs_families)  # seqid = ensembl peptide
#head(hs_string_centralities) # ids = ensembl peptide
dim(HS_CODING)
dim(HS_COUNT)
dim(hs_CU)
dim(hs_pepstats)
dim(HS_PFAM)
dim(hs_d2p2)
dim(hs_string_centralities)
dim(HS_SUPFAM)
dim(HS_KEGG)
dim(HS_FOLD)
dim(HS_COMPLEX)


HS_DATA = left_join(HS_CODING,HS_COUNT,by=c('uniprot')) %>%
  left_join(hs_CU,by=c('uniprot'='ID')) %>%
  left_join(hs_pepstats,by=c('uniprot'='UNIPROT')) %>%
  left_join(HS_PFAM,by=c('uniprot'='seq_id')) %>%
  left_join(hs_d2p2,by=c('uniprot')) %>%
  left_join(HS_SUPFAM,by=c('ensp'='seqid')) %>%
  left_join(hs_string_centralities,by=c('ensp'='ids')) %>%
  left_join(HS_KEGG,by=c('uniprot'='id')) %>%
  left_join(HS_FOLD,by=c('uniprot'='UNIPROT')) %>%
#  left_join(hs_stab,by=c('uniprot'='protein_id')) %>%
  left_join(HS_COMPLEX,by=c('uniprot')) %>%
  left_join(hs_r4s,by=c('uniprot'='id')) %>%
  ungroup %>%
  distinct() #%>%
  #relocate(uniprot,UP.is_uniref,ensg,ensp,gname, GENENAME, gene_biotype)


# Replace missing values  ------------------------------------------------------

# Replace NA
HS_DATA_nona = HS_DATA %>%
  mutate( across( where(is.logical) & starts_with('kegg.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('pfam.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('superfamily.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('PEPSTATS.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('meldal2019.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('jarzab2020.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('leuenberger2017.LIP'), ~replace_na(., F)) ) %>%
  #mutate( across( where(is.logical) & starts_with('paxdb.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('ENS.'), ~replace_na(., F)) ) %>%
  mutate( ENS.canonical = fct_explicit_na(ENS.canonical,'ensp') ) %>%
  mutate( ENS.is_canonical = fct_explicit_na(ENS.canonical,'0.5') ) %>%
  mutate( across( where(is.logical) & starts_with('string.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('UP.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('HGNC.'), ~replace_na(., F)) ) %>%
  distinct()


IDCOLS = c("uniprot","ensp","ensg","GENENAME","gene_biotype")
#source(here::here("analysis","function_evorate_fitting.R"))
#test = HS_DATA %>% mutate( across(where(is.logical), .fns = ~replace_na(.,F) )) %>% distinct()

miss0 = check_missing_var(HS_DATA)
miss1 = check_missing_var(HS_DATA_nona)

#HS_DATA_nona[,IDCOLS] %>% summarize(across(everything(), sum.na))

# Features selection (and fixing missing values) --------------------------
# Remove rare variables (less than 2 occurrences in orthologs)
HS_DATA_pred1 = remove_rare_vars(df=HS_DATA_nona,min_obs=2)
miss_pred1 = check_missing_var(HS_DATA_pred1)
# Use lower stringencies for string centrality + random forest imputation for proteins with unknown interactions
HS_DATA_pred2 = fix_missing_centrality(df=HS_DATA_pred1, id='ensp', col_prefix="STRING.", taxon=9606)
HS_DATA_pred2 = HS_DATA_pred2 %>% mutate(across(starts_with("STRING.cent_"), .fns = min_))
miss_pred2 = check_missing_var(HS_DATA_pred2)

HS_DATA_pred3 = HS_DATA_pred2 %>%
                mutate(across(starts_with('jarzab2020'), fn=replace_na(.,replace = mean_(.))))

miss_pred3 = check_missing_var(HS_DATA_pred3)


# Abundance -------------------------------------------------------------------
# hs_ppm_uni = readRDS(here::here("data","paxdb_integrated_human.rds"))
hs_uni = hs_map_uni %>% filter(extdb=='UniProtKB-ID') %>% dplyr::select(uni,extid)
hs_paxdb_datasets = find_paxdb_datasets(9606)

# raw paxdb datasets
hs_ppm = load.paxdb(9606,rm.zero=T) %>%
            pivot_wider(id_cols = c(protid,id_uniprot), names_from = c('id'),
                       values_from=ppm, values_fn = mean_) %>%
            left_join(hs_uni, by=c('id_uniprot'='extid')) %>%
            left_join(hs_r4s,by=c('uni'='id')) %>%
            mutate(log10.rate = log10(r4s_mammals.rate_norm))

# integrated paxdb datasets
 hs_paxdb = get.paxdb(tax=9606,abundance = 'integrated',rm.zero=T)
HS_PPM = hs_paxdb %>%
          group_by(protid,id_uniprot) %>%
          summarize( PPM_MIN_ORGAN = min_(ppm_int),
                   PPM_MAX_ORGAN = max_(ppm_int),
                   PPM_AVG_ORGAN = mean_(ppm_int),
                   PPM_MD_ORGAN = median_(ppm_int),
                   PPM_geomAVG_ORGAN = geomean(ppm_int)) %>%
          mutate( across(starts_with('PPM_'), log10) ) %>%
          left_join(hs_uni, by=c('id_uniprot'='extid')) %>%
          left_join( pivot_wider(hs_paxdb, id_cols = c(protid,id_uniprot,n_data,n_int),
                           names_from = 'organ', names_prefix = 'PPM_',
                           values_from = 'ppm_int', values_fn = log10) )
#           left_join(hs_r4s,by=c('uni'='id')) %>%
#           mutate(log10.rate = log10(r4s_mammals.rate_norm)) %>%
#           relocate(protid,id_uniprot,uniprot=uni,n_data,n_int)
#
#
# ppm_ER = HS_PPM %>% dplyr::select(starts_with('PPM_')) %>% colnames %>%
#   bind_cols( map_dfr(., function(x){ spearman.toplot(X=hs_ER$log10.rate, Y=hs_ER[[x]]) }) ) %>%
#   rename(ppm_dataset ='...1') %>% arrange(desc(abs(estimate)),N) %>%
#   dplyr::select(ppm_dataset,N,r=estimate,p.value,) %>% mutate(R2 = 100*r^2)


save(list = ls(pattern = '^hs_'), file = here('output','hs_datasets.rdata'))
save(list = ls(pattern = '^HS_'), file = here('output','hs_integrated_datasets.rdata'))

# Orthologs dataset (with/out abundance) ----------------------------------
all_orthologs = HS_DATA_pred2 %>%
  filter(!is.na(r4s_mammals.rate) & uniprot %in% hs_r4s$id) %>%
  distinct()
dim(all_orthologs)

validation = all_orthologs %>% filter( !(uniprot %in% HS_PPM$uni) )
orthologs =  all_orthologs %>% filter(uniprot %in% HS_PPM$uni)
dim(validation)
dim(orthologs)

# Evo Rate vs. Expression -------------------------------------------------
all_ppm_er_cor = map_dfr(hs_paxdb_datasets$id ,
                         function(x){ spearman.toplot(X=hs_ppm$log10.rate, Y=hs_ppm[[x]]) }
)  %>% add_column(paxdb = hs_paxdb_datasets$id) %>%
  arrange(desc(abs(estimate))) %>%
  left_join(hs_paxdb_datasets, by=c('paxdb'='id')) %>%
  dplyr::select(organ,is_integrated,estimate,N,p.value,filename,cov,yr) %>%
  mutate(filename = ifelse(is_integrated,
                           str_remove_all(filename,"(9606\\-|\\-integrated\\.txt)"),
                           str_remove_all(filename,"(9606\\-|\\.txt|_Maxquant|_gene|SEQUEST)")),
         cor_range = cut(estimate,breaks=c(-0.5,-0.2,-0.1,0.1,0.2,0.5)),
         filename_n = paste0(filename," (N ",N,")")) %>%
  print(n=200)

p_ER_int = ggplot(all_ppm_er_cor %>% filter(is_integrated),
                  aes(y=reorder(filename_n,abs(estimate)),fill=organ,x=estimate)) +
  geom_col() +
  #geom_text(aes(label=N,x=1.2*max(estimate))) +
  geom_text(aes(label=round(estimate,3),x=estimate),size=3,hjust='inward') +
  theme(legend.position='none')

ggsave(p_ER_int, filename=here::here('plots','hs-abundance-evolution-tissues-integrated.pdf'))

p_ER = ggplot(all_ppm_er_cor %>% filter(!is_integrated),
              aes(y=reorder(filename_n,abs(estimate)),fill=organ,x=estimate)) +
  geom_col() +
  #geom_text(aes(label=N,x=1.2*max(estimate)),size=3) +
  geom_text(aes(label=round(estimate,3),x=estimate),size=2,hjust='inward') +
  facet_grid(organ~cor_range,scales = 'free_y',drop = T) +
  theme(legend.position='none',axis.text = element_text(size=6), strip.text = element_text(size=5))

ggsave(p_ER, width = 20, height=40, filename=here::here('plots','hs-abundance-evolution-tissues-experimental.pdf'))


# Dataset for prediction --------------------------------------------------
hs_predictors = HS_DATA_pred2 %>% dplyr::select(where(is.numeric) | where(is.logical)) %>% colnames
nonpred = setdiff(colnames(HS_DATA_pred2),predictor_vars)
check_missing_var(orthologs[,hs_predictors]) #%>% drop_na(-leuenberger2017.tm_protein))
# dim(orthologs)
mean(colnames(orthologs) %in% predictor_vars )

#save.image(here::here('output','hs_predictors_proteome.rdata'))
#predictor_vars = predictor_vars()

HS_DATA_pred2$PPM = log10(HS_DATA_pred2$paxdb.ppm_wholeorg+min_ppm)
PREDICTORS= HS_DATA_pred2
PREDICTORS_nona = PREDICTORS %>% drop_na
check_missing_var(PREDICTORS_nona)

XCOL="PPM"
YCOL="r4s_mammals.rate"
ZCOL=""

fit0 = fit_m0(orthologs,XCOL,YCOL,PREDICTORS_nona,ZCOL,IDCOLS)
formula_M0=reformulate(response = YCOL, termlabels = XCOL, intercept = T )
LM0 = lm(data=PREDICTORS_nona, formula_M0)
df0=decompose_variance(LM0,T)
fit0=list(LM0=LM0,P=PREDICTORS_nona)
spearman.toplot(PREDICTORS$r4s_mammals.rate,PREDICTORS$PPM)

n_pred_vars = colnames(fit0$P) %>% str_extract("cat_.+") %>% setdiff(NA) %>% n_distinct()

P_best= select_variable(fit0,response = '.resid', min_ess = 1, min_ess_frac = 1)
best_pred = fit0$P[,c(XCOL, YCOL, y_resid, ZCOL, P_best$variable)]
n_best = n_distinct(P_best$variable)

hs_evo_ppm = left_join(hs_r4s,HS_PPM, by=c('id'='uniprot')) %>% drop_na
hs_evo_ppm %>% filter(is.dup(id))
spearman.toplot(hs_evo_ppm$r4s_mammals.rate_norm,hs_evo_ppm$paxdb.ppm_int)































# Additional datasets (human specific) ------------------------------------
load.mcshane2016.data = function(){
  # Cell 2016 mcshane E.
  # Kinetic Analysis of Protein Stability Reveals Age-Dependent Degradation
  #Table S4. All AHA Pulse-Chase and Related Data for the Human RPE-1 Cells
  S4="https://www.cell.com/cms/10.1016/j.cell.2016.09.015/attachment/ed36d0dc-573c-4f65-ba99-6ce9713d2584/mmc4.xlsx"
  temp = tempfile()
  download.file(S4,destfile = temp)
  deg_rate = readxl::read_xlsx(path = S4)
  rio::import(S4)
  unlink(temp)
}

load.kristensen2013.data = function(){
  # Mol. Sys. Bio. 2013
  # Protein synthesis rate is the predominant regulator of protein expression during differentiation
  S3 = "https://www.embopress.org/action/downloadSupplement?doi=10.1038%2Fmsb.2013.47&file=msb201347-sup-0004.xlsx"
  temp = tempfile()
  download.file(S3,destfile = temp)

  open.url(S3)
  options(internet.info = 0)
  x=rvest::session("https://www.embopress.org/doi/full/10.1038/msb.2013.47")
  rvest::session_follow_link(S3)
  deg_rate = rio::import(S3)
}
load.bossi2009.data = function(){
  # Mol. Sys. Bio. 2009
  # Tissue specificity and the human protein interaction network
  # https://doi.org/10.1038/msb.2009.17
  interactome="https://www.embopress.org/action/downloadSupplement?doi=10.1038%2Fmsb.2009.17&file=msb200917-sup-0002.zip"
  temp = tempfile()
  test = download.file(interactome)

  tmp <- tempfile()
  curl::curl_download(interactome, tmp)
}


load.soonhengtan2018.data=function(){
  # Science 2018
  # Thermal proximity coaggregation for system-wide profiling of protein complex dynamics in cells
  # DOI: 10.1126/science.aan0346
  Stab = "https://www.science.org/doi/suppl/10.1126/science.aan0346/suppl_file/tabless1_to_s27.zip"
  download_file(Stab,'~/test')
}
load.schroeter2018.data= function(){
  # FASEB 2018, C.B.Schroeter et al
  # Protein half-life determines expression of proteostatic networks in podocyte differentiation
  # https://doi.org/10.1096/fj.201701307R
  halflife = "https://faseb.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1096%2Ffj.201701307R&file=fsb2fj201701307r-sup-0010.xlsx"
  read.xlsx( xlsxFile = halflife, sheet = 1)
}


