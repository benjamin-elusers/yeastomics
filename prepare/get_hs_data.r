#load(here::here('output','hs_datasets.rdata'))
#load(here::here('output','hs_integrated_datasets.rdata'))
#load(here::here('output','hs_predictors_proteome.rdata'))

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
            get.d2p2.diso(as.df=T) %>%
            summarize.d2p2() %>%
            mutate(uniprot = str_extract(d2p2.id, UNIPROT.nomenclature())) %>%
            dplyr::select(-d2p2.id)
# Peptide stats ----------------------------------------------------------------
hs_pepstats = load.dubreuil2019.data(8) %>%
              dplyr::select(UNIPROT,PEPSTATS.netcharge=netcharge,PEPSTATS.mw=MW,PEPSTATS.pI=pI,UP.prot_len=prot.size) %>%
              group_by(UNIPROT) %>% mutate(PEPSTATS.mean_MW = PEPSTATS.mw/UP.prot_len) %>%
              dplyr::filter(UP.prot_len > 50) %>%
              mutate(PEPSTATS.AA_costly = PEPSTATS.mean_MW > 118, PEPSTATS.AA_cheap=PEPSTATS.mean_MW <= 105)

# Dubreuil et al 2019 (Disorder/Stickiness) ------------------------------------
IUP_cols = c("L.IUP20" = 'IUP20_L',"L.IUP30" = 'IUP30_L',"L.IUP40" = 'IUP40_L',
             "f.IUP20" = 'IUP20_f',"f.IUP30" = 'IUP30_f',"f.IUP40" = 'IUP40_f')
subset_cols = c('standard'='iupred_standard','medium'='iupred_medium','high'='iupred_high',
                'average'='uniprot_average','large'='uniprot_long')

iseq = 1:length(hs_prot)
full_score <- pbmcapply::pbmclapply(
  X=iseq,  FUN = function(x){ get_aa_score(string = as.character(hs_prot[[x]]) ) },
  mc.cores = 14,mc.cleanup = T)

hs_fullscore = tibble(uniprot=hs_uniref, length=widths(hs_prot)) %>%
                bind_cols(bind_rows(full_score)) %>%
                group_by(uniprot,length) %>%
                summarize( across(aggrescan:wimleywhite, ~ sum(.x)/length ) ) %>%
                distinct() %>%
                dplyr::rename_with(.cols = -c(uniprot,length), .fn = xxS, 'full',s='_' )


hs_dubreuil = load.dubreuil2019.data(8) %>%
         dplyr::select('UNIPROT',
                      contains('IUP'),
                      ends_with(c('dom','iup40')),
                      c('standard','medium','high','small','average','large'),
                      ) %>%
         dplyr::rename(set_names(names(IUP_cols),IUP_cols),
                        set_names(names(subset_cols),subset_cols)) %>%
         dplyr::rename_with(.cols= ends_with(c('.dom','.iup40')), .fn = str_replace_all, "\\.", "_" ) %>%
         right_join(hs_fullscore,by=c('UNIPROT'='uniprot')) %>%
         dplyr::select(-contains(c('ratioR_RK'))) %>%
         dplyr::rename_with(.cols=-UNIPROT, .fn = Pxx, 'dubreuil2019', s='.')

# Domains ----------------------------------------------------------------------
hs_pfam=load.pfam(tax = 9606) %>%
             filter(seq_id %in% hs_uniref) %>%
             group_by(clan) %>% mutate(pfam.clansize = n_distinct(seq_id)) %>%
              group_by(seq_id) %>% mutate(pfam.ndom=n_distinct(hmm_acc)) %>%
              add_count(hmm_acc,name="pfam.repeat")

#dup_dom = hs_pfam %>% filter(duplicated(seq_id,hmm_acc)) %>% pull(seq_id)
# pfam_dup = hs_pfam %>% filter(is.dup(seq_id)) %>% arrange(seq_id,alignment_start,alignment_end)
# pf_ol_list = pbmcapply::pbmclapply(X= unique(pfam_dup$seq_id),FUN = calculate_pfam_overlap,  pf=pfam_dup, verbose=F,
#                               mc.cores =  14)
#
# pf_ol = pf_ol_list %>% bind_rows()
#summary(pf_ol)
#pf_ol %>% filter( seq_id %in% pf_ol$seq_id[pf_ol$overlap & pf_ol$no_overlap])

hs_pfam_count = left_join(hs_ref,hs_pfam, by=c('AC'='seq_id')) %>%
             group_by(AC) %>%
             summarize( pfam.HMM_none= is.na(pfam.ndom) | pfam.ndom==0,
                        pfam.HMM_single= !is.na(pfam.ndom) & pfam.ndom==1,
                        pfam.HMM_pair= !is.na(pfam.ndom) & pfam.ndom==2,
                        pfam.HMM_multi=!is.na(pfam.ndom) & pfam.ndom>=3) %>% distinct()

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

HS_PFAM = left_join(hs_pfam_count,hs_pfam_dom,by=c('AC'='seq_id')) %>%
          left_join(hs_pfam_clan,by=c('AC'='seq_id')) %>%
          distinct()

hs_supfam=load.superfamily(tax = 'hs') %>%
          dplyr::rename(seqid = "sequence_id") %>%
          janitor::clean_names() %>%
          dplyr::filter(seqid %in% hs_ref$ensp) %>%
          group_by(seqid) %>% mutate(superfamily.ndom=n_distinct(superfamily_id)) %>%
          add_count(superfamily_id,name="pfam.repeat")

hs_supfam_count = left_join(hs_ref,hs_supfam, by=c('AC'='seqid')) %>% group_by(AC) %>%
                 summarize(superfamily.supfam_none= is.na(superfamily.ndom) | superfamily.ndom==0 ,
                 superfamily.supfam_single= !is.na(superfamily.ndom) & superfamily.ndom==1,
                 superfamily.supfam_pair= !is.na(superfamily.ndom) & superfamily.ndom==2,
                 superfamily.supfam_multi = !is.na(superfamily.ndom) & superfamily.ndom>=3)
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

HS_SUPFAM = left_join(hs_supfam_count,hs_superfamilies,by=c('AC'='seqid')) %>%
            left_join(hs_families,by=c('AC'='seqid')) %>%
            distinct()

# Folding energy and stability -------------------------------------------------
hs_stab = load.leuenberger2017.data("Human HeLa Cells",rawdata = F) %>%
          add_count(protein_id,name='npep') %>%
          distinct() %>%
          mutate( nres = round(0.01*protein_coverage*length),
                  Tm_stable = protinfo=="Stable",
                  Tm_medium = protinfo=="Medium",
                  Tm_unstable = protinfo=="Unstable") %>%
          dplyr::select(-c(protinfo,protein_coverage,length,
                           measured_domains,theoretical_number_of_domain,nres)) %>%
          rename_with(.fn=Pxx, px='leuenberger2017.LIP',s='_',.cols=-protein_id)

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

# Biological functions (GO) ----------------------------------------------------
hs_unigo = get.uniprot.go(hs_uniref,taxon = 9606) %>%  # Uniprot-based GO annotation
        mutate( go = paste0(ONTOLOGY,"_", str_replace_all(goterm,pattern="[^A-Za-z0-9\\.]+","_"))) %>%
        filter(!obsolete & shared > 10 & ONTOLOGY != "CC") %>%
        arrange(go)

# Convert to a matrix format (columns = GO term, rows=proteome)
HS_UNIGO  = hs_unigo %>%  mutate(seen=T) %>%
           pivot_wider( id_cols = c('UNIPROT'), names_from = 'go', values_from = 'seen',
               values_fn=list(seen = unique),values_fill = list(seen=F)) %>%
           dplyr::rename_with(.cols = -UNIPROT,.fn = Pxx, 'go.', s='')

# Subcellular locations (Uniprot) ----------------------------------------------
hs_uniloc = query_uniprot_subloc(uniprot = hs_uniref,todf=T) # as a wide dataframe
HS_UNILOC = hs_uniloc %>%
            group_by(id) %>%
            dplyr::select(id, where(~is.numeric(.x) && sum(.x) >50 ) ) %>%
            dplyr::rename_with(-id,.fn=Pxx, 'UP.loc', s='_') %>%
            mutate( across(.cols=everything(),.fns= as.logical))

# Biological pathways (KEGG) ---------------------------------------------------
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

# Based on uniprot accession
head(hs_r4s) # id = uniprot AC
hs_orthologs = unique(hs_r4s$id)
dim(HS_CODING)
colnames(HS_CODING)
dim(HS_COUNT)
dim(hs_dubreuil)
dim(HS_PFAM)
dim(HS_SUPFAM)
dim(hs_CU)
dim(hs_pepstats)
dim(hs_d2p2)
dim(HS_UNIGO)
dim(HS_UNILOC)
dim(hs_string_centralities)
dim(HS_KEGG)
dim(HS_FOLD)
dim(HS_COMPLEX)

HS_DATA = left_join(HS_CODING,HS_COUNT,by=c('uniprot')) %>%
  left_join(hs_dubreuil,by=c('uniprot'='UNIPROT')) %>%
  left_join(HS_PFAM,by=c('uniprot'='AC')) %>%
  left_join(HS_SUPFAM,by=c('ensp'='AC')) %>%
  left_join(hs_CU,by=c('uniprot'='ID')) %>%
  left_join(hs_pepstats,by=c('uniprot'='UNIPROT','UP.prot_len')) %>%
  left_join(hs_d2p2,by=c('uniprot')) %>%
  left_join(HS_UNIGO,by=c('uniprot'='UNIPROT')) %>%
  left_join(HS_UNILOC,by=c('uniprot'='id')) %>%
  left_join(hs_string_centralities,by=c('ensp'='ids')) %>%
  left_join(HS_KEGG,by=c('uniprot'='id')) %>%
  left_join(HS_FOLD,by=c('uniprot'='UNIPROT')) %>%
  left_join(HS_COMPLEX,by=c('uniprot')) %>%
  left_join(hs_r4s,by=c('uniprot'='id')) %>%
  ungroup %>%
  distinct() #%>%
  #relocate(uniprot,UP.is_uniref,ensg,ensp,gname, GENENAME, gene_biotype)

# Replace missing values  ------------------------------------------------------
miss0 = check_missing_var(HS_DATA)


#### 1. Fix logical variables (NA replaced by FALSE) ####
HS_FEATURES.1 = HS_DATA %>%
  mutate( across( where(is.logical) & starts_with('kegg.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('pfam.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('superfamily.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('PEPSTATS.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('dubreuil2019.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('meldal2019.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('jarzab2020.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('leuenberger2017.LIP'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('go.MF'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('go.BP'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('UP.loc'), ~replace_na(., F)) ) %>%
  #mutate( across( where(is.logical) & starts_with('paxdb.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('ENS.'), ~replace_na(., F)) ) %>%
  mutate( ENS.canonical = fct_explicit_na(ENS.canonical,'ensp') ) %>%
  mutate( ENS.is_canonical = fct_explicit_na(ENS.canonical,'0.5') ) %>%
  mutate( across( where(is.logical) & starts_with('string.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('UP.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('HGNC.'), ~replace_na(., F)) ) %>%
  distinct()

miss1 = check_missing_var(HS_FEATURES.1)

#### 2. Remove rare variables  (less than 2 occurrences in proteome) ####
HS_FEATURES.2 = remove_rare_vars(df=HS_FEATURES.1,min_obs=2)
miss2 = check_missing_var(HS_FEATURES.2)

#### 3. Fix network centrality (lower confidence + random forest) ####
# Take a long time (@! ~50mn for human !@)
HS_FEATURES.3 = fix_missing_centrality(df=HS_FEATURES.2, id='ensp', col_prefix="STRING.", taxon=9606)
HS_FEATURES.3 = HS_FEATURES.3 %>% mutate(across(starts_with("STRING.cent_"), .fns = min_))

miss3 = check_missing_var(HS_FEATURES.3)

#### 4. Fix thermodynamcis stability (replace NA by average) ####
HS_FEATURES.4 = HS_FEATURES.3 %>%
                 mutate(across(starts_with(c('jarzab2020','leuenberger2017')),
                               ~coalesce(as.numeric(.x), mean_(get(cur_column())))))

miss4 = check_missing_var(HS_FEATURES.4)

#### 5. Fix PEPSTATS missing values (i.e not included in Dubreuil et al. 2019) ####
HS_FEATURES.5 = fix_missing_peptide_stats(df=HS_FEATURES.4, id='uniprot', taxon=9606,
                                         col_len='UP.prot_len',
                                         col_mw='PEPSTATS.mw',
                                         col_mw_avg='PEPSTATS.mean_MW',
                                         col_charge='PEPSTATS.netcharge',
                                         col_pi='PEPSTATS.pI')
miss5 = check_missing_var(HS_FEATURES.5)

#### 6. Fix D2P2 missing values (fetch from d2p2 or use consensus prediction from mobiDB) ####

# Get all uniprot or ensembl reference identifiers
uni_na_d2p2 = HS_FEATURES.5 %>% filter(is.na(D2P2.diso_len)) %>%  pull(uniprot)
ensp_na_d2p2 = HS_FEATURES.5 %>% filter(is.na(D2P2.diso_len) & !is.na(ensp)) %>% pull(ensp)
#n_distinct(uni_na_d2p2)
#n_distinct(ensp_na_d2p2)
id_na_d2p2 = c(uni_na_d2p2,ensp_na_d2p2)

#  Fetch d2p2 missing identifiers
d2p2_na =  load.d2p2(id_na_d2p2, 'output/hs-missing-d2p2.rds') %>%
        get.d2p2.diso(as.df=T) %>%
        summarize.d2p2() %>%
        mutate(uniprot = str_extract(d2p2.id, UNIPROT.nomenclature())) %>%
        mutate(ensp = str_extract(d2p2.id, ENSEMBL.nomenclature()),
               has_d2p2=T) %>%
        dplyr::select(-d2p2.id)

# Get mobiDB data for missing identifiers
mobidb_d2p2 = fetch.mobidb(id_na_d2p2)
# Filter for consensus disorder prediction (50% agreement) and without d2p2 data
mobidb_d2p2_na = mobidb_d2p2 %>%
        filter(feature=='disorder' & source=='th_50' & !(acc %in% c(d2p2_na$ensp,d2p2_na$uniprot)) ) %>%
        group_by(acc) %>%
        # To combine d2p2 and mobidb they must have absolutely identical column names
        mutate(D2P2.diso_seg.count = n(), D2P2.diso_segmax.len = max(feature_len),
             D2P2.diso_len=sum(feature_len), D2P2.diso_frac=sum(feature_len)/length,
             uniprot=acc, has_d2p2=F) %>%
        ungroup() %>%
        dplyr::select(-c(S,E,feature_len,evidence,feature,source,content_fraction,content_count,length,acc)) %>%
        distinct()
# Combine d2p2 and mobidb data
d2p2_nona = bind_rows(d2p2_na,mobidb_d2p2_na) %>%
          mutate( uniprot=ifelse(is.na(uniprot) | uniprot=="P00000",NA,uniprot))
uni_d2p2_nona =  d2p2_nona %>% filter(!is.na(uniprot)) %>% dplyr::select(-ensp)
ensp_d2p2_nona =  d2p2_nona %>% filter(!is.na(ensp)) %>% dplyr::select(-uniprot)

# Replace the missing d2p2 values by the d2p2/mobidb fixed data
HS_FEATURES.6 = coalesce_join(HS_FEATURES.5,uni_d2p2_nona,by=c('uniprot')) %>%
                coalesce_join(ensp_d2p2_nona,by=c('ensp')) %>%
                mutate(has_d2p2 = ifelse(is.na(has_d2p2),T,has_d2p2) )

miss6 = check_missing_var(HS_FEATURES.6)

#save.image(here('output','checkpoint-hs-data.rdata'))
#load(here::here('output','checkpoint-hs-data.rdata'))
#### 7. Fix missing amino acid propensities ####
uni_na_aascore = HS_FEATURES.6 %>%
                    dplyr::select(AC,starts_with('dubreuil2019')) %>%
                    dplyr::filter(row_number() %in% find_na_rows(.,as.indices = T)) %>%
                    pull(AC) %>% unique
seq_aa_score = hs_prot[uni_na_aascore]
na_aascore = retrieve_missing_aascore(uni_na_aascore,seq_aa_score)
df_na_aascore = na_aascore %>%
                dplyr::rename_with(.cols=-acc, .fn=str_replace_all, 'pawar_', 'pawar_ph7_' ) %>%
                dplyr::select(-contains(c('voronoi'))) %>%
                dplyr::rename_with(.cols=-acc, .fn=Pxx, 'dubreuil2019', s='.' ) %>%
                dplyr::rename(uniprot='acc')


iup40_cols = str_subset(colnames(HS_FEATURES.6),'^dubreuil.+iup40$') %>% sort
dom_cols = str_subset(colnames(HS_FEATURES.6),'^dubreuil.+dom$') %>% sort
full_cols = str_subset(colnames(HS_FEATURES.6),'^dubreuil.+full$') %>% sort %>% setdiff("dubreuil2019.voronoi_sickiness_full")

HS_FEATURES.7 = coalesce_join(HS_FEATURES.6,df_na_aascore,by='uniprot') %>%
                rowwise() %>%
                mutate(dubreuil2019.stickiness_iup40  = coalesce(dubreuil2019.stickiness_iup40   , dubreuil2019.stickiness_full)) %>%
                mutate(dubreuil2019.aggrescan_iup40   = coalesce(dubreuil2019.aggrescan_iup40    , dubreuil2019.aggrescan_full     )) %>%
                mutate(dubreuil2019.foldamyloid_iup40 = coalesce(dubreuil2019.foldamyloid_iup40  , dubreuil2019.foldamyloid_full   )) %>%
                mutate(dubreuil2019.pawar_ph7_iup40   = coalesce(dubreuil2019.pawar_ph7_iup40    , dubreuil2019.pawar_ph7_full     )) %>%
                mutate(dubreuil2019.camsol_iup40      = coalesce(dubreuil2019.camsol_iup40       , dubreuil2019.camsol_full        )) %>%
                mutate(dubreuil2019.wimleywhite_iup40 = coalesce(dubreuil2019.wimleywhite_iup40  , dubreuil2019.wimleywhite_full   )) %>%
                mutate(dubreuil2019.roseman_iup40     = coalesce(dubreuil2019.roseman_iup40      , dubreuil2019.roseman_full       )) %>%
                mutate(dubreuil2019.kytedoolittle_iup40  = coalesce(dubreuil2019.kytedoolittle_iup40, dubreuil2019.kytedoolittle_full )) %>%
                mutate(dubreuil2019.stickiness_dom    = coalesce(dubreuil2019.stickiness_dom     , dubreuil2019.stickiness_full    )) %>%
                mutate(dubreuil2019.aggrescan_dom     = coalesce(dubreuil2019.aggrescan_dom      , dubreuil2019.aggrescan_full     )) %>%
                mutate(dubreuil2019.foldamyloid_dom   = coalesce(dubreuil2019.foldamyloid_dom    , dubreuil2019.foldamyloid_full   )) %>%
                mutate(dubreuil2019.pawar_ph7_dom     = coalesce(dubreuil2019.pawar_ph7_dom      , dubreuil2019.pawar_ph7_full     )) %>%
                mutate(dubreuil2019.camsol_dom        = coalesce(dubreuil2019.camsol_dom         , dubreuil2019.camsol_full        )) %>%
                mutate(dubreuil2019.wimleywhite_dom   = coalesce(dubreuil2019.wimleywhite_dom    , dubreuil2019.wimleywhite_full   )) %>%
                mutate(dubreuil2019.roseman_dom       = coalesce(dubreuil2019.roseman_dom        , dubreuil2019.roseman_full       )) %>%
                mutate(dubreuil2019.kytedoolittle_dom = coalesce(dubreuil2019.kytedoolittle_dom  , dubreuil2019.kytedoolittle_full )) %>%
                mutate(dubreuil2019.IUP20_L = coalesce(dubreuil2019.IUP20_L  , D2P2.diso_len )) %>%
                mutate(dubreuil2019.IUP30_L = coalesce(dubreuil2019.IUP30_L  , D2P2.diso_len )) %>%
                mutate(dubreuil2019.IUP40_L = coalesce(dubreuil2019.IUP40_L  , D2P2.diso_len )) %>%
                mutate(dubreuil2019.IUP20_f = coalesce(dubreuil2019.IUP20_f  , D2P2.diso_frac )) %>%
                mutate(dubreuil2019.IUP30_f = coalesce(dubreuil2019.IUP30_f  , D2P2.diso_frac )) %>%
                mutate(dubreuil2019.IUP40_f = coalesce(dubreuil2019.IUP40_f  , D2P2.diso_frac)) %>%
    ungroup %>%
    dplyr::select(-contains(c('content_count','pawar_ph7')))

miss7 = check_missing_var(HS_FEATURES.7)

#### 8. Fix missing cDNA sequences ####

#cdna_na = HS_FEATURES.7 %>% dplyr::select( all_of(miss7$variable) )
#missing_cdna = HS_FEATURES.7 %>% filter( is.na(UP_cdna.AAA_Lys_K) ) %>% pull(uniprot)
#get.uniprot.proteome(9606,DNA=T)

HS_FEATURES.8 = HS_FEATURES.7 %>%
                mutate( ENS.cds_len = coalesce( ENS.cds_len, UP.cdna_len),
                        ensp = coalesce(ensp, ensp_canonical),
                        GENENAME = coalesce(GENENAME),GN)
miss8 = check_missing_var(HS_FEATURES.8 %>% ungroup)

HS_PROTEOME_DATA = HS_FEATURES.8 %>%
                     # Get out the rows with NA
                     drop_na(starts_with(c('UP_cdna.','elek2022.','D2P2.','ENS.','dubreuil2019.'))) %>%
                     dplyr::select(-GENENAME) %>%
                     # RENAME VARIABLES NON-ALPHANUMERIC CHARACTERS (NOT VALID FOR FORMULA)
                     dplyr::rename_with(.cols = everything(), .fn = str_replace_all, pattern="[^A-Za-z0-9\\.]+", replacement="_")


miss = check_missing_var(HS_PROTEOME_DATA)
dim(HS_PROTEOME_DATA)

# Retrieve human abundance -----------------------------------------------------
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

save(list = ls(pattern = '^hs_'), file = here('output','hs_datasets.rdata'))
save(list = ls(pattern = '^HS_'), file = here('output','hs_integrated_datasets.rdata'))

# Orthologs dataset (with/out abundance) ---------------------------------------
# Full proteome N=20722 - after removing NA = 19417
HS_ORTHOLOGS = HS_PROTEOME_DATA %>%
                filter(!is.na(r4s_mammals.rate) & uniprot %in% hs_r4s$id) %>%
                filter(!is.na(uniprot) & !is.dup(uniprot)) %>%
                # Remove rare variables among orthlogs (min 3 orthologs must share the feature)
                remove_rare_vars(df=.,min_obs=3) %>%
                distinct()

dim(HS_ORTHOLOGS)
# All orthologs N = 12977

# Define two predictions datasets with/out abundance (training/validation)
hs_validation = HS_ORTHOLOGS %>% filter( !(uniprot %in% HS_PPM$uni) )
hs_orthologs =  HS_ORTHOLOGS %>% filter(uniprot %in% HS_PPM$uni)
dim(hs_validation) # n=263
dim(hs_orthologs) # n=12714

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
              aes(y=reorder(filename_n,abs(estimate)),x=organ, fill=estimate)) +
  geom_raster() + scale_fill_viridis_c() + theme(axis.text.x = element_text(angle=90,hjust=1))
  #geom_text(aes(label=N,x=1.2*max(estimate)),size=3) +
  #geom_text(aes(label=round(estimate,3),x=estimate),size=2,hjust='inward') +
  #facet_grid(organ~cor_range,scales = 'free_y',drop = T) +
  #theme(legend.position='none',axis.text = element_text(size=6), strip.text = element_text(size=5))
ggsave(p_ER,scale=3,filename=here::here('plots','hs-abundance-evolution-tissues-experimental.pdf'))


# Explaining variance in evolutionary rate -------------------------------------
#IDCOLS = c("uniprot","ensp","ensg","GENENAME","gene_biotype")
ID_COLS = c('uniprot','GN','ensp')
XCOL="PPM"
YCOL="ER"
ZCOL=NULL

hs_ortho_predictors = hs_orthologs %>% group_by(uniprot,GN,ensp) %>%
  dplyr::select(where(is.numeric) | where(is.logical)) %>%
  dplyr::select(-c('SV', starts_with("r4s_mammals"),"has_unique_id") )

HS_PREDICTORS = hs_ortho_predictors %>% ungroup

HS_PPM$PPM = HS_PPM$PPM_WHOLE_ORGANISM
HS_ER = left_join(hs_r4s,HS_PPM, by=c('id'='uni'))
HS_ER$ER = log10(HS_ER$r4s_mammals.rate_norm)

HS_LMDATA = left_join(HS_PREDICTORS,HS_ER,by=c('uniprot'='id')) %>% drop_na(PPM) %>%
            filter(!is.na(uniprot) & !is.dup(uniprot)) %>%
            distinct()

fit0 = fit_m0(INPUT_LM = HS_LMDATA,PREDICTORS = HS_PREDICTORS,
              XCOL = XCOL,YCOL = YCOL, ZCOL=NULL, IDCOLS=ID_COLS,
              MAX_XCOR = 0.5, MAX_YCOR=0.5, MIN_N = 2)
df0=decompose_variance(fit0$LM0,T)
spearman.toplot(fit0$P$ER,fit0$P$PPM)

# Filter Dataset for prediction --------------------------------------------------
hs_evo_all= select_variable(fit0,response='.resid', raw=T)
saveRDS(hs_evo_all,here("output","lm-evo-hs.rds"))

hs_best= hs_evo_all %>% filter(pc_ess > 0.1 & variable != YCOL)
hs_best= hs_evo_all %>% filter(pc_ess > 0.5 & variable != YCOL)
hs_best= hs_evo_all %>% filter(pc_ess > 1 & variable != YCOL)


hs_best_pred = fit0$P[,c(XCOL, YCOL, '.resid', ZCOL, hs_best$variable)]
hs_nbest = n_distinct(hs_best$variable)

formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM_ER = lm(data=hs_best_pred, formula_null)
decompose_variance(LM_ER)

#formula_hs_best_ER=paste0(YCOL," ~ ",paste0( P_evo$variable,collapse=" + "))
formula_hs_best = reformulate(hs_best$variable, response =  '.resid')

m_best_ER = step(object=LM_ER, scope = as.formula(formula_hs_best), direction = 'forward',k=log(hs_nbest)*2,trace=0)
decompose_variance(m_best_ER,to.df = T)

hs_best_lm = lm(reformulate(response = YCOL, termlabels = labels(m_best_ER)),data=hs_best_pred)
print(labels(hs_best_lm))
print(hs_nbest)
decompose_variance(hs_best_lm,to.df = T)

## Validate
hs_valid = left_join(hs_validation,HS_ER,by=c('uniprot'='id')) %>%
  filter(!is.na(uniprot) & !is.dup(uniprot)) %>%
  distinct()

LM_validation_ER = lm(data=hs_valid, formula_null)
decompose_variance(LM_validation_ER,to.df = T)

hs_validation_lm = lm(reformulate(response = YCOL, termlabels = labels(hs_best_lm)),data=hs_valid)
print(labels(hs_validation_lm))
print(coefficients(hs_validation_lm)[is.na(coefficients(hs_validation_lm))])

decompose_variance(hs_validation_lm,to.df = T)

sc_evo_var = readRDS(here::here('output','sc_features_ess_over_0.5.rds'))
library(stringdist)
sim_var = stringsimmatrix(sc_evo_var,colnames(HS_LMDATA), method='jw',p=0.1,useNames='strings')
maxS= apply(sim_var,1,max_)
i_maxS= apply(sim_var,1,which.max)
matched = tibble(sc_var=rownames(sim_var), hs_var = colnames(sim_var)[i_maxS], similarity=maxS)

write_delim(matched,here::here('output','sc-hs-features.tsv'),delim = '\t')


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


