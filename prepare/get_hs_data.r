load(here::here('output','hs_datasets.rdata'))
load(here::here('output','hs_integrated_datasets.rdata'))
load(here::here('output','hs_predictors_proteome.rdata'))

source(here::here("src","__setup_yeastomics__.r"))
source(here::here("analysis","function_evorate_fitting.R"))

col_ens = setNames(nm =c('uniprot','ensp','ensg'),
                   object=c('uniprotswissprot','ensembl_peptide_id','ensembl_gene_id'))

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


# Sequences --------------------------------------------------------------------
hs_prot = get.uniprot.proteome(9606,DNA = F)
hs_cdna = get.uniprot.proteome(9606,DNA = T)
hs_uniref = names(hs_prot)

hs_aa_freq = (letterFrequency(hs_prot,as.prob = T,letters = get.AA1()) * 100) %>% bind_cols( id=names(hs_prot))
hs_aa_count = letterFrequency(hs_prot,as.prob = F,letters = get.AA1()) %>% bind_cols( id=names(hs_prot))
hs_cdna_gc = (100*rowSums(letterFrequency(hs_cdna, letters="CG",as.prob = T))) %>% round(digits = 2)
hs_codon_freq = Biostrings::trinucleotideFrequency(hs_cdna,step = 3,as.prob = T) *100

hs_map_uni = get.uniprot.mapping(9606)
hs_uni2ensp = hs_map_uni %>%
       filter(extdb %in% c('Gene_Name','Ensembl','Ensembl_PRO') ) %>%
       mutate(row = dense_rank(uni)) %>%
       tidyr::pivot_wider(id_cols = c(uni,row), names_from='extdb', values_from='extid', values_fn = list ) %>%
       unnest_longer(col = 'Gene_Name') %>%
       unnest_longer(col = 'Ensembl') %>%
       unnest_longer(col = 'Ensembl_PRO') %>%
       separate(col = Ensembl_PRO, sep = '\\.',into=c('ensp','ensp_v')) %>%
       separate(col = Ensembl, sep = '\\.',into=c('ensg','ensg_v')) %>%
       dplyr::select(-row,-ensp_v,-ensg_v) %>%
       dplyr::rename(uniprot=uni, gname=Gene_Name) %>%
       mutate(is_uniref = uniprot %in% hs_uniref) %>%
       distinct()

# Proteome of reference  (Ensembl and Uniprot) ---------------------------------
hs_ref = hs_uni2ensp %>%
         left_join(get_ensembl_hsprot() %>% dplyr::rename(all_of(col_ens[1:2])) , by=c('ensp','uniprot')) %>%
         filter(is_uniref & gene_biotype == 'protein_coding') %>%
         group_by(uniprot) %>% mutate(n_ensg = n_distinct(ensg), n_ensp = n_distinct(ensp)) %>%
         arrange(desc(n_ensg), ensg, desc(n_ensp), ensp, uniprot, gname)

# Genomics (%GC, and chromosome number) ----------------------------------------
hs_gc = get_hs_GC() %>% as_tibble %>% dplyr::rename(all_of(col_ens)) %>%
        rename(ensembl.GC_gene=percentage_gene_gc_content) %>%
        left_join( tibble(uniprot=names(hs_cdna),uniprot.GC_cdna = hs_cdna_gc), by=c('uniprot'))

hs_chr = get_hs_chr(remove_patches = F) %>% dplyr::rename(all_of(col_ens)) %>%
          dplyr::select(-ensp) %>% distinct()

# Sequences Length -------------------------------------------------------------

hs_transcript = get_ensembl_hs(longest_transcript = T) %>% dplyr::rename(all_of(col_ens)) %>%
  dplyr::select(ensg,ensp,uniprot, cds_length,transcript_length,
                n_exons,n_exons_mini,has_introns) %>%
  distinct()


HS_CODING = left_join(hs_ref,hs_transcript) %>%
            left_join(hs_gc) %>%
            left_join(hs_chr) %>%
            relocate(uniprot,ensp,is_uniref,ensg,gene_biotype,ensembl.GC_gene,uniprot.GC_cdna,
                     cds_length,transcript_length,
                     has_introns,n_exons,n_exons_mini,
                     paste0('chr_',c(1:22,'X','Y','MT'))) %>%
            group_by(uniprot) %>% filter( row_number() == nearest(value=F,x=max(transcript_length),y=transcript_length,n = 1)) %>%
            relocate(uniprot,is_uniref,gname,ensg,ensp,gene_biotype) %>%
            dplyr::rename_with(-c(uniprot:gene_biotype, starts_with('uniprot.'), starts_with('ensembl.')),.fn = Pxx, 'ensembl')

# Codons -----------------------------------------------------------------------
#hs_codons = read_delim("/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/codonR/CODON-COUNTS/9606_hs-uniprot.ffn")
# library(coRdon)
hs_len = get.width(hs_prot) %>% rename(uniprot=orf, uniprot.prot_len = len ) %>%
  left_join(get.width(hs_cdna), by=c('uniprot'='orf')) %>% rename(uniprot.cdna_len = len )

codon_table = get_codon_table()
hs_codon = bind_cols(uniprot=names(hs_cdna),hs_codon_freq) %>%
  rename(all_of(set_names(codon_table$CODON,codon_table$codon_aa))) %>%
  rename_with(-uniprot, .fn=Pxx, 'ensembl')

# Add amino acid with its associated codons
hs_CU=load.codon.usage(cds=hs_cdna,with.counts=F,sp = 'hsa') %>%
   dplyr::rename_with(starts_with('CU_'),.fn=Pxx,'elek2022')

# Single amino-acid frequencies ------------------------------------------------
hs_aa = hs_aa_freq %>% rename_with(.fn = Pxx, px="uniprot.f",s='_',.cols=-id)
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
              rename_with(.fn = Pxx, px="uniprot.f",s='_',.cols=-id)

HS_COUNT = left_join(hs_aa,hs_aa_class) %>% left_join(hs_codon,by=c('id'='uniprot')) %>%
            dplyr::rename(uniprot=id) %>% left_join(hs_len,by='uniprot')
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
        summarise(d2p2.L = sum_(d2p2.diso>=7),
                  d2p2.f = mean_(d2p2.diso>=7),
                  d2p2.nseg = n_distinct(d2p2.seg),
                  d2p2.Lsegmax = max(d2p2.seglen)) %>%
        mutate(id = str_extract(d2p2.id,UNIPROT.nomenclature())) %>%
        dplyr::select(-d2p2.id)

hs_pepstats = load.dubreuil2019.data(8) %>%
              dplyr::select(UNIPROT,pepstats.netcharge=netcharge,pepstats.mw=MW,pepstats.pI=pI,uniprot.prot_size=prot.size) %>%
              group_by(UNIPROT) %>% mutate(pepstats.mean_MW = pepstats.mw/uniprot.prot_size) %>%
              dplyr::filter(uniprot.prot_size > 50) %>%
              mutate(pepstats.AA_costly = pepstats.mean_MW > 118, pepstats.AA_cheap=pepstats.mean_MW <= 105)

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
  dplyr::select(members,where(~ is.logical(.x) && sum(.x) > 9 )) %>%
  dplyr::rename(meldal2019.unknown_oligomer = meldal2019.NA)

hs_assembly = pivot_wider(hs_complex, id_cols = members,
                          names_from = 'CPLX_ASSEMBLY', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum,values_fill = F) %>%
  mutate(across(where(is.integer), as.logical)) %>%
  dplyr::select(members,where(~ is.logical(.x) && sum(.x) > 9 )) %>%
  dplyr::rename(meldal2019.unknown_assembly = meldal2019.NA)

hs_complexes = pivot_wider(hs_complex, id_cols = members,
                           names_from = 'CPLX_NAME', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum, values_fill = F) %>%
                mutate(across(where(is.integer), as.logical)) %>%
                dplyr::select(members,where(~ is.logical(.x) && sum(.x) > 4 ))
HS_COMPLEX = left_join(hs_oligomers,hs_assembly) %>% left_join(hs_complexes)

# Functional interactions ------------------------------------------------------

hs_string = load.string(tax="9606",phy=F, ful=T, min.score = 900) %>%
  mutate(ens1 = str_extract(protein1,ENSEMBL.nomenclature()),
         ens2 = str_extract(protein2,ENSEMBL.nomenclature())
  ) %>% relocate(ens1,ens2) %>% dplyr::select(-c(protein1,protein2))

hs_string_centralities = network.centrality(hs_string %>% dplyr::select(ens1,ens2))
hs_string_centralities$megahub_func = hs_string_centralities$cent_deg >= 300
hs_string_centralities$superhub_func = between(hs_string_centralities$cent_deg,100,300)
hs_string_centralities = hs_string_centralities %>% type_convert() %>%
  dplyr::rename_with(.cols = -ids, .fn = str_replace_all, pattern='string.',replacement="") %>%
  dplyr::rename_with(-ids,.fn=Pxx, 'string')
#Biological pathways -----------------------------------------------------------
hs_pathways=get.KEGG(sp='hsa',type='pathway',as.df=T,to_uniprot = T) %>%
            mutate(desc=str_replace_all(tolower(desc),' ','_'))

kegg_pathways = hs_pathways %>% mutate(pathway_val=T) %>%
                pivot_wider(id_cols=id,names_from = 'desc',names_prefix = 'kegg.path_',
                     values_from = 'pathway_val', values_fill = F, values_fn = sum ) %>%
                mutate(across(where(is.integer), as.logical)) %>%
                dplyr::select(id,where(~ is.logical(.x) && sum(.x) >9 ))

hs_modules=get.KEGG(sp='hsa',type='module',as.df=T,to_uniprot = T) %>%
           mutate(desc=str_replace_all(tolower(desc),' ','_'))

kegg_modules = hs_modules %>% mutate(module_val=T) %>%
  pivot_wider(id_cols=id,names_from = 'desc',names_prefix = 'kegg.mod_',
              values_from = 'module_val', values_fill = F, values_fn = sum ) %>%
          mutate(across(where(is.integer), as.logical)) %>%
          dplyr::select(id,where(~ is.logical(.x) && sum(.x) >9 ))

HS_KEGG = left_join(kegg_pathways,kegg_modules)

# Integrate all datasets -------------------------------------------------------
save(list = ls(pattern = '^hs_'), file = here('output','hs_datasets.rdata'))
save(list = ls(pattern = '^HS_'), file = here('output','hs_integrated_datasets.rdata'))
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
dim(HS_PFAM)
dim(HS_SUPFAM)
dim(HS_COUNT)
dim(hs_CU)
dim(hs_pepstats)
dim(HS_KEGG)
dim(HS_FOLD)
dim(HS_COMPLEX)


HS_DATA = left_join(HS_COUNT,HS_CODING,by=c('uniprot')) %>%
  left_join(hs_pepstats,by=c('uniprot'='UNIPROT')) %>%
  left_join(HS_PFAM,by=c('uniprot'='seq_id')) %>%
  left_join(HS_SUPFAM,by=c('ensp'='seqid')) %>%
  left_join(hs_d2p2,by=c('uniprot'='id')) %>%
  left_join(hs_string_centralities,by=c('ensp'='ids')) %>%
  left_join(hs_CU,by=c('uniprot'='ID')) %>%
  left_join(HS_FOLD,by=c('uniprot'='UNIPROT')) %>%
#  left_join(hs_stab,by=c('uniprot'='protein_id')) %>%
  left_join(HS_COMPLEX,by=c('uniprot'='members')) %>%
  left_join(HS_KEGG,by=c('uniprot'='id')) %>%
  left_join(hs_r4s,by=c('uniprot'='id')) %>%
  distinct() %>%
  relocate(uniprot,is_uniref,ensg,ensp,ensembl.gname, GENENAME, gene_biotype)


# Replace missing values  ------------------------------------------------------

# Replace NA
HS_DATA_nona = HS_DATA %>%
  mutate( across( starts_with('kegg.'), ~replace_na(., F)) ) %>%
  mutate( across( starts_with('pfam.'), ~replace_na(., F)) ) %>%
  mutate( across( starts_with('superfamily.'), ~replace_na(., F)) ) %>%
  mutate( across( starts_with('pepstats.'), ~replace_na(., F)) ) %>%
  mutate( across( starts_with('meldal2019.'), ~replace_na(., F)) ) %>%
  mutate( across( starts_with('jarzab2020.'), ~replace_na(., F)) ) %>%
  mutate( across( starts_with('leuenberger2017.LIP'), ~replace_na(., F)) ) %>%
  #mutate( across( where(is.logical) & starts_with('paxdb.'), ~replace_na(., F)) ) %>%
  mutate( across( starts_with('ensembl.'), ~replace_na(., F)) ) %>%
  mutate( across( where(is.logical) & starts_with('string.'), ~replace_na(., F)) ) %>%
  mutate( is_uniref=replace_na(is_uniref,F)) %>%
  distinct()

#source(here::here("analysis","function_evorate_fitting.R"))
test = HS_DATA %>% mutate( across(where(is.logical), .fns = ~replace_na(.,F) )) %>% distinct()

# Orthologs dataset (with/out abundance) ----------------------------------

all_orthologs = HS_DATA_nona %>%
  filter(!is.na(r4s_mammals.rate) & uniprot %in% hs_orthologs) %>%
  distinct()

## Abundance -------------------------------------------------------------------
# hs_ppm_uni = readRDS(here::here("data","paxdb_integrated_human.rds"))
# HS_PPM = hs_ppm_uni %>%
#   dplyr::select(ensp=protid,uniprot,GENENAME,is_uniref,is_duplicated,
#                 ppm_wholeorg,ppm_max,ppm_int_max,ppm_int) %>%
#   dplyr::rename_with(-c(ensp,uniprot,GENENAME), .fn = Pxx, 'paxdb')
hs_ppm = get.paxdb(tax=9606,abundance = 'integrated',rm.zero=T)

hs_organ = hs_ppm %>% group_by(protid) %>%
           summarize( PPM_MIN_ORGAN = min_(ppm_int),
                   PPM_MAX_ORGAN = max_(ppm_int),
                   PPM_AVG_ORGAN = mean_(ppm_int),
                   PPM_MD_ORGAN = median_(ppm_int),
                   PPM_geomAVG_ORGAN = geomean(ppm_int)) %>%
          mutate( across(starts_with('PPM_'), log10) )


hs_ppm_organ = pivot_wider(hs_ppm, id_cols = c(protid,id_uniprot,n_data,n_int),
                           names_from = 'organ', names_prefix = 'PPM_',
                           values_from = 'ppm_int', values_fn = log10) %>%
               left_join(hs_organ,by='protid')
hs_map_uni = get.uniprot.mapping(9606,'UniProtKB-ID') %>% dplyr::select(uni,extid)
HS_PPM = left_join(hs_ppm_organ,hs_map_uni, by=c('id_uniprot'='extid'))

hs_ER = left_join(hs_r4s,HS_PPM, by=c('id'='uni')) %>%
        mutate(log10.rate = log10(r4s_mammals.rate))

ppm_ER = hs_ER %>% dplyr::select(starts_with('PPM_')) %>% colnames %>%
          bind_cols( map_dfr(., function(x){ spearman.toplot(X=hs_ER$log10.rate, Y=hs_ER[[x]]) }) ) %>%
          rename(ppm_dataset ='...1') %>% arrange(estimate,N) %>%
  dplyr::select(ppm_dataset,N,r=estimate,p.value,) %>% mutate(R2 = 100*r^2)

ppm_ER %>% print(n=40)


validation = all_orthologs %>% filter( !(uniprot %in% HS_PPM$uni) )
orthologs =  all_orthologs %>% filter(uniprot %in% HS_PPM$uni)
dim(validation)
dim(orthologs)


# miss0 = check_missing_var(HS_DATA)
miss1 = check_missing_var(HS_DATA_nona)
miss1.1 = check_missing_var(test)
miss2 = check_missing_var(all_orthologs)

# Features selection (and fixing missing values) --------------------------

# Remove rare variables (less than 2 occurrences in orthologs)
HS_DATA_pred1 = remove_rare_vars(df=orthologs,min_obs=2)
# Use lower stringencies for string centrality + random forest imputation for proteins with unknown interactions
HS_DATA_pred2 = fix_missing_centrality(df=HS_DATA_pred1, id='ensp', col_prefix="string.", taxon=9606)
miss_pred1 = check_missing_var(HS_DATA_pred1)
miss_pred2 = check_missing_var(HS_DATA_pred2)

predictor_vars = HS_DATA_pred2 %>% dplyr::select(where(is.numeric) | where(is.logical)) %>% colnames
nonpred = setdiff(colnames(HS_DATA_pred2),predictor_vars)
#save.image(here::here('output','hs_predictors_proteome.rdata'))

predictor_vars = predictor_vars()

min_ppm  = min(HS_DATA_pred2$paxdb.ppm_wholeorg[HS_DATA_pred2$paxdb.ppm_wholeorg!=0])

HS_DATA_pred2$PPM = log10(HS_DATA_pred2$paxdb.ppm_wholeorg+min_ppm)
PREDICTORS= HS_DATA_pred2
PREDICTORS_nona = PREDICTORS %>% drop_na
check_missing_var(PREDICTORS_nona)

XCOL="PPM"
YCOL="r4s_mammals.rate"
ZCOL=""
IDCOLS = c("uniprot","ensp","ensg","GENENAME","gene_biotype")

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
