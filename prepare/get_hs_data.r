source(here::here("src","__setup_yeastomics__.r"))
col_ens = setNames(nm =c('uniprot','ensp','ensg'),
                   object=c('uniprotswissprot','ensembl_peptide_id','ensembl_gene_id'))

# Sequences --------------------------------------------------------------------
hs_prot = get.uniprot.proteome(9606,DNA = F)
hs_cdna = get.uniprot.proteome(9606,DNA = T)
hs_gc_cdna = (100*rowSums(letterFrequency(hs_cdna, letters="CG",as.prob = T))) %>% round(digits = 2)
hs_uniref = names(hs_prot)

hs_map_uni = get.uniprot.mapping(9606)
hs_uni2ensp = hs_map_uni %>% filter(extdb == 'Ensembl_PRO') %>%
              dplyr::rename(uniprot=uni,ensp=extid) %>%
              dplyr::select(-sp,-upid,-extdb) %>%
              mutate(is_uniref = uniprot %in% hs_uniref)

# Proteome of reference --------------------------------------------------------
hs_ref = hs_uni2ensp %>% separate(col='ensp',into = c('ensp','vers'), sep='\\.') %>%
  dplyr::select(-vers) %>% distinct() %>%
  inner_join(get_ensembl_hsprot() %>% dplyr::rename(all_of(col_ens)) , by=c('ensp','uniprot')) %>%
  filter(chromosome_name %in% c(1:23,"X","Y","MT"))

# Genomics (%GC, and chromosome number) ----------------------------------------
hs_gc = get_hs_GC() %>% as_tibble %>% dplyr::rename(all_of(col_ens)) %>%
        rename(ensembl.GC_gene=percentage_gene_gc_content) %>%
        left_join( tibble(uniprot=names(hs_cdna),uniprot.GC_cdna = hs_gc_cdna), by=c('uniprot'))

hs_chr = get_hs_chr() %>%dplyr::rename(all_of(col_ens)) %>% dplyr::select(-ensp) %>% distinct()

# Reference identifiers for human proteome (Ensembl and Uniprot) ---------------

hs_len = get.width(hs_prot) %>% rename(uniprot=orf, uniprot.prot_len = len ) %>%
         left_join(get.width(hs_cdna), by=c('uniprot'='orf')) %>% rename(uniprot.cdna_len = len )

hs_transcript = get_ensembl_hs(longest_transcript = T) %>% dplyr::rename(all_of(col_ens)) %>%
  dplyr::select(ensg,ensp,uniprot, cds_length,transcript_length,
                n_exons,n_exons_mini,has_introns) %>%
  distinct() %>% left_join(hs_len, by=c('uniprot'))



# Codons -----------------------------------------------------------------------
#hs_codons = read_delim("/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/codonR/CODON-COUNTS/9606_hs-uniprot.ffn")
# library(coRdon)
codon_table = get_codon_table()
hs_codon = Biostrings::trinucleotideFrequency(hs_cdna,step = 3) %>% as_tibble %>%
  rename(all_of(set_names(codon_table$CODON,codon_table$codon_aa)))
# Add amino acid with its associated codons
hs_CU=load.codon.usage(cds=hs_cdna,with.counts=F,sp = 'hsa') %>%
   dplyr::rename_with(.fn=Pxx,px='coRdon',s='.',.cols=starts_with('CU_'))

# Single amino-acid frequencies ------------------------------------------------
AA.FR = letterFrequency(hs_prot,as.prob = T,letters = get.AA1()) %>% bind_cols( id=names(hs_prot))
AA.COUNT = letterFrequency(hs_prot,as.prob = F,letters = get.AA1())%>% bind_cols( id=names(hs_prot))
hs_aa = AA.FR %>% rename_with(.fn = Pxx, px="uniprot.f",s='_',.cols=-id)

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
# Abundance --------------------------------------------------------------------
hs_ppm_uni = readRDS(here::here("data","paxdb_integrated_human.rds"))

# Conservation -----------------------------------------------------------------
#hs_r4s = load.evorate(resdir = "/data/benjamin/Evolution/HUMAN",ref = NULL,ext.r4s = '.r4s')
hs_r4s = read_rds(here("data","RESIDUE-EVORATE-HUMAN.rds")) %>% group_by(id) %>%
             summarize(lmsa=max(msa_pos), lref=max(ref_pos),
                       pid=mean(matched==total-1), pdiv=mean(mismatched>0),pgap=mean(indel>0),
                       nsp = mean(total),
                       r4s_mammals = abs(min(r4s_rate)) + mean_(r4s_rate)) %>%
             dplyr::filter(!is.na(r4s_mammals))

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
        mutate(id = str_extract(d2p2.id,UNIPROT.nomenclature()))

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
             add_count(seq_id,name="pfam.ndom") %>%
             mutate(
                   pfam.HMM_none=pfam.ndom==0,
                   pfam.HMM_single=pfam.ndom==1,
                   pfam.HMM_pair=pfam.ndom==2,
                   pfam.HMM_multi=pfam.ndom>=3)

hs_pfam_dom = pivot_wider(hs_pfam %>% mutate(pfam_val=T), id_cols=seq_id,
                          names_from = 'hmm_name',names_prefix = 'pfam.dom_',
                          values_from = 'pfam_val', values_fill = F, values_fn = sum ) %>%
              dplyr::select(where(~ is.numeric(.x) && sum(.x) >9 ))
hs_pfam_clan = pivot_wider(hs_pfam %>% mutate(pfam_val=T), id_cols=seq_id,
                          names_from = 'clan_name',names_prefix = 'pfam.clan_',
                          values_from = 'pfam_val', values_fill = F, values_fn = sum ) %>%
  dplyr::select(where(~ is.numeric(.x) && sum(.x) >9 ))


hs_supfam=load.superfamily(tax = 'hs') %>%
          dplyr::rename(seqid = "sequence_id") %>%
          janitor::clean_names() %>%
          dplyr::filter(seqid %in% hs_ref$ensp) %>%
          dplyr::add_count(seqid,name="superfamily.ndom") %>%
          mutate(superfam.supfam_none=superfamily.ndom==0,
                 superfam.supfam_single=superfamily.ndom==1,
                 superfam.supfam_pair=superfamily.ndom==2,
                 superfam.supfam_multi =superfamily.ndom>=3)

hs_superfamilies = pivot_wider(hs_supfam %>% mutate(supfam_val=T), id_cols=seqid,
                          names_from = 'superfamily_description',names_prefix = 'superfamily.SF_',
                          values_from = 'supfam_val', values_fill = F, values_fn = sum ) %>%
                  dplyr::select(where(~ is.numeric(.x) && sum(.x) >9 ))

hs_families = pivot_wider(hs_supfam %>% mutate(supfam_val=T), id_cols=seqid,
                               names_from = 'family_description',names_prefix = 'superfamily.F_',
                               values_from = 'supfam_val', values_fill = F, values_fn = sum ) %>%
  dplyr::select(where(~ is.numeric(.x) && sum(.x) >9 ))


# Folding energy and stability -------------------------------------------------
hs_stab = load.leuenberger2017.data("Human HeLa Cells",rawdata = F) %>%
          add_count(protein_id,name='npep') %>%
          distinct() %>%
          mutate( nres = round(0.01*protein_coverage*length),
                  Tm_stable = protinfo=="Stable",
                  Tm_medium = protinfo=="Medium",
                  Tm_unstable = protinfo=="Unstable") %>%
          rename_with(.fn=Pxx, px='leuenberger2017.LIP',s='_',.cols=-protein_id)


hs_tm = load.jarzab2020.data(org = "H.sapiens") %>%
        dplyr::select(UNIPROT,GENENAME,Tm_celsius,Tm_type,AUC) %>%
        mutate(Tm_low = Tm_type == 'Low-Tm',
               Tm_med = Tm_type == 'Medium-Tm',
               Tm_hi = Tm_type == 'High-Tm',
               Tm_nomelt =  Tm_type == 'Non-melter' ) %>%
        rename_with(.fn=Pxx, px='jarzab2020.TPP',s='_',.cols=-c(UNIPROT,GENENAME))


# Complexes --------------------------------------------------------------------

nmers= paste0(c('mono','di','tri','tetra','penta','hexa','septa','octa','nona','deca'),"mer")
hs_complex = load.meldal.2019.data(species = 'human') %>%
  filter(is_uniprot) %>% # BASED ON UNIPROT
  mutate(oligomers = cut(n_members, breaks = c(1:10,20,81),
                         labels = paste0("meldal2019.CPX_",c(nmers[2:10],'high_oligomer','molecular_machine'))),
         CPLX_ASSEMBLY = gsub("(.+)(\\.$)","\\1",CPLX_ASSEMBLY) ) %>%
  mutate(oligomers_val=T)

hs_oligomers = pivot_wider(hs_complex, id_cols = members, names_from = 'oligomers',
                           values_from = oligomers_val, values_fn=sum,values_fill = F)
hs_assembly = pivot_wider(hs_complex, id_cols = members,
                          names_from = 'CPLX_ASSEMBLY', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum,values_fill = F)
hs_complexes = pivot_wider(hs_complex, id_cols = members,
                           names_from = 'CPLX_NAME', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum,values_fill = F) %>%
               dplyr::select(where(~ is.numeric(.x) && sum(.x) >9 ))

# Functional interactions ------------------------------------------------------

hs_string = load.string(tax="9606",phy=F, ful=T, min.score = 900) %>%
  mutate(ens1 = str_extract(protein1,ENSEMBL.nomenclature()),
         ens2 = str_extract(protein2,ENSEMBL.nomenclature())
  ) %>% relocate(ens1,ens2) %>% dplyr::select(-c(protein1,protein2))

hs_string_centralities = network.centrality(hs_string %>% dplyr::select(ens1,ens2))
hs_string_centralities$string.megahub_func = hs_string_centralities$cent_deg >= 300
hs_string_centralities$string.superhub_func = between(hs_string_centralities$cent_deg,100,300)

#%>%
#   mutate(string.megahub_func = )
# INTERACTIONS$string.megahub_func = cent.STRING %>% dplyr::filter(cent_deg>300) %>% pull(ids)
# INTERACTIONS$string.superhub_func = cent.STRING %>% dplyr::filter(between(cent_deg,100,300)) %>% pull(ids)


#Biological pathways -----------------------------------------------------------
hs_pathways=get.KEGG(sp='hsa',type='pathway',as.df=T,to_uniprot = T)
hs_modules=get.KEGG(sp='hsa',type='module',as.df=T,to_uniprot = T)


# Integrate all datasets -------------------------------------------------------
save(list = ls(pattern = '^hs_'), file = here('output','hs_datasets.rdata'))
load(here('output','hs_datasets.rdata'))

hs_ref
#hs_ens
hs_ens2uni
hs_uni2ens
hs_chr
hs_gc
hs_transcript

# Based on uniprot accession
head(hs_r4s) # id = uniprot AC
head(hs_codon) # ID = uniprot AC
head(hs_CU) # ID = uniprot AC
head(hs_aa) # id = uniprot AC
head(hs_aa_class) # id = uniprot AC
head(hs_d2p2) # id = uniprot AC
head(hs_pfam) # seq_id = uniprot AC
head(hs_pfam_dom) # seq_id = uniprot AC
head(hs_pfam_clan) # seq_id = uniprot AC

head(hs_stab) # protein_id = uniprot AC
head(hs_tm) # UNIPROT = uniprot AC; jarzab2020.TPP_GENENAME = gene symbol
head(hs_oligomers) # members = uniprot AC
head(hs_assembly) # members = uniprot AC
head(hs_complexes) # members = uniprot AC
head(hs_pathways) # uniprot = uniprot AC
head(hs_modules)  # uniprot = uniprot AC
# Based on Ensembl peptide identifiers
head(hs_ppm_uni) # protid = ensembl peptide ; uniprot = uniprot AC ; id_uniprot = uniprot NAME ; GENENAME = gene symbol
head(hs_supfam) # seqid = ensembl peptide
head(hs_superfamilies)  # seqid = ensembl peptide
head(hs_families)  # seqid = ensembl peptide
head(hs_string_centralities) # ids = ensembl peptide
