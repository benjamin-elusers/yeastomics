source(here::here("src","__setup_yeastomics__.r"))
# Sequences --------------------------------------------------------------------
hs_prot = get.uniprot.proteome(9606,DNA = F)
hs_cdna = get.uniprot.proteome(9606,DNA = T)

# Codons -----------------------------------------------------------------------
#hs_codons = read_delim("/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/codonR/CODON-COUNTS/9606_hs-uniprot.ffn")
library(coRdon)
hs_codons = codonTable(hs_cdna)
# Add amino acid with its associated codons
codon_table = seqinr::SEQINR.UTIL$CODON.AA %>%
  as_tibble() %>%
  mutate(CODON=toupper(CODON), AA=str_to_title(AA),
         L=str_replace(string = L, "\\*", "STOP"),
         codon_aa = paste0(CODON,"_",AA,"_",L))
codon_count = Biostrings::trinucleotideFrequency(hs_cdna,step = 3) %>% as_tibble
col_codons=codon_table$codon_aa[codon_table$CODON == colnames(codon_count)]
names(codon_count) = col_codons
hs_CU=load.codon.usage(cds=hs_cdna,with.counts=F,sp = 'hsa') %>% dplyr::rename_with(.fn=Pxx,px='coRdon',s='.',.cols=starts_with('CU_'))

# Single amino-acid frequencies ----------------------------------------------
AA.FR = letterFrequency(hs_prot,as.prob = T,letters = get.AA1()) %>% bind_cols( id=names(hs_prot))
AA.COUNT = letterFrequency(hs_prot,as.prob = F,letters = get.AA1())%>% bind_cols( id=names(hs_prot))
hs_aa = AA.FR %>% rename_with(.fn = Pxx, px="uniprot.f",s='_',.cols=-id)

# Grouped amino-acid frequencies ---------------------------------------------
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

# Conservation ---------------------------------------------------------------
#hs_r4s = load.evorate(resdir = "/data/benjamin/Evolution/HUMAN",ref = NULL,ext.r4s = '.r4s')
hs_r4s = read_rds(here("data","RESIDUE-EVORATE-HUMAN.rds")) %>% group_by(id) %>%
             summarize(lmsa=max(msa_pos), lref=max(ref_pos),
                       pid=mean(matched==total-1), pdiv=mean(mismatched>0),pgap=mean(indel>0),
                       nsp = mean(total),
                       r4s_mammals = abs(min(r4s_rate)) + mean_(r4s_rate)) %>%
             dplyr::filter(!is.na(r4s_mammals))

#plot(density(x = hs_r4s_prot$r4s_mammals))

# Disorder -------------------------------------------------------------------
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
                  d2p2.Lsegmax = max(d2p2.seglen))

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

# Domains --------------------------------------------------------------------

hs_pfam=load.pfam(tax = 9606) %>%
             group_by(clan) %>% mutate(pfam.clansize = n_distinct(seq_id)) %>%
             add_count(seq_id,name="pfam.ndom") %>%
             mutate(pfam.HMM_none=pfam.ndom==0,
                   pfam.HMM_single=pfam.ndom==1,
                   pfam.HMM_pair=pfam.ndom==2,
                   pfam.HMM_multi=pfam.ndom>=3)

hs_supfam=load.superfamily(tax = 'hs') %>%
          dplyr::rename(seqid = "sequence_id") %>%
          janitor::clean_names() %>%
          dplyr::add_count(seqid,name="superfamily.ndom") %>%
          mutate(superfam.supfam_none=superfamily.ndom==0,
                 superfam.supfam_single=superfamily.ndom==1,
                 superfam.supfam_pair=superfamily.ndom==2,
                 superfam.supfam_multi =superfamily.ndom>=3)


# Folding energy and stability -----------------------------------------------
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
        rename_with(.fn=Pxx, px='jarzab2020.TPP',s='_',.cols=-UNIPROT)


# Complexes ------------------------------------------------------------------

nmers= paste0(c('mono','di','tri','tetra','penta','hexa','septa','octa','nona','deca'),"mer")
hs_complex = load.meldal.2019.data(species = 'human') %>%
  filter(is_uniprot) %>% # BASED ON UNIPROT
  mutate(oligomers = cut(n_members, breaks = c(1:10,20,81),
                         labels = paste0("meldal2019.CPX_",c(nmers[2:10],'high_oligomer','molecular_machine'))),
         CPLX_ASSEMBLY = gsub("(.+)(\\.$)","\\1",CPLX_ASSEMBLY) )


# Functional interactions ----------------------------------------------------

hs_string = load.string(tax="9606",phy=F, ful=T, min.score = 900) %>%
  mutate(ens1 = str_extract(protein1,ENSEMBL.nomenclature()),
         ens2 = str_extract(protein2,ENSEMBL.nomenclature())
  ) %>% relocate(ens1,ens2) %>% dplyr::select(-c(protein1,protein2))

hs_string_centralities = network.centrality(hs_string %>% dplyr::select(ens1,ens2))
#%>%
#   mutate(string.megahub_func = )
# INTERACTIONS$string.megahub_func = cent.STRING %>% dplyr::filter(cent_deg>300) %>% pull(ids)
# INTERACTIONS$string.superhub_func = cent.STRING %>% dplyr::filter(between(cent_deg,100,300)) %>% pull(ids)


#Biological pathways --------------------------------------------------------

hs_pathways=get.KEGG('hsa','pathway',as.df=T)
hs_modules=get.KEGG('hsa','module',as.df=T)


# Integrate all datasets -----------------------------------------------------
save(list = ls(pattern = '^hs_'), file = here('output','hs_datasets.rdata'))
load(here('output','hs_datasets.rdata'))

head(hs_r4s)
head(hs_codons)
head(hs_CU)
head(hs_aa)
head(hs_aa_class)
head(hs_d2p2)
head(hs_pfam)
head(hs_supfam)
head(hs_stab)
head(hs_tm)
head(hs_complex)
head(hs_string_centralities)
head(hs_pathways)
head(hs_modules)
