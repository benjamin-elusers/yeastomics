source(here::here("src","__setup_yeastomics__.r"))
hs_prot = get.uniprot.proteome(9606,DNA = F)
hs_cdna = get.uniprot.proteome(9606,DNA = T)
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
df.CU=load.codon.usage(cds=hs_cdna,with.counts=F,sp = 'hsa')#%>% dplyr::rename_with(.fn=Pxx,px='coRdon',s='.',.cols=starts_with('CU_'))

#hs_r4s = load.evorate(resdir = "/data/benjamin/Evolution/HUMAN",ref = NULL,ext.r4s = '.r4s')
hs_r4s = read_rds(here("data","RESIDUE-EVORATE-HUMAN.rds"))
hs_d2p2 =  read_rds(here("data","d2p2-human-uniprotKB.rds"))
#hs_codons = read_delim("/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/codonR/CODON-COUNTS/9606_hs-uniprot.ffn")
hs_diso = get.d2p2.(hs_d2p2)
D2P2 = get.d2p2.diso(hs_d2p2,as.df = T) %>%
  mutate(d2p2.seg = find.consecutive(d2p2.diso>=7, TRUE, min=3),
         d2p2.gap = find.consecutive(d2p2.diso>=7, FALSE, min=1)) %>%
  group_by(d2p2.seg) %>% mutate( d2p2.seglen = sum_(d2p2.seg!=0)) %>%
  group_by(d2p2.gap) %>% mutate( d2p2.gaplen = sum_(d2p2.gap!=0))
df.d2p2 = D2P2 %>%
  dplyr::filter(has.d2p2) %>%
  dplyr::select(-c(has.d2p2,d2p2.size)) %>%
  group_by(d2p2.id) %>%
  summarise(d2p2.L = sum_(d2p2.diso>=7),
            d2p2.f = mean_(d2p2.diso>=7),
            d2p2.nseg = n_distinct(d2p2.seg),
            d2p2.Lsegmax = max(d2p2.seglen))
#df.pepstats = DUB %>% dplyr::select(UNIPROT,pepstats.netcharge=netcharge,pepstats.mw=MW,pepstats.pI=pI)

AA.FR = letterFrequency(hs_prot,as.prob = T,letters = get.AA1()) %>% bind_cols( id=names(hs_prot))
AA.COUNT = letterFrequency(hs_prot,as.prob = F,letters = get.AA1())%>% bind_cols( id=names(hs_prot))
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
### Single amino-acid frequencies ------------------------------------------------
df.aa = AA.FR %>%
  rename_with(.fn = Pxx, px="uniprot.f",s='_',.cols=-id)


### Grouped amino-acid frequencies -----------------------------------------------
df.aa_class =  AACLASS.FR %>% purrr::reduce(full_join,by='id') %>% rename_with(~aa.class, starts_with('fr')) %>%
  rename_with(.fn = Pxx, px="uniprot.f",s='_',.cols=-id)
names(df.aa_class)

PFAM=load.pfam(tax = 9606)
SUPERFAM=load.superfamily(tax = 'hs') %>% dplyr::rename(seqid = "sequence_id")

DELTAG = load.leuenberger2017.data("Human HeLa Cells",rawdata = F) %>%
  add_count(protein_id,name='npep') %>%
  distinct() %>%
  mutate( nres = round(0.01*protein_coverage*length))

TM = load.jarzab2020.data(org = "H.sapiens") # error: only NAs

STRING = load.string(tax="9606",phy=F, ful=T, min.score = 900) %>%
  mutate(ORF1 = str_extract(protein1,UNIPROT.nomenclature()),
         ORF2 = str_extract(protein2,UNIPROT.nomenclature())
  ) %>% relocate(ORF1,ORF2) %>% dplyr::select(-c(protein1,protein2))

org='H.sapiens'
message("REF: A. Jarzab et al., 2020, Nature Methods")
message("Meltome atlas—thermal proteome stability across the tree of life")
#https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-020-0801-4/MediaObjects/41592_2020_801_MOESM7_ESM.xlsx
F2_url = "https://figshare.com/ndownloader/files/21653313"
species = c("T.thermophilus", "P.torridus", "G.stearothermophilus",
            "E.coli", "B.subtilis",
            "S.cerevisiae", "M.musculus",  "H.sapiens", "C.elegans", "D.melanogaster", "D.rerio",
            "A.thaliana", "O.antarctica")
organisms = match.arg(org,species,several.ok = T)

F2a = rio::import(F2_url,col_names = TRUE,skip=1,sheet=1) %>%
  as_tibble %>%
  dplyr::filter(Species %in% organisms) %>%
  janitor::clean_names()

F2b = rio::import(F2_url,col_names = TRUE,skip=1,sheet=2) %>%
  as_tibble %>%
  mutate(Species = str_replace_all(Dataset," +","") ) %>%
  dplyr::select(-Dataset) %>%
  dplyr::filter(Species %pin% organisms) %>%
  janitor::clean_names()

F2=left_join(F2a,F2b) %>%
  separate(col = protein_id, sep = '_',into = c('UNIPROT','GENENAME')) %>%
  dplyr::select(species,UNIPROT,GENENAME,
                Tm_celsius=melting_point_c,
                Tm_type=protein_classification,
                AUC=area_under_the_melting_curve)

nmers= paste0(c('mono','di','tri','tetra','penta','hexa','septa','octa','nona','deca'),"mer")
CPX = load.meldal.2019.data(species = 'human') %>%
  filter(is_uniprot) %>% # BASED ON UNIPROT
  mutate(oligomers = cut(n_members, breaks = c(1:10,20,81),
                         labels = paste0("meldal2019.CPX_",c(nmers[2:10],'high_oligomer','molecular_machine'))),
         CPLX_ASSEMBLY = gsub("(.+)(\\.$)","\\1",CPLX_ASSEMBLY) )

PATHWAY=get.KEGG('hsa','pathway',as.df=T)
MODULE=get.KEGG('hsa','module',as.df=T)

