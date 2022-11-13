#load(here::here('output','sc_datasets.rdata'))
#load(here::here('output','sc_integrated_datasets.rdata'))
#load(here::here('output','sc_predictors_proteome.rdata'))
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("analysis","function_evorate_fitting.R"))

uniprot_cerevisiae = find.uniprot_refprot(c('cerevisiae','s288c','eukaryota'))
TAXID = 559292 # uniprot_cerevisiae$tax_id

# Sequences --------------------------------------------------------------------
# from sgd
s288c_prot = load.sgd.proteome(withORF=T)
s288c_cdna = load.sgd.CDS(withORF = T)
# from uniprot
sc_prot = get.uniprot.proteome(TAXID,DNA = F)
sc_cdna = get.uniprot.proteome(TAXID,DNA = T)

# sgd-uniprot
orf2uni = tibble(orf=names(s288c_prot), uni_matched = names(sc_prot)[match(s288c_prot,sc_prot)])
orf_na = orf2uni$orf[is.na(orf2uni$uni_matched)]
#s288c_prot[orf_na] # 680 ORFs with no known uniprot

uni_with_orf = orf2uni %>% drop_na %>% pull(uni_matched)
sc_uniref = names(sc_prot)

# Reference Identifiers --------------------------------------------------------
sc_uniprot = get_uniprot_reference(TAXID) %>% arrange(AC)
sc_sgd = load.sgd.features(by.chr=T)

sc_map_uniref = get.uniprot.mapping(TAXID) # for reference proteome
sc_map_uni = read_delim("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping.dat.gz",
                        col_names = c('uni','extdb','extid'))
#SC_UNI = get_uniprot_ids(sc_uniref) # Not necessary if all IDs are within the reference proteome

sc_uni2ens = sc_map_uni %>%
              filter(extdb == "EnsemblGenome_PRO" ) %>%
              dplyr::select(AC=uni,ensp=extid)

sc_uni2sgd = sc_map_uni %>%
  filter(extdb == 'SGD') %>%
  dplyr::select(AC=uni,sgd=extid)

sc_uni_mapped = left_join(sc_uni2ens,sc_uni2sgd,by='AC')

sc_reference = inner_join(sc_sgd %>% dplyr::select(sgdid,type,qual,name,gname,chr),coalesce_join(sc_uni_mapped,sc_uniprot, by='AC'),by=c('sgdid'='sgd'))
#sc_orfref = drop_na(sc_sgd,name) %>% filter(name %in% names(s288c_prot)) %>% dplyr::select(sgdid,orf=name)

## REFERENCE PROTEOME SEQUENCES ------------------------------------------------
UNIPROT = sc_reference %>%
  group_by(AC) %>%
  mutate( is_uniref = AC %in% sc_uniref, one2one = (n()==1) ) %>%
  dplyr::rename(ORF = name, GNAME=gname, CHR=chr)

sc_len = get.width(sc_prot) %>% rename(prot_len=len) %>%
  left_join(get.width(sc_cdna), by=c('orf')) %>% rename(uniprot=orf,cdna_len=len)

sc_aa_freq = (letterFrequency(s288c_prot,as.prob = T,letters = get.AA1()) * 100) %>% bind_cols( id=names(s288c_prot))
sc_aa_count = letterFrequency(s288c_prot,as.prob = F,letters = get.AA1()) %>% bind_cols( id=names(s288c_prot))
sc_cdna_gc = (100*rowSums(letterFrequency(s288c_cdna, letters="CG",as.prob = T))) %>% round(digits = 2)
sc_codon_freq = Biostrings::trinucleotideFrequency(s288c_cdna,step = 3,as.prob = T) *100

# Proteome of reference  (Ensembl and Uniprot) ---------------------------------
sc_ref = UNIPROT %>%
         arrange(AC,ensp,GN) %>%
         mutate(
           ensg = ORF,
           uniprot=AC,
           is_uniref = AC %in% sc_uniref,
           has_orf = ORF %in% names(s288c_prot),
           has_ensg = !is.na(ensg) | ensg == "",
           has_ensp = ensp %in% sc_uni_mapped$ensp,
           has_cdna = !is.na(id_cdna),
           ) %>%
        left_join(sc_len, by='uniprot') %>%
        left_join(tibble(ORF=names(s288c_cdna),UP.GC_cdna = sc_cdna_gc ),by='ORF')

sc_biomart = get_ensembl_dataset('fungi_mart','cerevisiae')
# Genomics (%GC, and chromosome number) ----------------------------------------
sc_gc =  get_ensembl_gc(ENSG = sc_ref$ensg, BIOMART = sc_biomart) %>%
         rename(ENS.GC_gene=percentage_gene_gc_content)

chromosomes = c(LETTERS[1:2],"2mu",LETTERS[3:16],"mito")
sc_chrom = sc_sgd %>%
  filter(!is.na(start) & !is.na(end) ) %>%
  mutate(chromosome_name = factor(chr, levels=unique(chr), labels = chromosomes)) %>%
  dplyr::select(sgdid,name,type,qual,chromosome_name) %>%
  mutate( chr_val=T ) %>%
  pivot_wider(id_cols=c('sgdid','name','type','qual'),
              names_from=chromosome_name, names_prefix='chr_',
              values_from = chr_val, values_fill = F )  %>%
  distinct()

sc_chr = sc_chrom %>% filter(type == 'ORF') %>% rename(ensg=name)

# Transcriptomics --------------------------------------------------------------
sc_transcript = get_ensembl_tx(ENSG=sc_ref$ensg, ENSP=sc_ref$ensp, BIOMART=sc_biomart)
id_cols = c("OS","sgdid","ORF","uniprot","AC","ID","GN","GNAME","NAME","SV","ensg","enst","ensp","ensp_canonical")
uni_col  = c("DB","OX","PE","one2one","CHR","has_cdna","id_cdna","is_uniref","has_ensp","has_ensg","cdna_len","prot_len")
# hgnc_col = c("symbol","name","location","gene_group","locus_group","locus_type")
sgd_col = c("qual","type",paste0('chr_',chromosomes),'has_orf')
ens_col  = c("canonical","is_canonical","is_enspref","has_ensp_canonical",
             "has_introns",'n_proteins','n_transcripts','n_exons','n_exons_mini',
             "cds_len","transcript_len","gene_len")

SC_CODING = left_join(sc_ref,sc_transcript, by=c('ensg'), suffix=c('','_canonical')) %>%
            left_join(sc_gc) %>%
            left_join(sc_chr) %>%
            relocate(all_of(id_cols),all_of(uni_col),all_of(sgd_col),all_of(ens_col)) %>%
            dplyr::rename_with(all_of(uni_col),.fn = Pxx, 'UP') %>%
            dplyr::rename_with(all_of(sgd_col),.fn = Pxx, 'SGD') %>%
            dplyr::rename_with(all_of(ens_col),.fn = Pxx, 'ENS') %>%
            arrange(desc(UP.is_uniref))


# Codons -----------------------------------------------------------------------
library(coRdon)
codon_table = get_codon_table()
sc_codon = bind_cols(orf=names(s288c_cdna),sc_codon_freq) %>%
  rename(all_of(set_names(codon_table$CODON,codon_table$codon_aa))) %>%
  rename_with(-orf, .fn=Pxx, 'SGD_cdna')

# Add amino acid with its associated codons
yeast_codon_usage.rds = here::here('output','sc-codon-usage.rds')
sc_CU = preload(yeast_codon_usage.rds, load.codon.usage(cds=s288c_cdna,with.counts=F,sp = 'sce') ) %>%
        dplyr::rename_with(starts_with('CU_'),.fn=Pxx,'elek2022')

# Single amino-acid frequencies ------------------------------------------------
sc_aa = sc_aa_freq %>% rename_with(.fn = Pxx, px="SGD_prot.f",s='_',.cols=-id)
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
AACLASS.FR = map(aa.prop, sum.aa.fr, BS=s288c_prot )
sc_aa_class = AACLASS.FR %>% purrr::reduce(full_join,by='id') %>%
              rename_with(~aa.class, starts_with('fr')) %>%
              rename_with(.fn = Pxx, px="SGD_prot.f",s='_',.cols=-id)

SC_COUNT = left_join(sc_aa,sc_aa_class) %>% left_join(sc_codon,by=c('id'='orf')) %>%
            dplyr::rename(ORF=id) %>% relocate(ORF)

# Conservation -----------------------------------------------------------------
 # sc_evo = load.evorate(alndir = "/media/WEXAC_data/FUNGI/fasta/",
 #              resdir="/media/WEXAC_data/FUNGI/",
 #              ref='Saccharomyces_cerevisiae',ext.r4s = 'raw.r4s')
sc_r4s = read_rds(here::here('data','evorate-fungi-orthologs.rds')) %>%
         group_by(id) %>% #= #read_rds(here("data","RESIDUE-EVORATE.rds")) %>% group_by(id) %>%
         summarize(lmsa=max(msa_pos), lref=max(ref_pos),
                   pid=mean(matched==total-1), qid=mean(mismatched>0), pgap=mean(indel>0),
                   nsp = mean(total),
                   rate = mean_(r4s_rate)) %>%
                   #rate_norm = abs(2*min_(r4s_rate)) + mean_(r4s_rate)) %>% # Only to use if the rates have negative values
             dplyr::filter(!is.na(rate)) %>%
             dplyr::rename_with(-id, .fn = Pxx, 'r4s_fungi')
#plot(density(x = sc_r4s$r4s_fungi.rate))

# Disorder ---------------------------------------------------------------------
sc_d2p2 =  read_rds(here("data","d2p2-yeast-uniprotKB.rds")) %>%
            get.d2p2.diso(as.df=T) %>%
            summarize.d2p2() %>%
            mutate(uniprot = str_extract(d2p2.id, UNIPROT.nomenclature())) %>%
            dplyr::select(-d2p2.id)
# Peptide stats ----------------------------------------------------------------
sc_pepstats = load.dubreuil2019.data(4) %>%
              dplyr::select(UNIPROT,PEPSTATS.netcharge=netcharge,PEPSTATS.mw=MW,PEPSTATS.pI=pI,UP.prot_len=prot.size) %>%
              group_by(UNIPROT) %>% mutate(PEPSTATS.mean_MW = PEPSTATS.mw/UP.prot_len) %>%
              dplyr::filter(UP.prot_len > 50) %>%
              mutate(PEPSTATS.AA_costly = PEPSTATS.mean_MW > 118, PEPSTATS.AA_cheap=PEPSTATS.mean_MW <= 105)

# Dubreuil et al 2019 (Disorder/Stickiness) ------------------------------------
IUP_cols = c("L.IUP20" = 'IUP20_L',"L.IUP30" = 'IUP30_L',"L.IUP40" = 'IUP40_L',
             "f.IUP20" = 'IUP20_f',"f.IUP30" = 'IUP30_f',"f.IUP40" = 'IUP40_f')
subset_cols = c('standard'='iupred_standard','medium'='iupred_medium','high'='iupred_high',
                'average'='uniprot_average','large'='uniprot_long','small'='uniprot_small',
                'isMB'='topcons_membrane')

iseq = 1:length(sc_prot)
full_score <- pbmcapply::pbmclapply(
  X=iseq,  FUN = function(x){ get_aa_score(string = as.character(sc_prot[[x]]) ) },
  mc.cores = 14,mc.cleanup = T)

sc_fullscore = tibble(id=names(sc_prot), length=widths(sc_prot)) %>%
                bind_cols(bind_rows(full_score)) %>%
                left_join(sc_uni2ens,by=c('id'='AC')) %>% # adding all known matched ORF/uniprot (even outside of reference proteome)
                group_by(id,ensp,length) %>%
                summarize( across(aggrescan:wimleywhite, ~ sum(.x)/length ) ) %>%
                distinct() %>%
                dplyr::rename_with(.cols = -c(id,ensp,length), .fn = xxS, 'full',s='_' )

sc_diso = load.dubreuil2019.data(4) %>%
          dplyr::select('UNIPROT','STRING', contains('IUP'), ends_with(c('dom','iup20','iup40')),
                        c('standard','medium','high','small','average','large','isMB'),
          ) %>%
          dplyr::rename(set_names(names(IUP_cols),IUP_cols),
                        set_names(names(subset_cols),subset_cols)) %>%
          dplyr::rename_with(.cols= ends_with(c('.dom','.iup20','.iup40')), .fn = str_replace_all, "\\.", "_" )

sc_dubreuil = left_join(sc_diso,sc_fullscore, by=c('UNIPROT'='id','STRING'='ensp')) %>%
              dplyr::rename_with(.cols=-c('UNIPROT','STRING'), .fn = Pxx, 'dubreuil2019', s='.')

# Domains ----------------------------------------------------------------------
sc_pfam=load.pfam(tax = TAXID) %>%
             filter(seq_id %in% sc_uniref) %>%
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

sc_pfam_count = left_join(sc_ref,sc_pfam, by=c('AC'='seq_id')) %>%
             group_by(AC) %>%
             summarize( pfam.HMM_none= is.na(pfam.ndom) | pfam.ndom==0,
                        pfam.HMM_single= !is.na(pfam.ndom) & pfam.ndom==1,
                        pfam.HMM_pair= !is.na(pfam.ndom) & pfam.ndom==2,
                        pfam.HMM_multi=!is.na(pfam.ndom) & pfam.ndom>=3) %>% distinct()

sc_pfam_dom = pivot_wider(sc_pfam %>% mutate(pfam_val=T), id_cols=seq_id,
                          names_from = 'hmm_name',names_prefix = 'pfam.dom_',
                          values_from = 'pfam_val', values_fill = F, values_fn = sum ) %>%
               mutate(across(where(is.integer), as.logical)) %>%
              dplyr::select(seq_id,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()

sc_pfam_clan = pivot_wider(sc_pfam %>% mutate(pfam_val=T), id_cols=seq_id,
                          names_from = 'clan_name',names_prefix = 'pfam.clan_',
                          values_from = 'pfam_val', values_fill = F, values_fn = sum ) %>%
              mutate(across(where(is.integer), as.logical)) %>%
              dplyr::select(seq_id,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()

SC_PFAM = left_join(sc_pfam_count,sc_pfam_dom,by=c('AC'='seq_id')) %>%
          left_join(sc_pfam_clan,by=c('AC'='seq_id')) %>%
          ungroup() %>%
          distinct()

sc_supfam=load.superfamily(tax = 'sc') %>%
          dplyr::rename(seqid = "sequence_id") %>%
          janitor::clean_names() %>%
          dplyr::filter(seqid %in% sc_ref$ensp) %>%
          group_by(seqid) %>% mutate(superfamily.ndom=n_distinct(superfamily_id)) %>%
          add_count(superfamily_id,name="pfam.repeat")

sc_supfam_count = left_join(sc_ref,sc_supfam, by=c('ensp'='seqid')) %>% group_by(ensp) %>%
                 summarize(superfamily.supfam_none= is.na(superfamily.ndom) | superfamily.ndom==0 ,
                 superfamily.supfam_single= !is.na(superfamily.ndom) & superfamily.ndom==1,
                 superfamily.supfam_pair= !is.na(superfamily.ndom) & superfamily.ndom==2,
                 superfamily.supfam_multi = !is.na(superfamily.ndom) & superfamily.ndom>=3)
sc_superfamilies = pivot_wider(sc_supfam %>% mutate(supfam_val=T), id_cols=seqid,
                          names_from = 'superfamily_description',names_prefix = 'superfamily.SF_',
                          values_from = 'supfam_val', values_fill = F, values_fn = sum ) %>%
                  mutate(across(where(is.integer), as.logical)) %>%
                  dplyr::select(seqid,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()

sc_families = pivot_wider(sc_supfam %>% mutate(supfam_val=T), id_cols=seqid,
                               names_from = 'family_description',names_prefix = 'superfamily.F_',
                               values_from = 'supfam_val', values_fill = F, values_fn = sum ) %>%
              mutate(across(where(is.integer), as.logical)) %>%
              dplyr::select(seqid,where(~ is.logical(.x) && sum(.x) >9 ))%>% distinct()

SC_SUPFAM = left_join(sc_supfam_count,sc_superfamilies,by=c('ensp'='seqid')) %>%
            left_join(sc_families,by=c('ensp'='seqid')) %>%
            ungroup() %>%
            distinct()

# Folding energy and stability -------------------------------------------------
sc_stab = load.leuenberger2017.data("S. cerevisiae",rawdata = F) %>%
          add_count(protein_id,name='npep') %>%
          distinct() %>%
          mutate( nres = round(0.01*protein_coverage*length),
                  Tm_stable = protinfo=="Stable",
                  Tm_medium = protinfo=="Medium",
                  Tm_unstable = protinfo=="Unstable") %>%
          dplyr::select(-c(protinfo,protein_coverage,length,
                           measured_domains,theoretical_number_of_domain,nres)) %>%
          rename_with(.fn=Pxx, px='leuenberger2017.LIP',s='_',.cols=-protein_id)

sc_tm = load.jarzab2020.data(org = "S.cerevisiae") %>%
        dplyr::select(UNIPROT,GENENAME,Tm_celsius,Tm_type,AUC) %>%
         group_by(UNIPROT,GENENAME) %>%
         mutate(Tm_mean=mean_(Tm_celsius),Tm_max=max_(Tm_celsius))

tm_low_celsius = max_(sc_tm$Tm_celsius[sc_tm$Tm_type == 'Low-Tm'])
tm_med_celsius = max_(sc_tm$Tm_celsius[sc_tm$Tm_type == 'Medium-Tm'])
tm_nomelt_celsius =max_(sc_tm$Tm_celsius[sc_tm$Tm_type =='Non-melter' ])

sc_tm = sc_tm %>% mutate(Tm_low = Tm_max < tm_low_celsius,
                         Tm_med = between(Tm_max,tm_low_celsius,tm_med_celsius),
                         Tm_high = Tm_max > tm_med_celsius,
                         Tm_nomelt = is.na(Tm_max)) %>%
        dplyr::select(-c(Tm_type,AUC,Tm_celsius)) %>% distinct() %>%
        rename_with(.fn=Pxx, px='jarzab2020.TPP',s='_',.cols=-c(UNIPROT,GENENAME))

SC_FOLD = left_join(sc_tm,sc_stab, by=c('UNIPROT'='protein_id'))
# Complexes --------------------------------------------------------------------
nmers= paste0(c('mono','di','tri','tetra','penta','hexa','septa','octa','nona','deca'),"mer")
sc_complex = load.meldal.2019.data(species = 'yeast') %>%
  filter(is_uniprot) %>% # BASED ON UNIPROT
  mutate(oligomers = cut(n_members, breaks = c(1:10,20,81),
                         labels = paste0("meldal2019.CPX_",c(nmers[2:10],'high_oligomer','molecular_machine'))),
         CPLX_ASSEMBLY = gsub("(.+)(\\.$)","\\1",CPLX_ASSEMBLY),
         CPLX_ASSEMBLY = str_replace_all(CPLX_ASSEMBLY," ","_"),
         CPLX_NAME = str_replace_all(CPLX_NAME," ","_")) %>%
  mutate(oligomers_val=T)

sc_oligomers = pivot_wider(sc_complex, id_cols = members, names_from = 'oligomers',
                           names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum,values_fill = F) %>%
  mutate(across(where(is.integer), as.logical)) %>%
  dplyr::select(uniprot=members,where(~ is.logical(.x) && sum(.x) > 9 )) %>%
  dplyr::rename(meldal2019.unknown_oligomer = meldal2019.NA)

sc_assembly = pivot_wider(sc_complex, id_cols = members,
                          names_from = 'CPLX_ASSEMBLY', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum,values_fill = F) %>%
  mutate(across(where(is.integer), as.logical)) %>%
  dplyr::select(uniprot=members,where(~ is.logical(.x) && sum(.x) > 9 )) %>%
  dplyr::rename(meldal2019.unknown_assembly = meldal2019.NA)

sc_complexes = pivot_wider(sc_complex, id_cols = members,
                           names_from = 'CPLX_NAME', names_prefix = 'meldal2019.',
                           values_from = oligomers_val, values_fn=sum, values_fill = F) %>%
                mutate(across(where(is.integer), as.logical)) %>%
                dplyr::select(uniprot=members,where(~ is.logical(.x) && sum(.x) > 4 ))
SC_COMPLEX = left_join(sc_oligomers,sc_assembly) %>% left_join(sc_complexes)

# Functional interactions ------------------------------------------------------
sc_string = load.string(tax="4932",phy=F, ful=T, min.score = 900) %>%
  mutate(ens1 = str_extract(protein1,SGD.nomenclature()),
         ens2 = str_extract(protein2,SGD.nomenclature())
  ) %>% relocate(ens1,ens2) %>% dplyr::select(-c(protein1,protein2))

yeast_centrality.rds = here::here('output','sc-string-centralities.rds')
sc_string_centralities = preload(yeast_centrality.rds, network.centrality(sc_string %>% dplyr::select(ens1,ens2)))

sc_string_centralities$megahub_func = sc_string_centralities$cent_deg >= 300
sc_string_centralities$superhub_func = between(sc_string_centralities$cent_deg,100,300)
sc_string_centralities = sc_string_centralities %>% type_convert() %>%
  dplyr::rename_with(.cols = -ids, .fn = str_replace_all, pattern='string.',replacement="") %>%
  dplyr::rename_with(-ids,.fn=Pxx, 'STRING')

# Biological functions (GO) ----------------------------------------------------
unigo.rds = here::here('output','sc-uniprot-go.rds')

unigo = preload(unigo.rds, get.uniprot.go(sc_uniref,taxon = 4932) )
sc_unigo = unigo %>%  # Uniprot-based GO annotation
           mutate( go = paste0(ONTOLOGY,"_", str_replace_all(goterm,pattern="[^A-Za-z0-9\\.]+","_"))) %>%
           filter(!obsolete & shared > 10 & ONTOLOGY != "CC") %>%
           arrange(go)

# Convert to a matrix format (columns = GO term, rows=proteome)
SC_UNIGO  = sc_unigo %>%  mutate(seen=T) %>%
           pivot_wider( id_cols = c('UNIPROT'), names_from = 'go', values_from = 'seen',
               values_fn=list(seen = unique),values_fill = list(seen=F)) %>%
           dplyr::rename_with(.cols = -UNIPROT,.fn = Pxx, 'go.', s='')

# Subcellular locations (Uniprot) ----------------------------------------------
sc_uniloc = query_uniprot_subloc(uniprot = sc_uniref,todf=T) # as a wide dataframe
SC_UNILOC = sc_uniloc %>%
            group_by(id) %>%
            dplyr::select(id, where(~is.numeric(.x) && sum(.x) >50 ) ) %>%
            dplyr::rename_with(-id,.fn=Pxx, 'UP.loc', s='_') %>%
            mutate( across(.cols=everything(),.fns= as.logical))

# Biological pathways (KEGG) ---------------------------------------------------
yeast_pathways.rds = here('output','sc-kegg-pathways.rds')
sc_pathways=preload(yeast_pathways.rds,
                    get.KEGG(sp='sce',type='pathway',as.df=T,to_uniprot = T)) %>%
            mutate(desc=str_replace_all(tolower(desc),' ','_'))

kegg_pathways = sc_pathways %>% mutate(pathway_val=T) %>%
                pivot_wider(id_cols=id,names_from = 'desc',names_prefix = 'KEGG.path_',
                     values_from = 'pathway_val', values_fill = F, values_fn = sum ) %>%
                mutate(across(where(is.integer), as.logical)) %>%
                dplyr::select(id,where(~ is.logical(.x) && sum(.x) >9 ))

yeast_modules.rds = here('output','sc-kegg-modules.rds')
sc_modules=preload(yeast_modules.rds,
                   get.KEGG(sp='sce',type='module',as.df=T,to_uniprot = T)) %>%
  mutate(desc=str_replace_all(tolower(desc),' ','_'))

kegg_modules = sc_modules %>% mutate(module_val=T) %>%
  pivot_wider(id_cols=id,names_from = 'desc',names_prefix = 'KEGG.mod_',
              values_from = 'module_val', values_fill = F, values_fn = sum ) %>%
          mutate(across(where(is.integer), as.logical)) %>%
          dplyr::select(id,where(~ is.logical(.x) && sum(.x) >9 ))

SC_KEGG = left_join(kegg_pathways,kegg_modules)

# Integrate all datasets -------------------------------------------------------

# Based on uniprot accession
head(sc_r4s) # id = uniprot AC
sc_orthologs = unique(sc_r4s$id)
dim(SC_CODING)
colnames(SC_CODING)
dim(SC_COUNT)
dim(sc_CU)
dim(sc_d2p2)
dim(sc_dubreuil)
dim(SC_PFAM)
dim(SC_SUPFAM)
dim(sc_pepstats)
dim(SC_UNIGO)
dim(sc_string_centralities)
dim(SC_UNILOC)
dim(SC_FOLD)
dim(SC_KEGG)
dim(SC_COMPLEX)

SC_DATA = left_join(SC_CODING,SC_COUNT,by=c('ORF')) %>%
  left_join(sc_CU,by=c('ORF'='ID')) %>%
  left_join(sc_d2p2,by=c('uniprot')) %>%
  left_join(sc_dubreuil,by=c('uniprot'='UNIPROT')) %>%
  left_join(SC_PFAM,by=c('uniprot'='AC')) %>%
  left_join(SC_SUPFAM,by=c('ensp'='ensp')) %>%
  left_join(sc_pepstats,by=c('uniprot'='UNIPROT','UP.prot_len')) %>%
  left_join(SC_UNIGO,by=c('uniprot'='UNIPROT')) %>%
  left_join(sc_string_centralities,by=c('ORF'='ids')) %>%
  left_join(SC_UNILOC,by=c('uniprot'='id')) %>%
  left_join(SC_FOLD,by=c('uniprot'='UNIPROT')) %>%
  left_join(SC_KEGG,by=c('uniprot'='id')) %>%
  left_join(SC_COMPLEX,by=c('uniprot')) %>%
  left_join(sc_r4s,by=c('ORF'='id')) %>%
  ungroup %>%
  distinct() #%>%
  #relocate(uniprot,UP.is_uniref,ensg,ensp,gname, GENENAME, gene_biotype)

# Replace missing values  ------------------------------------------------------
miss0 = check_missing_var(SC_DATA)


#### 1. Fix logical variables (NA replaced by FALSE) ####
SC_FEATURES.1 = SC_DATA %>%
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
  mutate( across( where(is.logical) & starts_with('SGD.chr'), ~replace_na(., F)) ) %>%
  distinct()

miss1 = check_missing_var(SC_FEATURES.1)

#### 2. Remove rare variables  (less than 2 occurrences in proteome) ####
SC_FEATURES.2 = remove_rare_vars(df=SC_FEATURES.1,min_obs=2)
miss2 = check_missing_var(SC_FEATURES.2)

#### 3. Fix network centrality (lower confidence + random forest) ####
# Take a long time (@! ~50mn for human !@)
SC_FEATURES.3 = fix_missing_centrality(df=SC_FEATURES.2, id='ORF', col_prefix="STRING.", taxon=4932)
SC_FEATURES.3 = SC_FEATURES.3 %>% mutate(across(starts_with("STRING.cent_"), .fns = min_))

miss3 = check_missing_var(SC_FEATURES.3)

#### 4. Fix thermodynamcis stability (replace NA by average) ####
SC_FEATURES.4 = SC_FEATURES.3 %>%
                 mutate(across(starts_with(c('jarzab2020','leuenberger2017')),
                               ~coalesce(as.numeric(.x), mean_(get(cur_column())))))

miss4 = check_missing_var(SC_FEATURES.4)

#### 5. Fix PEPSTATS missing values (i.e not included in Dubreuil et al. 2019) ####
SC_FEATURES.5 = fix_missing_peptide_stats(df=SC_FEATURES.4, id='ORF', taxon=4932,
                                         col_len='UP.prot_len',
                                         col_mw='PEPSTATS.mw',
                                         col_mw_avg='PEPSTATS.mean_MW',
                                         col_charge='PEPSTATS.netcharge',
                                         col_pi='PEPSTATS.pI')
miss5 = check_missing_var(SC_FEATURES.5)

#### 6. Fix D2P2 missing values (fetch from d2p2 or use consensus prediction from mobiDB) ####

# Get all uniprot or ensembl reference identifiers
uni_na_d2p2 = SC_FEATURES.5 %>% filter(is.na(D2P2.diso_len)) %>%  pull(uniprot)
ensp_na_d2p2 = SC_FEATURES.5 %>% filter(is.na(D2P2.diso_len) & !is.na(ensp)) %>% pull(ensp)
#n_distinct(uni_na_d2p2)
#n_distinct(ensp_na_d2p2)
id_na_d2p2 = c(uni_na_d2p2,ensp_na_d2p2)

#  Fetch d2p2 missing identifiers
d2p2_na =  load.d2p2(id_na_d2p2, 'output/sc-missing-d2p2.rds') %>%
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
SC_FEATURES.6 = coalesce_join(SC_FEATURES.5,uni_d2p2_nona,by=c('uniprot')) %>%
                coalesce_join(ensp_d2p2_nona,by=c('ensp')) %>%
                mutate(has_d2p2 = ifelse(is.na(has_d2p2),T,has_d2p2) )

miss6 = check_missing_var(SC_FEATURES.6)

#save.image(here('output','checkpoint-6-hs-data.rdata'))
#load(here::here('output','checkpoint-6-hs-data.rdata'))
#### 7. Fix missing amino acid propensities ####
id_na_aascore = SC_FEATURES.6 %>%
                    dplyr::select(AC,starts_with('dubreuil2019')) %>%
                    dplyr::filter(row_number() %in% find_na_rows(.,as.indices = T)) %>%
                    pull(AC) %>% unique  %>% intersect(names(sc_prot))
seq_aa_score = sc_prot[id_na_aascore]
na_aascore = retrieve_missing_aascore(id_na_aascore,seq_aa_score)
df_na_aascore = na_aascore %>%
                dplyr::rename_with(.cols=-acc, .fn=str_replace_all, 'pawar_', 'pawar_ph7_' ) %>%
                dplyr::select(-contains(c('voronoi'))) %>%
                dplyr::rename_with(.cols=-acc, .fn=Pxx, 'dubreuil2019', s='.' ) %>%
                dplyr::rename(uniprot='acc')


iup40_cols = str_subset(colnames(SC_FEATURES.6),'^dubreuil.+iup40$') %>% sort
dom_cols = str_subset(colnames(SC_FEATURES.6),'^dubreuil.+dom$') %>% sort
full_cols = str_subset(colnames(SC_FEATURES.6),'^dubreuil.+full$') %>% sort %>% setdiff("dubreuil2019.voronoi_sickiness_full")

SC_FEATURES.7 = coalesce_join(SC_FEATURES.6,df_na_aascore,by='uniprot') %>%
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

miss7 = check_missing_var(SC_FEATURES.7)

#### 8. Fix missing cDNA sequences ####

#cdna_na = HS_FEATURES.7 %>% dplyr::select( all_of(miss7$variable) )
#missing_cdna = HS_FEATURES.7 %>% filter( is.na(UP_cdna.AAA_Lys_K) ) %>% pull(uniprot)
#get.uniprot.proteome(9606,DNA=T)
### 18.09.2022 ==> CURRENT CHECKPOINT (I STOPPED HERE WHEN TESTING THE SCRIPT)

SC_FEATURES.8 = SC_FEATURES.7 %>%
                mutate( ENS.cds_len = coalesce( ENS.cds_len, UP.cdna_len),
                        ensp = coalesce(ensp, ensp_canonical),
                        GENENAME = coalesce(GENENAME),GN) %>%
                dplyr::select(-contains('ratioR_RK'))
miss8 = check_missing_var(SC_FEATURES.8 %>% ungroup)

SC_PROTEOME_DATA = SC_FEATURES.8 %>%
                     # Get out the rows with NA
                     drop_na(starts_with(c('UP_cdna.','elek2022.','ENS.','dubreuil2019.'))) %>%
                     dplyr::select(-GENENAME) %>%
                     # RENAME VARIABLES NON-ALPHANUMERIC CHARACTERS (NOT VALID FOR FORMULA)
                     dplyr::rename_with(.cols = everything(), .fn = str_replace_all, pattern="[^A-Za-z0-9\\.]+", replacement="_")

miss = check_missing_var(SC_PROTEOME_DATA)
dim(SC_PROTEOME_DATA)

saveRDS(SC_PROTEOME_DATA,here::here('released-dataset','yeastOmics-v1.rds'))

# Retrieve yeast abundance -----------------------------------------------------
sc_uni = sc_map_uni %>% filter(extdb=='UniProtKB-ID') %>% dplyr::select(uni,extid)
# raw paxdb datasets + ho et al. 2018
ho_et_al_2018 = load.abundance() %>%
                set_names(c("orf","ho2018_MPC","ho2018_MDPC","ho2018_gfp","ho2018_ms")) %>%
                ungroup()

sc_ho_dataset = skimr::skim(ho_et_al_2018[,-1]) %>% dplyr::select(skim_variable,complete_rate) %>%
                as_tibble() %>% dplyr::rename(id=skim_variable,cov=complete_rate) %>%
                mutate(organ='WHOLE_ORGANISM',taxon="4932", filename=paste0(taxon,"-",id),
                       is_integrated=T, score=NA, w=NA, yr="2018",ndata=19, cov=as.character(100*round(cov,1)))

sc_paxdb_datasets = find_paxdb_datasets(4932) %>% bind_rows(sc_ho_dataset)


sc_ppm = load.paxdb(4932,rm.zero=T) %>%
            pivot_wider(id_cols = c(protid,id_uniprot), names_from = c('id'),
                       values_from=ppm, values_fn = mean_) %>%
            left_join(sc_uni, by=c('id_uniprot'='extid')) %>%
            left_join(sc_r4s,by=c('protid'='id')) %>%
            mutate(log10.rate = log10(r4s_fungi.rate)) %>%
            left_join(ho_et_al_2018,by = c('protid'='orf'))


# integrated paxdb datasets
sc_paxdb = get.paxdb(tax=4932,abundance = 'integrated',rm.zero=T)
SC_PPM = sc_paxdb %>%
          group_by(protid,id_uniprot) %>%
          summarize( PPM_MIN_ORGAN = min_(ppm_int),
                   PPM_MAX_ORGAN = max_(ppm_int),
                   PPM_AVG_ORGAN = mean_(ppm_int),
                   PPM_MD_ORGAN = median_(ppm_int),
                   PPM_geomAVG_ORGAN = geomean(ppm_int)) %>%
          mutate( across(starts_with('PPM_'), log10) ) %>%
          left_join(sc_uni, by=c('id_uniprot'='extid')) %>%
          left_join( pivot_wider(sc_paxdb, id_cols = c(protid,id_uniprot,n_data,n_int),
                           names_from = 'organ', names_prefix = 'PPM_',
                           values_from = 'ppm_int', values_fn = log10) )


save(list = ls(pattern = '^sc_'), file = here('output','sc_datasets.rdata'))
save(list = ls(pattern = '^SC_'), file = here('output','sc_integrated_datasets.rdata'))

# Orthologs dataset (with/out abundance) ---------------------------------------
# Full proteome N=6549 - after removing NA = 5729
n_distinct(SC_DATA$ORF)
n_distinct(SC_PROTEOME_DATA$ORF)

SC_ORTHOLOGS = SC_PROTEOME_DATA %>%
                filter(!is.na(r4s_fungi.rate) & ORF %in% sc_r4s$id) %>%
                filter(!is.na(uniprot) & !is.dup(uniprot)) %>%
                filter(!is.na(ensp) & !is.dup(ensp)) %>%
                filter(!is.na(ORF) & !is.dup(ORF)) %>%
                # Remove rare variables among orthlogs (min 3 orthologs must share the feature)
                remove_rare_vars(df=.,min_obs=3) %>%
                distinct()

dim(SC_ORTHOLOGS)
n_distinct(SC_ORTHOLOGS$ORF)
# All orthologs N = 3726

n_distinct(SC_PPM$protid)

# Define two predictions datasets with/out abundance (training/validation)
MPC_nona = ho_et_al_2018 %>% drop_na(ho2018_MPC)
sc_validation = SC_ORTHOLOGS %>% filter( !(ORF %in% MPC_nona$orf) )
sc_orthologs =  SC_ORTHOLOGS %>% filter(ORF %in% MPC_nona$orf)
dim(sc_validation) # n=86
dim(sc_orthologs) # n=3639

# Evo Rate vs. Expression -------------------------------------------------
all_ppm_er_cor = map_dfr(sc_paxdb_datasets$id ,
                         function(x){ spearman.toplot(X=sc_ppm$log10.rate, Y=sc_ppm[[x]]) } ) %>%
  add_column(paxdb = sc_paxdb_datasets$id) %>%
  arrange(desc(abs(estimate))) %>%
  left_join(sc_paxdb_datasets, by=c('paxdb'='id')) %>%
  dplyr::select(organ,is_integrated,estimate,N,p.value,filename,cov,yr) %>%
  mutate(filename = ifelse(is_integrated,
                           str_remove_all(filename,"(9606\\-|\\-integrated\\.txt)"),
                           str_remove_all(filename,"(9606\\-|\\.txt|_Maxquant|_gene|SEQUEST)")),
         cor_range = cut(estimate,breaks=c(-0.7,-0.5,-0.3,-0.2,-0.1,0.1)),
         filename_n = paste0(filename," (N ",N,")")) %>%
  print(n=200)

p_ER_int = ggplot(all_ppm_er_cor,
                  aes(y=reorder(filename_n,abs(estimate)),fill=organ,x=estimate)) +
  geom_col() +
  #geom_text(aes(label=N,x=1.2*max(estimate))) +
  geom_text(aes(label=round(estimate,3),x=estimate),size=3,hjust='inward') +
  theme(legend.position='none')

ggsave(p_ER_int, filename=here::here('plots','sc-abundance-evolution.pdf'))

# Explaining variance in evolutionary rate -------------------------------------
#IDCOLS = c("uniprot","ensp","ensg","GENENAME","gene_biotype")
ID_COLS = c('ORF','GN')
XCOL="ho2018_MPC"
YCOL="ER"
ZCOL=NULL

sc_ortho_predictors = sc_orthologs %>% group_by(ORF,GN,uniprot) %>%
  dplyr::select(where(is.numeric) | where(is.logical)) %>%
  dplyr::select(-c('SV', starts_with("r4s_fungi"),"has_unique_id") )

SC_PREDICTORS = sc_ortho_predictors %>% ungroup

SC_PPM$PPM = SC_PPM$PPM_WHOLE_ORGANISM
SC_ER = left_join(sc_r4s,SC_PPM, by=c('id'='protid')) %>%
        left_join(ho_et_al_2018,by=c('id'='orf')) %>%
        mutate(ER = log10(r4s_fungi.rate))
        #drop_na(ho2018_MPC,ER)


ggplot(SC_ER, aes(x=ho2018_MPC,y=ER)) + geom_point() + geom_smooth()
spearman.toplot(SC_ER$ho2018_MPC,SC_ER$ER)

SC_LMDATA = left_join(SC_PREDICTORS,SC_ER,by=c('ORF'='id')) %>%
            drop_na(ho2018_MPC) %>%
            filter(!is.na(ORF) & !is.dup(ORF)) %>%
            distinct()

fit0 = fit_m0(INPUT_LM = SC_LMDATA,PREDICTORS = SC_PREDICTORS,
              XCOL = XCOL,YCOL = YCOL, ZCOL=NULL, IDCOLS=ID_COLS,
              MAX_XCOR = 0.5, MAX_YCOR=0.5, MIN_N = 2)

df0=decompose_variance(fit0$LM0,T)
spearman.toplot(fit0$P$ER,fit0$P$ho2018_MPC)

#save.image(here('output','checkpoint-sc-orthologs-data.rdata'))
load(here::here('output','checkpoint-sc-orthologs-data.rdata'))

# Filter Dataset for prediction --------------------------------------------------
sc_evo_all= preload(here("output","lm-evo-sc.rds"), select_variable(fit0,response='.resid', raw=T))
sc_best= sc_evo_all %>% filter(pc_ess > 0.5 & variable != YCOL)
sc_best= sc_evo_all %>% filter(pc_ess > 1 & variable != YCOL)
sc_best= sc_evo_all %>% filter(pc_ess > 2 & variable != YCOL)
sc_best= sc_evo_all %>% filter(pc_ess > 5 & variable != YCOL)

sc_best_pred = fit0$P[,c(XCOL, YCOL, '.resid', ZCOL, intersect(sc_best$variable, colnames(fit0$P)))]
sc_nbest = n_distinct(sc_best$variable)

formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM_ER = lm(data=sc_best_pred, formula_null)
decompose_variance(LM_ER)

#formula_sc_best_ER=paste0(YCOL," ~ ",paste0( P_evo$variable,collapse=" + "))
formula_sc_best = reformulate(sc_best$variable, response =  '.resid')

m_best_ER = step(object=LM_ER, scope = as.formula(formula_sc_best), direction = 'forward',k=log(sc_nbest)*2,trace=0)
decompose_variance(m_best_ER,to.df = T)

sc_best_lm = lm(reformulate(response = YCOL, termlabels = c(XCOL,labels(m_best_ER)),data=sc_best_pred)
print(labels(sc_best_lm))
print(sc_nbest)
decompose_variance(sc_best_lm,to.df = T)

formula_m0 = reformulate(response=YCOL,termlabels = XCOL,intercept = T)
LM0 = lm(data=sc_best_pred, formula_m0)
decompose_variance(LM0,T)

## Validate
sc_valid = left_join(sc_validation,SC_ER,by=c('ORF'='id')) %>%
  filter(!is.na(ORF) & !is.dup(ORF) ) %>%
  distinct()

LM_validation_ER = lm(data=sc_valid, formula_null)
decompose_variance(LM_validation_ER,to.df = T)

sc_validation_lm = lm(reformulate(response = YCOL, termlabels = labels(sc_best_lm)),data=sc_valid)
print(labels(sc_validation_lm))
print(coefficients(sc_validation_lm)[is.na(coefficients(sc_validation_lm))])
decompose_variance(sc_validation_lm,to.df = T)

###
###
###
# #sc_evo_var = readRDS(here::here('output','sc_features_ess_over_0.2.rds'))
# sc_evo_var = readRDS(here::here('output','sc_features_ess_over_2.rds'))
#
# library(stringdist)
# sim_var = stringsimmatrix(sc_evo_var,colnames(HS_LMDATA), method='jw',p=1e-5,useNames='strings')
# maxS= apply(sim_var,1,max_)
# i_maxS= apply(sim_var,1,which.max)
# matched = tibble(sc_var=rownames(sim_var), hs_var = colnames(sim_var)[i_maxS], similarity=maxS)
#
# View(matched)
# write_delim(matched,here::here('output','sc-hs-features.tsv'),delim = '\t')
#
# ### YEAST TO HUMAN
# sc2hs = readr::read_delim('/home/benjamin/Downloads/sc-hs-features - Sheet1.tsv',na = 'NA')
# sc2hs$hs_var[sc2hs$sc_var=='sgd.pGC'] = 'UP.GC_cdna'
# sc2hs$hs_var[str_detect(sc2hs$sc_var,'pu2008')] = NA
# nomatch = sc2hs %>% filter( similarity < 1 & !is.na(hs_var) & hs_var == "" )
#
#
# sc_matched = sc2hs %>% filter(similarity == 1 | !is.na(hs_var) & hs_var != "" & hs_var %in% hs_evo_all$variable)
# nvar = n_distinct(sc_matched$hs_var)
# formula_sc2hs_best = reformulate(sc_matched$hs_var, response =  '.resid')
#
# m_best_ER = step(object=LM_ER, scope = as.formula(formula_sc2hs_best), direction = 'forward',k=log(nvar)*2,trace=0)
# decompose_variance(m_best_ER,to.df = T)
# sc_best = labels(m_best_ER)
#
# hs_best_lm = lm(reformulate(response = YCOL, termlabels = labels(m_best_ER)),data=hs_best_pred)
# print(intersect(labels(hs_best_lm),sc_best))
# print(hs_nbest)
# decompose_variance(hs_best_lm,to.df = T)
library(broom)

# STEPWISE
sc_perf = aov_model(sc_best_lm,min_pv = 1e-3)
nvar = nrow(sc_perf)
ess_pc = sum(sc_perf$ess)
p=aov_plot(aov_model(sc_best_lm,min_pv = 1e-3), name='best')


## ADD SIGN
## V1 = LM(ER ~ PPM + BEST_VAR)
## V2 = LM(ER ~ BEST_VAR)

## PRED EVO RATE WITH BEST VAR
cor(sc_best_pred$ER, predict(sc_best_lm))

ggsave(p, path = here::here('plots'),
       scale=1.1, width = 12,height=12, bg = 'white',
       filename = 'yeast-stepwise-best-model-ess_65pc_27var.pdf')


### HUMAN
