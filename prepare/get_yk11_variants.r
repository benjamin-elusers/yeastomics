#install.packages('rjson')
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("src","function_YK11.r"))

load_aligned_yk11 = function(){

  DIR_YK11_PAL2NAL = "/media/elusers/users/benjamin/A-PROJECTS/03_PostDoc/EvoRate-paper/pal2nal/"
  cds_yk11.rds=file.path(DIR_YK11_PAL2NAL,"cds_yk11.rds")
  prot_yk11.rds=file.path(DIR_YK11_PAL2NAL,"prot_yk11.rds")

  if(!file.exists(cds_yk11.rds) & !file.exists(prot_yk11.rds)){
    DIR_YK11_CDS = file.path(DIR_YK11_PAL2NAL,"output/")
    muscle_yk11_cds = list.files(DIR_YK11_CDS,pattern='\\.out',full.names = T)
    is_aligned = map(muscle_yk11_cds,~readDNAStringSet(.x) %>% widths %>% unique)
    not_aligned = sapply(is_aligned,length) > 1

    DIR_YK11_PROT = file.path(DIR_YK11_PAL2NAL,"aln_s288c/")
    muscle_yk11_prot = list.files(DIR_YK11_PROT,pattern='\\.fasta',full.names = T)

    orf_both = intersect(get.orf(muscle_yk11_prot),get.orf(muscle_yk11_cds[!not_aligned]))

    cds_aligned = file.path(DIR_YK11_CDS,paste0(orf_both,'_s288c.out'))
    prot_aligned = file.path(DIR_YK11_PROT,paste0(orf_both,'_strains_only.fasta'))

    cds_yk11 = map(cds_aligned,readDNAMultipleAlignment)
    names(cds_yk11) = orf_both
    prot_yk11 = map(prot_aligned,readAAMultipleAlignment)
    names(prot_yk11) = orf_both

    saveRDS(cds_yk11,cds_yk11.rds)
    saveRDS(prot_yk11,prot_yk11.rds)
  }

  YK11_ALIGNED = list(CDS = readRDS(cds_yk11.rds), PROT = readRDS(prot_yk11.rds))
  return(YK11_ALIGNED)
}

count_msa = function(MSA_TABLE, alpabet="[ARNDCQEGHILKMFPSTWYV]"){

  idref = unique(MSA_TABLE$id_4932)
  CLADE = unique(MSA_TABLE$clade_name)
  OG = unique(MSA_TABLE$OG)
  msa_aas = MSA_TABLE %>% dplyr::select(-clade_name, -OG, -id_4932, -POS, -{{idref}})
  MSA_TABLE$msa_col = apply(msa_aas,1,paste,collapse="")

  MSA_TABLE$ref_aa = MSA_TABLE[[idref]]
  MSA_TABLE$ref_gap = MSA_TABLE$ref_aa == "-"

  MSA_TABLE$p_gap = str_count(MSA_TABLE$msa_col,"-") + MSA_TABLE$ref_gap
  MSA_TABLE$p_defined = str_count(MSA_TABLE$msa_col,alpabet) + !MSA_TABLE$ref_gap

  MSA_TABLE$p_matched = rowSums(!MSA_TABLE$ref_gap & MSA_TABLE$ref_aa == msa_aas) + !MSA_TABLE$ref_gap
  MSA_TABLE$p_mismatched = rowSums(!MSA_TABLE$ref_gap & MSA_TABLE$ref_aa != msa_aas & msa_aas != "-")

  MSA_TABLE$nseq = ncol(msa_aas) + 1
  MSA_TABLE$len_msa = max_(MSA_TABLE$POS)
  MSA_TABLE$len_ref = sum_(!MSA_TABLE$ref_gap)

  msa_pos = tibble(
    CLADE = CLADE,
    OG = OG,
    ref_id = idref,
    pos_msa=MSA_TABLE$POS,
    len_msa=max_(MSA_TABLE$POS),
    len_ref=sum_(!MSA_TABLE$ref_gap),
    ref_aa=MSA_TABLE$ref_aa,
    ref_gap=MSA_TABLE$ref_gap,
    msa_col=MSA_TABLE$msa_col,
    p_defined=MSA_TABLE$p_defined,
    p_mismatched=MSA_TABLE$p_mismatched,
    p_matched=MSA_TABLE$p_matched,
    p_gap=MSA_TABLE$p_gap)

  return(msa_pos)
}


# Find Single Amino Acid Polymorphisms (SNP_AA) ----------------------------------
S288C = load.sgd.proteome(withORF=T,rm.stop=F) # Reference SGD protein sequences
YK11_SNP_AA.rds = here("data",'YK11-SNP_AA.rds')

# GET PROTEOME AMINO ACID POLYMORPHISMS
if( file.exists(YK11_SNP_AA.rds) ){
  PROT_SNP_AA=read_rds(YK11_SNP_AA.rds)
}else{
  # PROTEOME SEQUENCES
  WT = get.positions(S288C) %>%
    group_by(orf) %>%
    mutate(bin.pos=dplyr::ntile(wt_pos,100))

  # PROTEOME 1011 STRAINS SEQUENCES
  YK11_fastadir = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"
  YK11_prot = here("data","YK11-PROT.rds")
  YK11 = preload(YK11_prot,{load.1011.strains(seqdir = YK11_fastadir,seqtype='AA')},"loading proteome from 1011 strains...")

  # PROTEOME DATAFRAME
  P = tibble( orf=names(YK11),
              n_strains=lengths(YK11),
              len = widths(YK11)) %>%
    left_join(get.width(S288C), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
    mutate( match_wt = len == len.s288c )

  # FIND AMINO ACID VARIANTS
  YK11_aa_var = here("data","YK11-VAR_AA.rds")
  AA_VAR = preload(YK11_aa_var,{lapply(YK11, get.variants,verbose=F) %>% bind_rows()},"find amino acid variants...")

  # GET AMINO ACID POLYMORPHISMS
  SNP_AA = get.variants.to.SNP(AA_VAR) %>%
    add_count(id,ref_pos,name='nvar') %>%
    group_by(id,ref_pos) %>%
    mutate( alt_cumfr=sum(alt_fr)) %>%
    left_join(P,by=c('id'='orf'))
  PROT_SNP_AA = left_join(SNP_AA,WT, by=c('id'='orf','ref_pos'='wt_pos','len.s288c'='len')) %>%
    mutate(wt_low = alt_aa == wt_aa,
           wt_missing=is.na(wt_aa),
           early_stop = (alt_aa == "*" & ref_pos != len),
           dSTI.ref=get.score.mutation(ref_aa,alt_aa),
           dSTI.wt=get.score.mutation(wt_aa,alt_aa))
  # Save amino acid polymorphism
  saveRDS(PROT_SNP_AA,here("data",'YK11-SNP_AA.rds'))
}

# Counting amino acid polymorphisms per protein
snp_aa_per_orf = PROT_SNP_AA %>%
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(ref_pos))
snp_aa_per_orf %>% arrange(n_snp)
saveRDS(snp_aa_per_orf,here("data",'YK11-ORF-VAR_AA.rds'))

# Find Single Nucleotide Polymorphisms (SNP) -----------------------------------
CDS = load.sgd.CDS(withORF=T) # Reference SGD DNA coding sequences
cds_snp_nt.rds = here("data",'YK11-SNP_NT.rds')
# GET CODING GENOME NUCLEOTIDE POLYMORPHISMS
if( file.exists(cds_snp_nt.rds) ){
  cds_snp_nt=read_rds(cds_snp_nt.rds)
}else{
  # CDS GENOME SEQUENCES
  wt = get.positions(CDS) %>%
       group_by(orf) %>%
       mutate(bin.pos=dplyr::ntile(wt_pos,100)) %>%
       mutate(aa_pos = ceiling(wt_pos/3), codon_pos = (wt_pos%%3 + (wt_pos%%3==0)*3) ) %>%
       mutate(codon = ifelse(codon_pos==2, paste0(lag(wt_aa),wt_aa,lead(wt_aa)), NA)) %>%
       group_by(orf,aa_pos) %>% fill(codon,.direction = 'updown') %>%
       mutate(codon_aa = Biostrings::GENETIC_CODE[codon])

  # CDS GENOME 1011 STRAINS
  yk11_fastadir = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/transfer_1638744_files_c25fb55c/CDS_withAmbiguityRes/"
  yk11_cds = here("data","YK11-CDS.rds")
  yk11 = preload(yk11_cds,{load.1011.strains(seqdir = yk11_fastadir,seqtype='DNA')},"loading genome cds from 1011 strains...")
  aligned_yk11 = load_aligned_yk11()
  #zero_len = sapply(yk11,length) == 0
  #yk11= yk11[!zero_len]

  # GENOME DATAFRAME
  # G = tibble( orf=names(yk11),
  #             n_strains=lengths(yk11),
  #             len = widths(yk11)) %>%
    G = tibble( orf=names(yk11),
              n_strains=lengths(yk11),
              len = widths(yk11)) %>%
    left_join(get.width(CDS), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
    mutate( match_wt = len == len.s288c )

  # FIND NUCLEOTIDE VARIANTS
  yk11_nt_var = here("data","YK11-VAR_NT.rds")
  nt_var = preload(yk11_nt_var, {lapply(yk11, get.variants,verbose=F) %>% bind_rows()},"find nucleotide variants...")

  snp_nt = get.variants.to.SNP(nt_var) %>%
    #filter(alt_aa %in% c('A','T','C','G') & ref_aa %in% c('A','T','C','G'))  %>%
    add_count(id,ref_pos,name='nvar') %>%
    group_by(id,ref_pos) %>%
    mutate( alt_cumfr=sum(alt_fr)) %>%
    left_join(G,by=c('id'='orf'))

  all_codons = get_all_codons(shorten=T)

  cds_snp_nt = left_join(snp_nt,wt, by=c('id'='orf','ref_pos'='wt_pos','len.s288c'='len')) %>%
    mutate(wt_low = alt_aa == wt_aa, wt_missing=is.na(wt_aa)) %>%
    mutate(alt_codon =stringi::stri_sub_replace(codon,from=codon_pos,to=codon_pos,replacement=alt_aa)) %>%
    mutate(alt_codon_aa = Biostrings::GENETIC_CODE[alt_codon], synonymous = (codon_aa == alt_codon_aa) ) %>%
    # replace aa by nt in column names
    dplyr::rename(ref_nt=ref_aa,alt_nt=alt_aa,wt_nt=wt_aa,cds_pos=ref_pos, cds_fr=ref_fr) %>%
    left_join(all_codons,by=c('alt_codon'='fuzzy_codon')) %>%
    mutate(rate_synonymous = str_count(pattern = fixed(codon_aa), string=aa_ambiguous) / n_codons ) %>%
    group_by(id,aa_pos) %>%
      mutate(dn     = (rate_synonymous==0) * nvar,
             dn_amb = (rate_synonymous<1) * nvar,
             ds_amb = (rate_synonymous>0) * nvar,
             ds     = (rate_synonymous>0.5) * nvar) %>%
    group_by(id) %>%
      mutate(N_POS = n(), N_VAR = sum(nvar),
             N_VAR_STANDARD = sum(!is_ambiguous),
             N_VAR_AMBIGUOUS = sum(is_ambiguous),
             dN = sum(dn), dS=sum(ds),
             dN_amb = sum(dn_amb), dS_amb = sum(ds_amb)) %>%
      relocate(id,n_strains,len,len.s288c,match_wt,
             cds_pos, bin.pos, wt_nt, wt_low, wt_missing,
             ref_nt,cds_fr,alt_nt,alt_fr,alt_cumfr,
             codon_pos, n_codons, codon, alt_codon, is_ambiguous, codon_ambiguous,
             aa_pos, n_aa, codon_aa,alt_codon_aa, aa_ambiguous,
             synonymous, rate_synonymous,
             nvar, dn,dn_amb, ds,ds_amb,
             N_POS, N_VAR, N_VAR_STANDARD, N_VAR_AMBIGUOUS,
             dN,dN_amb, dS,dS_amb)

  # Save nucleotide polymorphism
  saveRDS(cds_snp_nt,cds_snp_nt.rds)
}

# Counting amino acid polymorphisms per protein
snp_nt_per_orf = cds_snp_nt %>%
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(cds_pos))
snp_nt_per_orf %>% arrange(n_snp)
saveRDS(snp_nt_per_orf,here("data",'YK11-ORF-VAR_NT.rds'))

head(cds_snp_nt)
head(PROT_SNP_AA)

all_snp = left_join(cds_snp_nt,PROT_SNP_AA, by=c('id', "aa_pos"='ref_pos', 'alt_codon_aa' = 'alt_aa'), suffix = c("_nt", "_aa"), )

### Must have run align_s288c_to_1011strains.r !!!
# Filter orf with best match to s288c from all strains
aligned_1011_prot=here("prepare","proteome_aligned_s288c_vs_1011strains.rds")
if( file.exists(aligned_1011_prot) ){
  df_strains_s288c = read_rds(aligned_1011_prot)
  #df_best = df_strains_s288c %>% group_by(ID1) %>% filter(SCORE.BLOSUM100 == max(SCORE.BLOSUM100))

  proteomics_strains = c('AMH','CQC','BTT','CPI','BED','BAN','BPL','CMP')
  best_strain = df_strains_s288c %>%
    group_by(strain) %>%
    summarize(aligned = sum(S),
              norf = n_distinct(ID1),
              proteome_len = sum(L),
              pid_proteome = 100*(aligned/proteome_len),
              nsnp = sum(NS),
              orf_with_snv = sum(NS != 0 & (I+D) == 0),
              orf_with_indel = sum((I+D) > 0)
              ) %>%
    mutate(is_proteomics = strain %in% proteomics_strains) %>%
    arrange(nsnp,desc(pid_proteome),desc(norf))%>%
    left_join(load.peter2018.data(1), by=c('strain'='standardized_name'))

  saveRDS(best_strain, here("prepare","best_strain_aligned_to_s288c_proteome.rds"))

  # Show
  median_snp = median(best_strain$nsnp)
  ggplot(best_strain,aes(fill=is_proteomics)) +
    geom_col(aes(y=reorder(strain,pid_proteome),x=nsnp),orientation='y') +
    geom_vline(xintercept=median_snp,linetype='dashed') +
    facet_wrap(~ecological_origins, scales = 'free_y') + xlab('# of SNP') +
    scale_y_discrete(guide=guide_axis(n.dodge=2))+
    scale_fill_manual(values=c("TRUE"='red',"FALSE"='gray50')) +
    theme(legend.position = 'none', axis.text =  element_text(size = 5))
}else{
  message("First, run align_s288c_to_1011strains.r to align all strains to S288C reference strain!")
}


#### 8 strains for proteomics/RNASeq/RiboSeq
S288C = load.sgd.proteome(withORF=T,rm.stop=F) # Reference SGD protein sequences
CDS = load.sgd.CDS(withORF=T) # Reference SGD DNA coding sequences
WT = get.positions(S288C) %>% group_by(orf) %>% mutate(bin.pos=dplyr::ntile(wt_pos,100))
wt = get.positions(CDS) %>% group_by(orf) %>% mutate(bin.pos=dplyr::ntile(wt_pos,100)) %>%
  mutate(aa_pos = ceiling(wt_pos/3), codon_pos = (wt_pos%%3 + (wt_pos%%3==0)*3) ) %>%
  mutate(codon = ifelse(codon_pos==2, paste0(lag(wt_aa),wt_aa,lead(wt_aa)), NA)) %>%
  group_by(orf,aa_pos) %>% fill(codon,.direction = 'updown') %>%
  mutate(codon_aa = Biostrings::GENETIC_CODE[codon])


YK11_prot = read_rds(here("data","YK11-PROT.rds"))
YK11_cds = read_rds(here("data","YK11-CDS.rds"))
riboseq_strains = c('AMH','BAN','BED','BPL','BTT','CMP','CPI','CQC') # Strains with riboseq data (on 14/01/21)
strains.info = load.peter2018.data(1) %>%  # strains info from supp mat of Science paper
  mutate( has_riboseq = standardized_name %in% riboseq_strains)

Y8 = lapply(YK11_prot, function(E){ E[get.strain_orf(E,"strains") %in% riboseq_strains] }) %>% purrr::compact()
saveRDS(Y8, here("data","Y8-PROT.rds"))
y8 = lapply(YK11_cds, function(E){ E[get.strain_orf(E,"strains") %in% riboseq_strains] }) %>% purrr::compact()
saveRDS(y8,here("data","Y8-CDS.rds"))


# snp nucleotides
y8_nt_var = preload(here("data","Y8-VAR_NT.rds"),
                    {lapply(y8, get.variants,verbose=F) %>% bind_rows()},"find nucleotide variants...")

# GENOME DATAFRAME
G = tibble( orf=names(y8),n_strains=lengths(y8),len = widths(y8)) %>%
  left_join(get.width(CDS), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
  mutate( match_wt = len == len.s288c )

y8_snp_nt = get.variants.to.SNP(y8_nt_var) %>% add_count(id,ref_pos,name='nvar') %>%
  group_by(id,ref_pos) %>% mutate( alt_cumfr=sum(alt_fr)) %>%
  left_join(G,by=c('id'='orf'))
y8_cds_snp_nt = left_join(y8_snp_nt,wt, by=c('id'='orf','ref_pos'='wt_pos','len.s288c'='len')) %>%
  mutate(wt_low = alt_aa == wt_aa, wt_missing=is.na(wt_aa)) %>%
  mutate(alt_codon =stringi::stri_sub_replace(codon,from=codon_pos,to=codon_pos,replacement=alt_aa)) %>%
  mutate(alt_codon_aa = Biostrings::GENETIC_CODE[alt_codon], synonymous = (codon_aa == alt_codon_aa) ) %>%
  mutate(ambiguous= !(alt_aa %in% Biostrings::DNA_BASES)) %>%
  # replace aa by nt in column names
  dplyr::rename(ref_nt=ref_aa,alt_nt=alt_aa,wt_nt=wt_aa,cds_pos=ref_pos, cds_fr=ref_fr)

saveRDS(y8_cds_snp_nt,here("data",'Y8-SNP_NT.rds'))

# Counting amino acid polymorphisms per protein
y8_snp_nt_per_orf = y8_cds_snp_nt %>%
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(cds_pos))
y8_snp_nt_per_orf %>% arrange(n_snp)
saveRDS(y8_snp_nt_per_orf,here("data",'Y8-ORF-VAR_NT.rds'))

# snp amino-acid
# FIND AMINO ACID VARIANTS
Y8_aa_var = here("data","Y8-VAR_AA.rds")
Y8_AA_VAR = preload(Y8_aa_var,{lapply(Y8, get.variants,verbose=F) %>% bind_rows()},"find amino acid variants...")

# PROTEOME DATAFRAME
P = tibble( orf=names(Y8),n_strains=lengths(Y8),len = widths(Y8)) %>%
  left_join(get.width(S288C), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
  mutate( match_wt = len == len.s288c )

# GET AMINO ACID POLYMORPHISMS
Y8_SNP_AA = get.variants.to.SNP(Y8_AA_VAR) %>%
  add_count(id,ref_pos,name='nvar') %>%
  group_by(id,ref_pos) %>%
  mutate( alt_cumfr=sum(alt_fr)) %>%
  left_join(P,by=c('id'='orf'))

Y8_PROT_SNP_AA = left_join(Y8_SNP_AA,WT, by=c('id'='orf','ref_pos'='wt_pos','len.s288c'='len')) %>%
  mutate(wt_low = alt_aa == wt_aa,
         wt_missing=is.na(wt_aa),
         early_stop = (alt_aa == "*" & ref_pos != len),
         dSTI.ref=get.score.mutation(ref_aa,alt_aa),
         dSTI.wt=get.score.mutation(wt_aa,alt_aa))
saveRDS(Y8_PROT_SNP_AA,here("data",'Y8-SNP_AA.rds'))

# Counting amino acid polymorphisms per protein
y8_snp_aa_per_orf = Y8_PROT_SNP_AA %>%
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(ref_pos))
y8_snp_aa_per_orf %>% arrange(n_snp)
saveRDS(y8_snp_aa_per_orf,here("data",'Y8-ORF-VAR_AA.rds'))


###
yk11_snp_nt_per_orf = readRDS(here("data",'YK11-ORF-VAR_NT.rds'))
yk11_snp_nt = readRDS(here('data','YK11-SNP_NT.rds')) %>% left_join(yk11_snp_nt_per_orf)
y8_snp_nt_per_orf = readRDS(here("data",'Y8-ORF-VAR_NT.rds'))
y8_snp_nt = readRDS(here('data','Y8-SNP_NT.rds')) %>% left_join(y8_snp_nt_per_orf)

# Yeast SNP (across 1011 or 8 species) ------------------------------------
sc_snp_nt = left_join(yk11_snp_nt,y8_snp_nt,
                      by=c('id', "len", "len.s288c",  "wt_nt", "wt_low", "wt_missing",
                           "cds_pos", "bin.pos", "ref_nt", "match_wt",
                           "codon_pos", "aa_pos", "codon", "codon_aa",
                           "alt_nt", "alt_codon", "alt_codon_aa"),
                      suffix = c(".yk11", ".y8") ) %>%
            mutate(is_y8 = !is.na(n_snp.y8))
saveRDS(sc_snp_nt, here('data','YEAST_VAR_NT.rds'))

yk11_snp_aa_per_orf = readRDS(here("data",'YK11-ORF-VAR_AA.rds'))
yk11_snp_aa = readRDS(here('data','YK11-SNP_AA.rds')) %>% left_join(yk11_snp_aa_per_orf)

y8_snp_aa_per_orf = readRDS(here("data",'Y8-ORF-VAR_AA.rds'))
y8_snp_aa = readRDS(here('data','Y8-SNP_AA.rds'))  %>% left_join(y8_snp_aa_per_orf)

sc_snp_aa = left_join(yk11_snp_aa,y8_snp_aa,
                       by=c('id', "len", "len.s288c", "wt_low", "wt_missing",
                            'ref_pos', "bin.pos", 'ref_aa', "match_wt", 'alt_aa'),
                       suffix =  c(".yk11", ".y8") ) %>%
            mutate(is_y8 = !is.na(n_snp.y8))
saveRDS(sc_snp_aa, here('data','YEAST_VAR_AA.rds'))

# Yeast AA and NT variants per ORF (1011 or 8 strains) --------------------

y8_orf_aavar = read_rds( here('data','Y8-ORF-VAR_AA.rds') )
y8_orf_ntvar = read_rds( here('data','Y8-ORF-VAR_NT.rds') )
yk11_orf_aavar = read_rds( here('data','YK11-ORF-VAR_AA.rds') )
yk11_orf_ntvar = read_rds( here('data','YK11-ORF-VAR_NT.rds') )

yk11_orf_var = left_join(yk11_orf_aavar,yk11_orf_ntvar, by=c('id','n_strains'),suffix=c('.aa','.nt'))
y8_orf_var = left_join(y8_orf_aavar,y8_orf_ntvar, by=c('id','n_strains'),suffix=c('.aa','.nt'))

yeast_orf_var = left_join(yk11_orf_var,y8_orf_var,
                          by=c('id','len.aa','len.s288c.aa','len.nt','len.s288c.nt'),
                          suffix=c('.yk11','.y8')) %>%
  relocate(id,n_strains.y8,n_strains.yk11,len.s288c.nt,len.s288c.aa, len.aa,len.nt )
saveRDS(yeast_orf_var, here('data','YEAST_ORF_VAR.rds'))


align_cds = function(ID,REF,STRAINS,OUTDIR="/data/benjamin/Evolution/YK11/aligned_cds_s288c/"){
  outfile=sprintf("%s/s288c_yk11_%s_cds.mu",OUTDIR,ID)
  if(!file.exists(outfile)){
    refseq = DNAStringSet(REF[[ID]])
    names(refseq) = paste0("S288C_",ID)
    newseq = c(refseq,STRAINS[[ID]])
    ali = muscle::muscle(newseq, quiet = T, diags=T, maxiters=3)
    Biostrings::writeXStringSet(as(ali,"DNAStringSet"),filepath = outfile)
  }else{
    ali = Biostrings::readDNAMultipleAlignment(outfile)
  }
  return(ali)
}

add_s288c_cds = function(ID,REF,STRAINS,OUTDIR="/data/benjamin/Evolution/YK11/cds_s288c"){
  fastafile=sprintf("%s/s288c_yk11_%s_cds.fa",OUTDIR,ID)
  refseq = DNAStringSet(REF[[ID]])
  names(refseq) = paste0("S288C_",ID)
  newseq = c(refseq,STRAINS[[ID]])


  if(length(refseq) == 1 & length(STRAINS[[ID]])>1){
    Biostrings::writeXStringSet(as(newseq,"DNAStringSet"),filepath = fastafile)
  }else{
    return(refseq)
  }
  return(newseq)
}

add_s288c_prot = function(ID,REF,STRAINS,OUTDIR="/data/benjamin/Evolution/YK11/prot_s288c"){
  fastafile=sprintf("%s/s288c_yk11_%s_prot.fa",OUTDIR,ID)
  refseq = AAStringSet(REF[[ID]])
  names(refseq) = paste0("S288C_",ID)
  newseq = c(refseq,STRAINS[[ID]])

  if(length(refseq) == 1 & length(STRAINS[[ID]])>1){
    Biostrings::writeXStringSet(as(newseq,"AAStringSet"),filepath = fastafile)
  }else{
    return(refseq)
  }
  return(newseq)
}

cds_s288c_yk11 = map(names(CDS), ~add_s288c_cds(.x,CDS,YK11_cds))
prot_s288c_yk11 = map(names(S288C), ~add_s288c_prot(.x,S288C,YK11))

sum( lengths(prot_s288c_yk11) == 1 )

align_prot = function(ID,REF,STRAINS,OUTDIR="/data/benjamin/Evolution/YK11/aligned_prot_s288c/"){
  outfile=sprintf("%s/s288c_yk11_%s_prot.mu",OUTDIR,ID)
  if(!file.exists(outfile)){
    ali = muscle::muscle(newseq, quiet = T, diags=T, maxiters=3)
    Biostrings::writeXStringSet(as(ali,"AAStringSet"),filepath = outfile)
  }else{
    ali = Biostrings::readAAMultipleAlignment(outfile)
  }
  return(ali)
}
msa_nt_yk11 = map(names(CDS), ~align_cds(.x,CDS,YK11_cds))
msa_aa_yk11 = map(names(PROT), ~align_prot(.x,PROT,YK11))

