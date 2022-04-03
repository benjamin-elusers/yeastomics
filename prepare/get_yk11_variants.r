#install.packages('rjson')
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("src","function_YK11.r"))


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
  write_rds(PROT_SNP_AA,here("data",'YK11-SNP_AA.rds'))
}

# Counting amino acid polymorphisms per protein
snp_aa_per_orf = PROT_SNP_AA %>%
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(ref_pos))
snp_aa_per_orf %>% arrange(n_snp)
write_rds(snp_aa_per_orf,here("data",'YK11-ORF-VAR_AA.rds'))

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
  #zero_len = sapply(yk11,length) == 0
  #yk11= yk11[!zero_len]

  # GENOME DATAFRAME
  G = tibble( orf=names(yk11),
              n_strains=lengths(yk11),
              len = widths(yk11)) %>%
    left_join(get.width(yk11), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
    mutate( match_wt = len == len.s288c )

  # FIND NUCLEOTIDE VARIANTS
  yk11_nt_var = here("data","YK11-VAR_NT.rds")
  nt_var = preload(yk11_nt_var, {lapply(yk11, get.variants,verbose=F) %>% bind_rows()},"find nucleotide variants...")

  snp_nt = get.variants.to.SNP(nt_var) %>%
    add_count(id,ref_pos,name='nvar') %>%
    group_by(id,ref_pos) %>%
    mutate( alt_cumfr=sum(alt_fr)) %>%
    left_join(G,by=c('id'='orf'))

  cds_snp_nt = left_join(snp_nt,wt, by=c('id'='orf','ref_pos'='wt_pos','len.s288c'='len')) %>%
    mutate(wt_low = alt_aa == wt_aa, wt_missing=is.na(wt_aa)) %>%
    mutate(alt_codon =stringi::stri_sub_replace(codon,from=codon_pos,to=codon_pos,replacement=alt_aa)) %>%
    mutate(alt_codon_aa = Biostrings::GENETIC_CODE[alt_codon], synonymous = (codon_aa == alt_codon_aa) ) %>%
    # replace aa by nt in column names
    dplyr::rename(ref_nt=ref_aa,alt_nt=alt_aa,wt_nt=wt_aa,cds_pos=ref_pos, cds_fr=ref_fr)

  # Save nucleotide polymorphism
  write_rds(cds_snp_nt,cds_snp_nt.rds)
}

# Counting amino acid polymorphisms per protein
snp_nt_per_orf = cds_snp_nt %>%
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(cds_pos))
snp_nt_per_orf %>% arrange(n_snp)
write_rds(snp_nt_per_orf,here("data",'YK11-ORF-VAR_NT.rds'))

head(cds_snp_nt)
head(PROT_SNP_AA)

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

  write_rds(best_strain, here("prepare","best_strain_aligned_to_s288c_proteome.rds"))

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
