#install.packages('rjson')
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("src","function_YK11.r"))

S288C = load.sgd.proteome(withORF=T,rm.stop=F) # Reference SGD protein sequences
WT = get.positions(S288C) %>%
  group_by(orf) %>%
  mutate(bin.pos=dplyr::ntile(wt_pos,100))

yk11_seq_file = here("data","proteome-1011-strains.rds")
if( file.exists(yk11_seq_file) ){
  YK11 = readRDS(yk11_seq_file)  # 1011 strains proteomes sequences
}else{
  yk11_seq_dir = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"
  YK11 = load.1011.strains(seqdir = yk11_seq_dir)
}
P = tibble( orf=names(YK11),
            n_strains=lengths(YK11),
            len = widths(YK11)) %>%
    left_join(get.width(S288C), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
    mutate( match_wt = len == len.s288c )

yk11_var_file = here("data","YK11_VAR.rds")
if( file.exists(yk11_var_file) ){
  VAR = readRDS(yk11_var_file)
}else{
  VAR = purrr::map_df(YK11, get.variants,verbose=F)
}

SNP = get.variants.to.SNP(VAR) %>%
  add_count(id,ref_pos,name='nvar') %>%
  group_by(id,ref_pos) %>%
  mutate( alt_cumfr=sum(alt_fr)) %>%
  left_join(P,by=c('id'='orf'))

PROT_SNP = left_join(SNP,WT, by=c('id'='orf','ref_pos'='wt_pos','len.s288c'='len')) %>%
  mutate(wt_low = alt_aa == wt_aa,
         wt_missing=is.na(wt_aa),
         early_stop = (alt_aa == "*" & ref_pos != len),
         dSTI.ref=get.score.mutation(ref_aa,alt_aa),
         dSTI.wt=get.score.mutation(wt_aa,alt_aa))

snp_count_per_orf = PROT_SNP %>%
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(ref_pos))
snp_count_per_orf %>% arrange(n_snp)

write_rds(PROT_SNP,here("data",'YK11-SNP.rds'))
write_rds(snp_count_per_orf,here("data",'YK11-ORF-VAR.rds'))

PROT_SNP=read_rds(here("data",'YK11-SNP.rds'))
snp_count_per_orf=read_rds(here("data",'YK11-ORF-VAR.rds'))


align_pair  = function(p1,p2,mat = "BLOSUM62",tomatrix=F, opening=10, extend=4 ){
  ali = pairwiseAlignment(p1, p2, substitutionMatrix = mat, gapOpening=opening, gapExtension=extend)
  if(tomatrix){
    return( pairwise.alignment.to.df(ali) )
  }else{
    return(ali)
  }
}

get_s288c_strain = function(orf){
  in_s288c = orf %in% names(S288C)
  in_yk11 = orf %in% names(YK11)
  has_orf  = in_s288c && in_yk11
  if(!has_orf){
    message(sprintf("missing orf %s (S288C %s YK11 %s)",orf,in_s288c,in_yk11))
    return(NA)
  }
  cat(orf,"\n")

  return( AAStringSetList(S288C[orf],YK11[[orf]]) )
}


align_s288c_strains = function(orf){
  #tic('future_map_dfr')
  #future::plan(future::multisession(workers = 14))
  #df_strain = furrr::future_map_dfr(strains, .f = ~score_ali(p1=refseq, p2=strain_seq[.x], mat='BLOSUM100'), .progress = T )
  #cat("\n")
  #toc()
  #tic('map_dfr')
  SEQ = get_s288c_strain(orf)
  df_strain = sapply(SEQ, function(x){ print(x); score_ali(p1=SEQ[[1]], p2=x, mat='BLOSUM100') })
  #toc()
  #identical(df_strain,df_strain2)
  df_closest = df_strain %>%
               group_by(ID1) %>%
               dplyr::filter( SCORE.BLOSUM100== max(SCORE.BLOSUM100))
  return(df_closest)
}


future::plan(future::multisession(workers = 14))
orfs = intersect(names(S288C),names(YK11))

pbmcapply::pbmcmapply(align_s288c_strains, orfs)
for( O in orfs){
  SEQ = get_s288c_strain(O)
  align_s288c_strains(SEQ[[1]],SEQ[[2]])
}

df =  furrr::future_map2( S288C, function(x){ orf=names(x)})
