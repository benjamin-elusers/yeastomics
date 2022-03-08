#install.packages('rjson')
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("src","function_YK11.r"))

S288C = load.sgd.proteome(withORF=T,rm.stop=T) # Reference SGD protein sequences
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

  yk11_seq_dir = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"
  yk11_orf = file.path(yk11_seq_dir,paste0(orf,'.fasta'))
  in_s288c = orf %in% names(S288C)
  in_yk11 = file.exists(yk11_orf)
  has_orf  = in_s288c && in_yk11
  if(!has_orf){
    message(sprintf("missing orf %s (S288C %s YK11 %s)",orf,in_s288c,in_yk11))
    return(NA)
  }
  cat(orf,"\n")
  yk11  = load.proteome(yk11_orf,nostop = F)

  return( c(S288C[orf],yk11) )
}


align_s288c_strains = function(orf){
  #tic('future_map_dfr')
  #future::plan(future::multisession(workers = 14))
  #df_strain = furrr::future_map_dfr(strains, .f = ~score_ali(p1=refseq, p2=strain_seq[.x], mat='BLOSUM100'), .progress = T )
  #cat("\n")
  #toc()
  #tic('map_dfr')
  SEQ = get_s288c_strain(orf)
  tic(sprintf('align %s from s288c against 1011 strains',orf))
  df_strain = map_dfr(2:length(SEQ), function(x){ score_ali(p1=SEQ[1], p2=SEQ[x], mat='BLOSUM100') })
  toc()
  #identical(df_strain,df_strain2)
  # df_closest = df_strain %>%
  #              group_by(ID1) %>%
  #              dplyr::filter( SCORE.BLOSUM100== max(SCORE.BLOSUM100))
  return(df_strain)
}


yk11_seq_dir = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"
yk11_strains = load.peter2018.data(1) %>% pull(standardized_name) %>% sort
yk11_orfs = get.orf.filename(list.files(path=yk11_seq_dir,pattern='.fasta'))

orfs = intersect(names(S288C),yk11_orfs)
#orfs_same_length = orfs[ widths(S288C[orfs]) == widths(YK11[orfs]) ]

res1file=here::here("prepare","s288c_vs_strains_0001-2000.rds")
res2file=here::here("prepare","s288c_vs_strains_2001-4000.rds")
res3file=here::here("prepare","s288c_vs_strains_4001-6000.rds")
res4file=here::here("prepare","s288c_vs_strains_6001-6574.rds")

if(!file.exists(res1file)){
  res1 = pbmcapply::pbmcmapply(align_s288c_strains, orfs[0001:2000], mc.cores = 14)
  write_rds(res1,res1file)
}
if(!file.exists(res2file)){
  res2 = pbmcapply::pbmcmapply(align_s288c_strains, orfs[2001:4000], mc.cores = 14)
  write_rds(res2,res2file)
}
if(!file.exists(res3file)){
  res3 = pbmcapply::pbmcmapply(align_s288c_strains, orfs[4001:6000], mc.cores = 14)
  write_rds(res3,res3file)
}
if(!file.exists(res4file)){
  res4 = pbmcapply::pbmcmapply(align_s288c_strains, orfs[6001:6574], mc.cores = 14)
  write_rds(res4,res4file)
}

# Each column is stored as a list

tmp1 = t(read_rds(res1file)) %>% as_tibble() %>% unnest(cols=everything())
tmp2 = t(read_rds(res2file)) %>% as_tibble() %>% unnest(cols=everything())
tmp3 = t(read_rds(res3file)) %>% as_tibble() %>% unnest(cols=everything())
tmp4 = t(read_rds(res4file)) %>% as_tibble() %>% unnest(cols=everything())

df_strains_s288c = bind_rows(tmp1,tmp2,tmp3,tmp4) %>%
                  mutate( strain = get.strain_orf(ID2,'strains') )

df_best = df_strains_s288c %>% group_by(ID1) %>% filter(SCORE.BLOSUM100 == max(SCORE.BLOSUM100))

best_strain = df_strains_s288c %>%
  group_by(strain) %>%
  summarize(aligned = sum(S),pid.total = mean_(PID1), nsnp = sum(NS), norf = n_distinct(ID1)) %>%
  arrange(desc(pid.total),desc(norf))

peter2018 = load.peter2018.data(1) %>% arrange(total_number_of_sn_ps)
top10 = peter2018 %>% slice_max(order_by=-total_number_of_sn_ps,n=10) %>%
        select(isolate_name,standardized_name,isolation,ecological_origins,geographical_origins,total_number_of_sn_ps,clades)

best_strain %>% slice_max(order_by=pid.total,n=10)
top10

#
# for( O in orfs){
#   SEQ = get_s288c_strain(O)
#   align_s288c_strains(SEQ[[1]],SEQ[[2]])
# }
#
# df =  furrr::future_map2( S288C, function(x){ orf=names(x)})


add_missing_strains = function(BS,all_strains){
  this_orf = get.strain_orf(what='orf',x=names(BS))
  these_strains = get.strain_orf(what='strains',x=names(BS))
  missing_strains = setdiff(all_strains,these_strains)
  n = length(missing_strains)
  if(n==0){
    #message('no missing strains')
    return(BS)
  }else{
    cat(sprintf('%10s: %3s missing strains\r',this_orf,n))
    L = max(unique(width(BS)))
    unknown_seq = paste0(rep("-",L),collapse="")
    missing_seq = setNames(rep(AAStringSet(x=unknown_seq),n),missing_strains)
    return( c(BS,missing_seq) )
  }
}

library(bio3d)
yk11 = sapply(orfs, function(x){ add_missing_strains( BS=get_s288c_strain(x),all_strains=yk11_strains) })

LL  = sapply(yk11,function(x){ unique(nchar(x))})
weird = names( LL[ lengths(LL) > 1] )
# Aligned those where the reference sequence has a different size from the strain sequences
for( w in weird){
  tofasta = sprintf("~/Desktop/%s.fasta",w)
  toaln = sprintf("~/Desktop/%s_aligned.fa",w)
  Biostrings::writeXStringSet(yk11[[w]],filepath=tofasta)
  tmp  = bio3d::read.fasta(tofasta)
  aln  = bio3d::seqaln(tmp,outfile = toaln)
  unlink(tofasta)
  tmp = Biostrings::readAAMultipleAlignment(toaln,format='fasta')
  yk11[[w]] = tmp
}

# Following should give empty results

LL  = sapply(yk11,function(x){ unique(nchar(x))})
weird = names( LL[ lengths(LL) > 1] )

final_dir = "/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/1011G/fasta_strain_with_s288c/"
for( o in names(yk11) ){
  tofasta = sprintf("%s/%s.fasta",final_dir,o)
  print(tofasta)
  as(yk11[[o]], "AAStringSet")
  Biostrings::writeXStringSet(,filepath=tofasta,format = 'fasta')
}

