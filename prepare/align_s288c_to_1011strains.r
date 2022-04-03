#### LOAD PACKAGES ####
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("src","function_YK11.r"))
library(bio3d)
#### FUNCTIONS ####
add_missing_strains = function(BS,all_strains){
  this_orf = get.strain_orf(what='orf',x=names(BS))
  these_strains = get.strain_orf(what='strains',x=names(BS))
  missing_strains = setdiff(all_strains,these_strains)
  n = length(missing_strains)
  if(n==0){
    #message('no missing strains')
    return(BS)
  }else{
    suffix=hutils::longest_suffix(x = names(BS)[-1])
    cat(sprintf('%10s: %3s missing strains (%s)\r',this_orf,n,suffix))
    L = max(width(BS))
    unknown_seq = paste0(rep("-",L),collapse="")
    missing_seq = setNames(rep(AAStringSet(x=unknown_seq),n),paste0(missing_strains,suffix))
    return( c(BS,missing_seq) )
  }
}

align_pair  = function(p1,p2,mat = "BLOSUM62",tomatrix=F, opening=10, extend=4 ){
  ali = pairwiseAlignment(p1, p2, substitutionMatrix = mat, gapOpening=opening, gapExtension=extend)
  if(tomatrix){
    return( pairwise.alignment.to.df(ali) )
  }else{
    return(ali)
  }
}

get_s288c_strain = function(orf,yk11_dir= "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"){
  yk11_orf = file.path(yk11_dir,paste0(orf,'.fasta'))
  in_s288c = orf %in% names(S288C)
  in_yk11 = file.exists(yk11_orf)
  has_orf  = in_s288c && in_yk11
  if(!has_orf){
    message(sprintf("missing orf %s (S288C %s YK11 %s)",orf,in_s288c,in_yk11))
    return(NA)
  }
  #cat(orf,"\n")
  yk11  = load.proteome(yk11_orf,nostop = F)
  wt = S288C[orf]
  names(wt) = paste0('S288C_',names(wt))
  return( c(wt,yk11) )
}

align_s288c_strains = function(orf){
  SEQ = get_s288c_strain(orf)
  tic(sprintf('align %s from s288c against 1011 strains',orf))
  df_strain = map_dfr(2:length(SEQ), function(x){ score_ali(p1=SEQ[1], p2=SEQ[x], mat='BLOSUM100') })
  toc()
  return(df_strain)
}

#### MAIN ####
S288C = load.sgd.proteome(withORF=T,rm.stop=F) # Reference SGD protein sequences
CDS = load.sgd.CDS()
yk11_strains = load.peter2018.data(1) %>% pull(standardized_name) %>% sort
snp_count_per_orf=read_rds(here("data",'YK11-ORF-VAR_AA.rds'))

yk11_cds = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/transfer_1638744_files_c25fb55c/CDS_withAmbiguityRes/"
cds = read.sequences(seqfiles = list.files(yk11_cds,pattern='.fasta',full.names = T),type='DNA',strip.fname = T)

yk11_prot = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"
yk11_orfs = get.orf.filename(list.files(path=yk11_prot,pattern='.fasta'))

orfs = intersect(names(S288C),yk11_orfs)
#orfs_same_length = orfs[ widths(S288C[orfs]) == widths(YK11[orfs]) ]

## Add S288C and gapped sequence for missing strains ##
yk11 = sapply(orfs, function(x){ add_missing_strains( BS=get_s288c_strain(x),all_strains=yk11_strains) })

# Remove the stop codon at the end of the sequences (keep the same length for all sequences)
yk11_nostop = sapply(yk11,stripR)
yk11=yk11_nostop

## Find the 1000 most variable orfs (e.g. containing the most genetic variations) ##
most_var = snp_count_per_orf %>% ungroup() %>%
  dplyr::mutate( f_var = n_var/len,
                 rk_var = dense_rank(-f_var),
                 pc_var = percent_rank(-f_var)) %>%
  arrange(rk_var) %>% dplyr::filter(rk_var<=1000) %>%
  mutate(total_nvar = sum(n_var), total_len = sum(len),total_fvar =100*total_nvar/total_len)
orf_var = sort(unique(most_var$id))
orf_var = names(yk11)

## Find orfs not identical to s288c ##
LL  = sapply(yk11,function(x){ unique(nchar(x))})
weird = names( LL[ lengths(LL) > 1] )
# Aligned those where the reference sequence has a different size from the strain sequences
s288c_alidir = "/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/1011G/strains_s288c_not_aligned"
realigned=list()
for( w in intersect(weird,orf_var) ){
  tofasta = sprintf("%s/%s.fasta",s288c_alidir,w)
  toaln = sprintf("%s/%s_aligned.fa",s288c_alidir,w)
  Biostrings::writeXStringSet(yk11[[w]],filepath=tofasta)
  tmp  = bio3d::read.fasta(tofasta)
  aln  = bio3d::seqaln(tmp,outfile = toaln )
  unlink(tofasta)
  tmp = Biostrings::readAAMultipleAlignment(toaln,format='fasta')
  realigned[[w]] = tmp
}

# Following should give empty results
LL_ali  = sapply(yk11,function(x){ unique(nchar(x))})
weird_ali = names( LL_ali[ lengths(LL_ali) > 1] )

## Realign the orf with not identical length to s288c reference proteome ##
check_realigned = c()
for( orf in names(realigned) ){
  BS = as(realigned[[orf]],'AAStringSet')
  same_length = length(unique(widths(BS))) == 1
  internal_gaps = as.character(BS[[orf]]) %>% str_extract_all(pattern = "[A-Z]-{1,}[A-Z]") %>% unlist
  len_indel = nchar(internal_gaps) - 2
  check_realigned[orf] = same_length & sum(len_indel>1) < 1
}

## Output the fasta files with s288c and all the 1011 strains
# Gapped sequence == missing ORF in strains
#final_dir = "/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/1011G/top1000_strain_with_s288c/"
final_dir = "/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/1011G/fasta_strain_with_s288c/"
dir.create(final_dir)
for( o in names(yk11) ){

  strains = get.strain_orf(names(yk11[[o]]),what = 'strain')
  strains[is.na(strains)] = 'S288C'

  if( o %in% names(realigned) ){
    # WEIRD ORF (NOT UNIQUE LENGTH)
    if( check_realigned[o] ){
      # REALIGNED SEEMS FINE (SAME LENGTH, NOT TOO MANY (LONG) INDELS)
      tofasta = sprintf("%s/%s_realigned.fasta",final_dir,o)
      final_ali = add_missing_strains(as(realigned[[o]],'AAStringSet'),all_strains=yk11_strains)
      #print(unique(widths(final_ali)))
    }else{
      print(sprintf("AVOIDING WEIRD ALIGNED ORF %s",o))
    }
  }else{
    # NORMAL ORFS
    tofasta = sprintf("%s/%s_strains_only.fasta",final_dir,o)
    final_ali = yk11[[o]]
  }

  names(final_ali) = strains

  Biostrings::writeXStringSet(final_ali,filepath=tofasta,format = 'fasta')
}

#### TEST ####
#test = read.sequences(seqfiles = list.files(final_dir,pattern='.fasta',full.names = T))
#tmp = lapply(test, function(x){ unique(widths(x)) })
#test[ lengths(test) <1012 ]

#dir1011="/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/1011G/"
#bigali.idx=read_tsv(file.path(dir1011,'superalignment.idx'))
#bigali.seq=load.proteome(file.path(dir1011,'superalignment.fa'))

# Compare aligned proteome from s288c vs from 1011 strains with key statistics
#res1file=here::here("prepare","s288c_vs_strains_0001-2000.rds")
#if(!file.exists(res1file)){
# res1 = pbmcapply::pbmcmapply(align_s288c_strains, orfs[0001:2000], mc.cores = 14)
# write_rds(res1,res1file)
#}
# Each column is stored as a list
#tmp1 = t(read_rds(res1file)) %>% as_tibble() %>% unnest(cols=everything())
# df_strains_s288c = bind_rows(tmp1,tmp2,tmp3,tmp4) %>%
#                   mutate( strain = get.strain_orf(ID2,'strains') )
# write_rds(df_strains_s288c, here("prepare","proteome_aligned_s288c_vs_1011strains.rds"))
