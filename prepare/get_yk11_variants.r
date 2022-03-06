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
            len = sapply(widths(YK11),unique)) %>%
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

orfs = names(S288C)
orf='YAL001C'
align_s288c_strains = function(orf){
  refseq = S288C[orf]
  strain_seq = YK11[[ orfs[2] ]]
  strain_names = get.strain_orf(names(strainseq),what = 'strains')

  df_ali = lapply(names(strainseq), function(x){ align_pair(refseq,strain_seq[[x]],tomatrix = T) })
}

tmp = align_s288c_strains('YAL001C')

length(tmp)
tmp[[1]]

function(p1,p2, mat="BLOSUM62", opening=10, extend=4){

  ali = pairwiseAlignment(p1, p2, substitutionMatrix = mat, gapOpening=opening, gapExtension=extend)

  get_pid = function(ali){
    sapply(paste0("PID",1:4),function(x){ pid(ali,x) })
  }

  get_ol = function(ali){
    p=unaligned(pattern(ali))
    s=unaligned(subject(ali))
    L1 = width(p)
    L2 = width(s)
    OL1 = L1/L2
    OL2 = L2/L1
    OL  = min(L1,L2) / max(L2,L1)
    return(c("OL1"=OL1,"OL2"=OL2,"OL"=OL))
  }

  get_char = function(seq){
    return( unlist(strsplit(as.character(seq),"")) )
  }

  get_matches = function(ali){
    S1 = get_char(alignedPattern(ali))
    S2 = get_char(alignedSubject(ali))
    L=nchar(ali)
    N=nmatch(ali)
    NS=nmismatch(ali)
    GAP = sum((S1=="-" & S2=="-"))
    INDEL = sum(xor(S1!='-',S2!='-'))
    IN = sum(S1!="-" & S2=="-")
    DEL = sum(S1=="-" & S2!="-")
    return(c('S'=N,'NS'=NS,'G'=INDEL,"I"=IN,"D"=DEL,'L'=L,"gapped"=GAP))
  }

  PID = get_pid(ali)
  OL = get_ol(ali)
  SCORE = score(ali)
  MATCH = get_matches(ali)
  names(p1)
  get.strain_orf(names(p2),"strain")

  resi = seq_along(S2)
  tibble(G=,S=nmatch(A),NS= nid =, pid= pid(), ol =

                s1 = names(P), s2 = names(S),
                aa.s1 = S1, gap.s1 = (aa.s1=='-'),
                aa.s2= S2, gap.s2 = (aa.s2=='-'),
                l.s1 = width(p), l.s2 = width(s), l=width(P),
                ol.12 = l.s1/l.s2, ol.21 = l.s2/l.s1,
                pid.aligned = pid(A,type = "PID2" ),
                pid.short = pid(A,type = "PID3" ),
                pid.long  = 100* (nid / max(l.s1,l.s2)) ) %>%
  group_by(gap.s2) %>% mutate( pos.s2 = ifelse(gap.s2, NA, no=row_number() ) ) %>%
  group_by(gap.s1) %>% mutate( pos.s1 = ifelse(gap.s1, NA, no=row_number() ) ) %>%
  ungroup() %>% distinct()
ali.list[[i]] = m.ali