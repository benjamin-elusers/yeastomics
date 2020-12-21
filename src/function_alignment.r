source("src/utils.r")
# Pairwise alignment -----------------------------------------------------------
get.pair.prot = function(prot, pair){
  seq.1 = prot[pair[[1]]]
  seq.2 = prot[pair[[2]]]
  return(list(s1=seq.1,s2=seq.2))
}

align.pair.prot  = function(p1,p2,mat = "BLOSUM62",tomatrix=F, opening=10, extend=4 ){
  if(length(p1) != length(p2))
    stop("p1 and p2 should have the same length!")
  ali = pairwiseAlignment(p1, p2, substitutionMatrix = mat, gapOpening=opening, gapExtension=extend)
  if(tomatrix){
    return( pairwise.alignment.to.df(ali) )
  }else{
    return(ali)
  }
}

pairwise.alignment.to.df  = function(ali){
  library(dplyr)
  N=length(ali)
  ali.list=list()
  for( i in 1:N ){
    prg = sprintf("%.1f %% [%5s/%5s]        \r", 100*(i/N), i,N)
    cat(prg)
    A = ali[i]
    P=alignedPattern(A)
    S=alignedSubject(A)
    S1=unlist(strsplit(x = toString(P),''))
    S2=unlist(strsplit(x = toString(S),''))
    resi = seq_along(S2)
    m.ali = tibble( pos=resi, aligned = (S1==S2), nid = nmatch(A),
                    s1 = names(P), s2 = names(S),
                    aa.s1 = S1, gap.s1 = (aa.s1=='-'),
                    aa.s2= S2, gap.s2 = (aa.s2=='-'),
                    l.s1 = width(P), l.s2 = width(S),
                    ol.12 = l.s1/l.s2, ol.21 = l.s2/l.s1,
                    pid.aligned = pid(A,type = "PID2" ),
                    pid.short = pid(A,type = "PID3" ),
                    pid.long  = 100* (nid / max(l.s1,l.s2)) ) %>%
      group_by(gap.s2) %>% mutate( pos.s2 = ifelse(gap.s2, NA, no=row_number() ) ) %>%
      group_by(gap.s1) %>% mutate( pos.s1 = ifelse(gap.s1, NA, no=row_number() ) ) %>%
      ungroup() %>% distinct()
    ali.list[[i]] = m.ali
  }
  return( bind_rows(ali.list) )
}

aa2str = function(aacol,by){
  aa.list = split(as.character(aacol), f = by)
  prot = Biostrings::AAStringSet( sapply(aa.list,paste0,collapse='') )
  return(prot)
}

list2str = function(L){
  if(!is.list(L))
    stop("Input must be a list of vector!")
  bs  = Biostrings::BString(sapply(L,function(x){ print(x); paste0(x,collapse='') }) )
}

compute.aascore.ali = function(df.ali,aascores=sticky){
  df.score = tibble( aa=names(aascores), sti=aascores)
  df.ali %>%
    left_join(df.score, by=c('s1'='aa')) %>% rename(sti.1 = sti ) %>%
    left_join(df.score, by=c('s2'='aa')) %>% rename(sti.2 = sti )
}


