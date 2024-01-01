#source("src/utils.r",local = T)
library(tidyverse)
library(Biostrings)
# Pairwise alignment -----------------------------------------------------------
get.pair.prot = function(prot, pair){
  seq.1 = prot[pair[[1]]]
  seq.2 = prot[pair[[2]]]
  return(list(s1=seq.1,s2=seq.2))
}

pairwise.alignment.to.df  = function(ali){
  library(dplyr)
  N=length(ali)
  ali.list=list()
  doing="Global alignnment of pairwise sequences..."
  message(doing)
  tic(doing)
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
  toc()
  return( bind_rows(ali.list) )
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

get_identity = function(ali){
  sapply(paste0("PID",1:4),function(x){ pid(ali,x) })
}

count_gaps = function(ali,bycol=T){

  if( !bycol ){
    pgaps = rep(0,nrow(ali))
    freqmat=alphabetFrequency(ali,as.prob = T)
    return(freqmat[, "-"])
  }else{
    pgaps = rep(0,ncol(ali))
    consmat = consensusMatrix(ali,as.prob = T)
    if( "-" %in% rownames(consmat) ){
      pgaps = consmat["-",]
    }
    return(pgaps)
  }
}

count_gaps_fromfile = function(ali,bycol=T){
  alignment = Biostrings::readAAMultipleAlignment(ali,'fasta')
  return( count_gaps(alignment,bycol) )
}


get_overlap = function(ali){
  p=unaligned(pattern(ali))
  s=unaligned(subject(ali))
  L1 = width(p)
  L2 = width(s)
  OL1 = L1/L2
  OL2 = L2/L1
  OL  = min(L1,L2) / max(L2,L1)

  overlap = as.double( c(OL,OL1,OL2)*100.0 )

  return(setNames(overlap,c("OL","OL1","OL2")))
}


seq2df = function(BS){
  is_string_set = class(BS) == 'AAStringSet'
  has_one_seq = length(BS)
  if( !is_string_set ) { stop("requires a 'AAStringSet' object!") }
  if( !has_one_seq   ) { stop("must contain 1 sequence at most!") }
  return(
    tibble(
      id = names(BS),
      resi = 1:width(BS),
      resn = unlist(strsplit( toString(BS), "" ))
    )
  )
}

count_pair_ali = function(ali){
  S1 = seq2char(alignedPattern(ali))
  S2 = seq2char(alignedSubject(ali))
  L=nchar(ali)
  L1= pattern(ali) |> unaligned() |> width()
  L2= subject(ali) |> unaligned() |> width()
  N=nmatch(ali)
  NS=nmismatch(ali)
  GAP = sum((S1=="-" & S2=="-"))
  INDEL = sum(xor(S1!='-',S2!='-'))
  IN = sum(S1!="-" & S2=="-")
  DEL = sum(S1=="-" & S2!="-")
  return(c('S'=N,'NS'=NS,'G'=INDEL,"I"=IN,"D"=DEL,'L'=L,'L1'=L1,'L2'=L2,"gapped"=GAP))
}

msa2df = function(MSA_SEQ,REF_NAME,ID=NULL,verbose=F){

  if(is(MSA_SEQ,"AAStringSet")){
    if( n_distinct(widths(MSA_SEQ)) != 1){
      warning('Sequences not aligned! Aligning them with muscle...')
      MSA_SEQ = muscle::muscle(MSA_SEQ,quiet = T,diags=T)
    }else{
      MSA_SEQ = AAMultipleAlignment(MSA_SEQ)
    }
  }

  NC = nchar(MSA_SEQ)
  NS = nrow(MSA_SEQ)
  NAMES = rownames(MSA_SEQ)
  IREF = 1 ## default to use the first sequence in the alignment

  if(is.null(ID)){ stop("The msa should have an identifier (e.g. filename)") }

  if( is.character(REF_NAME) ){
    IREF = grep(REF_NAME,NAMES)
    if(length(IREF) == 0){
      IREF = 1
      warning(sprintf('Reference name %s not found.\n',REF_NAME))
      warning(sprintf('First sequence used as reference %s...',NAMES[IREF]))
    }else if(length(IREF)>1){
      warning(sprintf('Refrence name matched %s sequence identifiers.\n',length(IREF)))
      warning(sprintf("%s\n"), toString(NAMES[IREF]))
      IREF = select.list(title = 'pick the reference sequence to use:',
                         multiple=F, preselect=NAMES[IREF][1], choices=NAMES[IREF])
    }
  }else if( is_integer(REF_NAME) && between(REF_NAME,1,NS) ){
    IREF = REF_NAME
  }else{
    IREF=1
  }
  REF = NAMES[IREF]

  if(verbose)
    message(paste("msa=",ID," ref ->", REF))

  MSAFULL = as.matrix(MSA_SEQ)
  REFSEQ = MSAFULL[IREF,]
  MSA_NOREF = MSAFULL[-IREF,]
  REFMAT = matrix(REFSEQ,nrow = NS-1, ncol=NC,byrow=T)

  CONSENSUS =  consensusMatrix(MSA_SEQ)
  iref_aa=match(REF,rownames(MSA_SEQ))

  df_msa = tibble( idfile = ID,
                   id     = REF,
                   ref_aa  = REFSEQ,
                   ref_gap = REFSEQ == "-",
                   msa_pos = 1:NC,
                   matched = CONSENSUS[iref_aa + (seq_along(REFSEQ)-1) * nrow(CONSENSUS)],
                   mismatched = colSums(REFMAT != MSA_NOREF & MSA_NOREF!="-" & REFMAT!="-"),
                   indel = colSums( xor(MSA_NOREF=='-',REFMAT=="-") ),
                   ins = colSums(REFMAT=="-" & MSA_NOREF!="-"),
                   del  = colSums(REFMAT!="-" & MSA_NOREF=="-")
    ) %>%
    group_by(ref_gap) %>% mutate(ref_pos = ifelse(ref_gap, NA, row_number() ) ) %>%
    rowwise() %>% mutate( total = sum(c_across(cols = c(indel,matched,mismatched))) )

  return(df_msa)
}


score_ali = function(p1,p2,s1,s2, mat="BLOSUM62", opening=10, extend=4){

  ali = pairwiseAlignment(p1, p2, substitutionMatrix = mat, gapOpening=opening, gapExtension=extend)

  if(missing(s1)){ s1=get.strain_orf(names(p1),'both') }
  if(missing(s2)){ s2=get.strain_orf(names(p2),'strain') }

  IDS=c(ID1=names(p1),ID2=names(p2))
  PID = get_identity(ali)
  OL = get_overlap(ali)
  SCORE = setNames(nm=paste0("SCORE.",mat),score(ali))
  MATCH = count_pair_ali(ali)

  # Best way to create a one-row dataframe from a list of named vectors
  row_data = tibble::lst(IDS,PID,OL,SCORE,MATCH) %>%
    purrr::flatten() %>%
    tibble::as_tibble()
  return(row_data)
}

