library(Biostrings)
library(tidyverse)
library(stringr)
library(hutils)
library(hablar)

### Strains --------------------------------------------------------------------

# Get strains with identical residues at a certain position in sequences
which.strains = function(orf,pos,letter){
  if(length(pos)>1){ warning("More than one position provided, only the first will be used!") }
  if(any(pos[1]>widths(orf))){ stop(sprintf("Position %s is out of bounds!",pos[1])) }
  matched =  names(orf)[hasLetterAt(orf,letter[1],at=pos[1])]
  matched_strains = get.strain_orf(matched,what='strains')
  return(matched_strains)
}

# Extract the isolate id (standardized name) or the orf from string in format "strain_orf_gene"
get.strain_orf= function(x,what=c('orf','strains','both')){
  #library(stringr)
  #library(hutils)
  if (!is.character(x) && !is(x, "XStringSet"))
    stop("'x' must be a character vector, or a XStringSet object")
  nms = x
  if ( is(x, "XStringSet")) { nms <- names(x) }

  genericName = hutils::longest_suffix(nms)
  orf = stringr::str_extract( string=genericName, pattern="([Y][A-P][LR][0-9]{3}[WC](?:-[A-Z])?)|(Q[0-9]{4})|(R[0-9]{4}[WC])")
  strains = stringr::str_replace(string = nms, pattern = fixed(genericName), replacement = "")
  res = data.frame( strain_orf = orf, strain_name = strains, stringsAsFactors = F)

  what = match.arg(arg=tolower(what),choices = c('orf','strains','both'), several.ok = F)

  if(what == 'orf'){ return(unique(res$strain_orf)) }
  else if(what == 'strains'){ return(res$strain_name)}
  return(res)
}

### Variants = alternative residues --------------------------------------------

# Get the number of variants per position which cumulative frequency is over pmin
get.variants.count = function(S){
  FREQ = consensusMatrix(S, as.prob = T)
  max.freq = apply(FREQ,2,max) # Find max AA frequency
  VARCOUNT = colSums(FREQ[,drop=F]>0) - 1*(max.freq==1) # counting variants per position (should be min 2)
  invariant = VARCOUNT == 1
  if( any(invariant) ){ warning(sprintf("The following positions are invariants (single aa): %s",toString(which(invariant)))) }
# q=1-PMIN
# VAR = max.freq <= q  # variants are AA with max frequency below q = (1-p)
# TEST = replace(VAR,VAR==TRUE,colSums(M[,VAR,drop=F]>0,2))
# if( !all(VARCOUNT==TEST) ){ warning("DEBUG: DID NOT PASSED THE CHECK FOR COUNTING THE VARIANTS") }
  return(VARCOUNT)
}
count.variants = function(S){ sum(get.variants.count(S)) }

# Get the positions of variants which cumulative frequency is over pmin
get.variants.positions = function(S,as.ind=T){
  POSVAR = get.variants.count(S) != 0
  if(as.ind) return(which(POSVAR))
  return(POSVAR)
}

filter.variants.freq = function(S,zero2na=T){
  FR = consensusMatrix(S, as.prob = T)
  varpos = get.variants.positions(S,as.ind=T)
  freq.var = FR[,varpos,drop=F] # drop=F to keep data.frame format even for single columns
  colnames(freq.var) = varpos
  if(zero2na) freq.var[ freq.var == 0 ] = NA
  return(freq.var)
}

# Get the frequencies of variants positions
get.variants.freq = function(S,as.vec=T){
  VARFREQ = filter.variants.freq(S,zero2na = T)
  if(as.vec) return(VARFREQ[!is.na(VARFREQ)])
  return(VARFREQ)
}

# Get the variants residues (positions, residue, frequency, counts)
get.variants= function(S,verbose=T){
  #library(tidyverse)
  id = get.strain_orf(S,what='orf')
  if(verbose){ message(id) }
  if( count.variants(S) > 0  ){
    varfreq = stack( get.variants.freq(S,as.vec=F) )
    variants = na.omit(varfreq)
    colnames(variants) = c('var_aa','var_pos','var_fr')
    VAR = as_tibble(variants) %>%
      mutate(var_pos=as.integer(as.character(var_pos)), var_aa=as.character(var_aa)) %>%
      group_by(var_pos) %>%
      mutate( id=id, maxfreq = var_fr == max(var_fr)) %>%
      dplyr::select(id,everything())
    return(VAR)
  }
}

# Get the early STOP with frequency above pmin
get.STOP = function(S){
  M = consensusMatrix(S, as.prob = T)
  STOP = (M["*",])
  if( M['*',ncol(M)] == 1  ){ STOP[ncol(M)] = FALSE }
  return(1*STOP)
}

# Count number of SNPS along the sequence using a fixed-window non-overlapping window
count.SNP.bypos = function(s,nbin=100,.id=NULL){
  snp = get.variants.count(s)
  earlySTOP = get.STOP(s)
  len =unique(width(s))

  snppos = tibble( orf=.id, bin.pos = dplyr::ntile(1:len, nbin), sites=snp!=0, snp=snp, earlystop=earlySTOP) %>%
    group_by(orf,bin.pos) %>%
    summarize( n_sites=sum(sites), n_snp=sum(snp), n_stop = sum(earlySTOP), n_pos = n() )
  return(snppos)
}

### SNP = segregating sites ----------------------------------------------------

get.maxfreq = function(VAR){
  MAXFREQ = VAR %>%
    dplyr::filter(maxfreq) %>%
    rename_with(~sub(x=.x,"var_","ref_")) %>%
    dplyr::select(-maxfreq)
  return(MAXFREQ)
}

get.lowfreq = function(VAR){
  LOWFREQ = VAR %>%
    dplyr::filter(!maxfreq) %>%
    rename_with(~sub(x=.x,"var_","alt_")) %>%
    dplyr::select(-maxfreq)
  return(LOWFREQ)
}

get.variants.to.SNP = function(VAR){
  ALT = get.lowfreq(VAR)
  REF = get.maxfreq(VAR)
  SNP = left_join(REF, ALT, by = c('id'='id','ref_pos'='alt_pos'))
  return(SNP)
}

get.ALT = function(S,VAR,verbose=F){
  if( count.variants(S) > 0  ){
    ALT = get.lowfreq(VAR = get.variants(S,verbose))
    return(ALT)
  }
}

### REF = most frequent variant
get.REF = function(S,verbose=F){
  if( count.variants(S) > 0  ){
    REF = get.maxfreq(VAR = get.variants(S,verbose))
    return(REF)
  }
}


get.SNP = function(S,verbose=F){
  if( count.variants(S) > 0  ){
    V = get.variants(S,verbose)
    SNP = get.variants.to.SNP(V)
    return(SNP)
  }
}

get.aa.pos = function(S,POS){
  OOB = POS>length(S)
  if(OOB){
    warning(sprintf("The following position is out of bounds: %s",toString(POS)))
  }
  validPOS = POS[-OOB]
  if(is.null(validPOS)){ return(NA) }
  else{
    return(Biostrings::letter(S,validPOS))
  }
}

get.SNPscore = function(S, scoremat = "BLOSUM100", masked=c()){
  library(msa)
  mat.sub = grep( data(package="Biostrings")$results[,'Item'], pattern="(BLOSUM)|(PAM)", value=TRUE)
  scoremat = match.arg(scoremat, choices=mat.sub, several.ok = F)
  ali = msa(S)
  msaConservationScore(ali,scoremat)
}

#http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz
#ref.file = 'data/YK11_REF.rds'
#alt.file = 'data/YK11_ALT.rds'
#snp.file = 'data/YK11_SNP.rds'

#ALT = preload(alt.file, loading.call = map_df(YK11, get.ALT,verbose=F),doing='Extract ALT amino-acid...')
#REF = preload(ref.file, loading.call = map_df(YK11, get.REF,verbose=F),doing='Extract REF amino-acid...')
#SNP = preload(snp.file, loading.call = map_df(YK11, get.SNP,verbose=T), doing='Find all SNPs...')
# WT  = preload(wt.file, loading.call = map2_df(unfactor(SNP$ref_pos),SNP$id,get.aa,S=S288C),
#               doing='Extract WT amino-acid...'
#              )

### Prepare data and plot ------------------------------------------------------
get_snp_orf = function(PSNP,pmin){
  SNPORF = PSNP %>%
            dplyr::filter(alt_fr >= pmin) %>%
            group_by(id) %>%
            summarize(sites=n_distinct(ref_pos), var_count = n(), PMIN=pmin) %>%
            ungroup() %>% mutate(md_sites=median(sites))
  return(SNPORF)
}

get_snp_orf_plot = function(psnp=PROT_SNP,BIN=5,FREQMIN=c(0.0001,1e-3,1e-2,5e-2,1e-1,2e-1,3e-1,4e-1)){
  require(ggplot2)
  require(ggthemes)

  tmp = map_df(.x = FREQMIN, .f = get_snp_orf, PSNP=psnp)

  P = ggplot(tmp,aes(label=sites,x=sites)) +
        geom_histogram(binwidth = BIN, fill=MAIN.COLOR,color='black') +
        geom_vline(aes(xintercept=md_sites),color='red') +
        geom_text(aes(x=md_sites,y=0,label=paste("median =",md_sites)),color='red',size=4,hjust=-0.1,vjust='inward',check_overlap = T) +
        xlab("Number of SNP") + ylab('# ORFS')
  if(length(FREQMIN)>1){ P = P + facet_wrap(~PMIN,scales='free') }
  return(P)
}

###
get_snp_protlen = function(PSNP,pmin){
  snp_by_len = PSNP %>%
    dplyr::filter(alt_fr >= pmin) %>%
    group_by(id) %>%
    summarize(var_count=sum(nvar), sites_count=n_distinct(ref_pos),len = mean(len)) %>%
    mutate(PMIN = factor(pmin)) %>%
    left_join(DESC,by=c('id'='ORF'))

  C=cor.sub.by(snp_by_len,XX='len',YY='sites_count', BY = 'PMIN',na.rm=T) %>%
    mutate( toshow=sprintf("r %.3f\np %.1e\nn %s ",r,p,N) )

  SNPLEN = left_join(snp_by_len,C,c('PMIN'='PMIN'))

  return(SNPLEN)
}

get_snp_protlen_plot = function(SNPLEN){
  require(ggplot2)
  require(ggthemes)
  P = ggplot(SNPLEN) +
       geom_point(aes(y=sites_count,x=len,text=id,func=FUNCTION),shape=21,size=0.8) +
       geom_text(aes(y=0,x=5000,label=toshow),hjust='inward',vjust='inward',check_overlap =T,size=4)+
       xlab("Protein Length (AA)") + ylab('SNP count')
  return(P)
}


get_snp_protlen_anim = function(animated=T, psnp = PROT_SNP, FREQMIN=c(0.0001,1e-3,1e-2,5e-2,1e-1,2e-1,3e-1,4e-1), ...){
  library(gganimate)
  library(gifski)
  #require(magick)
  #require(av)
  tmp = map_df(.x = FREQMIN, .f = get_snp_protlen, PSNP=psnp)
  anim =  get_snp_protlen_plot(tmp) +
    theme(legend.position='none') +
    geom_text(aes(y=0,x=5000,label=toshow),hjust='inward',vjust='inward',check_overlap =T,size=4)+
    transition_states(PMIN,wrap=F,transition_length = 2, state_length = 3) +
    ggtitle("Minimum SNP% > {closest_state}")+
    view_follow(fixed_y = FALSE)+
    ease_aes()
  if(!animated){ return(anim) }
  animate(plot = anim, ... )
}

###
get_snp_pos = function(PSNP,pmin){
  SNPPOS=snp_by_pos = PSNP %>% dplyr::filter(alt_fr > pmin)%>%
    group_by(bin.pos) %>%
    summarize( sites=n_distinct(id,ref_pos), var_count=n(), PMIN=pmin) %>%
    group_by(PMIN) %>%
    mutate(YMAX=max(var_count))
  return(SNPPOS)
}

get_snp_pos_plot =function(psnp=PROT_SNP,FREQMIN=c(0.0001,1e-3,1e-2,5e-2,1e-1,2e-1,3e-1,4e-1)){
  require(ggplot2)
  require(ggthemes)
  tmp = map_df(.x = FREQMIN, .f = get_snp_pos, PSNP=psnp)
  P = ggplot( tmp,aes(x=bin.pos)) +
    geom_line(col='red',aes(y=sites)) + #geom_point(aes(y=sites),col=MAIN.COLOR,size=0.3) +
    annotate('text',col='red',label='SNP count',y=0,x=100,vjust='inward',hjust='inward') +
    geom_line(aes(y=var_count),col='dodgerblue') +# geom_point(aes(y=var_count),col=MAIN.COLOR,size=0.3) +
    geom_text(col='dodgerblue',label='Variants count',aes(y=YMAX),x=100,vjust='inward',hjust='inward',check_overlap = T) +
    #geom_vline(xintercept=seq(5,95,by=5),size=0.1,linetype='dotted') +
    ylab('# SNP') + xlab('Sequence position (%)') +
    theme(legend.position='none')
    if( length(FREQMIN) > 1){ P = P + facet_wrap(~PMIN, scales = 'free_y') }
    return(P)
}

