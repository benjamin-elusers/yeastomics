library(Biostrings)

get.SNP.positions = function(S,pmin=0.05,as.ind=T){
  M = consensusMatrix(S, as.prob = T)
  q=1-pmin
  WT.freq = apply(M,2,max) # Find max AA frequency (likely wildtype)
  SNP = WT.freq <= q  # SNP are AA with max frequency below q = (1-p)
  if(as.ind) SNP = which(SNP)
  return(SNP)
}

count.SNP.bypos = function(s,p,nbin=100,.id=NULL){
  # if no id given, extract the yeast ORF (locus name) from the sequnece names
  if(is.null(.id)){ .id = get.strain_orf(s,what='orf') }
  snp = get.SNP.positions(s,p,as.ind=F)
  snppos = tibble( orf=.id, bin.pos = dplyr::ntile(seq_along(snp), nbin), snp=snp*1) %>%
    group_by(orf,bin.pos) %>%
    summarize( n_snp=sum(snp))
  return(snppos)
}

get.noSNP = function(S){ !get.SNP.positions(S,pmin= 1e-7,as.ind=F) } # the lowest aa frequency allowed is 1e-7 (7 digits precision)
count.SNP = function(s,p){ sum(get.SNP.positions(S=s,pmin=p,as.ind=F)) }
count.STOP = function(s,p){ sum(get.SNP.positions(S=s,pmin=p,as.ind=F)) }

freq.SNP = function(s,p,zero2na=T){
  M = consensusMatrix(s, as.prob = T)
  pos = get.SNP.positions(S=s,pmin=p,as.ind=T)
  freq.snp = M[,pos,drop=F] # drop=F to keep data.frame format even for single columns
  colnames(freq.snp) = pos
  if(zero2na) freq.snp = replace(freq.snp,freq.snp==0,NA)
  return(freq.snp)
}


### REF = Wildtype or most frequent amino-acid
get.REF.aa = function(s,p){
  FR=freq.SNP(s,p)
  AA = rownames(FR)
  REF = AA[apply(FR,2,which.max)]
  names(REF)=names(FR)
  return(REF)
}

get.REF.freq = function(s,p){
  FR=freq.SNP(s,p,zero2na = T)
  REF = apply(FR,2,max,na.rm=T)
  names(REF)=colnames(FR)
  return(REF)
}

is.ref = function(s,p){
  FR=freq.SNP(s,p,zero2na = F)
  is_ref = apply(FR,2,function(x){ x==max(x,na.rm=T) })
  return(is_ref)
}

get.strain_orf= function(x,what=c('orf','strains','both')){
  library(stringr)
  library(hutils)
  if (!is.character(x) && !is(x, "XStringSet"))
    stop("'x' must be a character vector, or a XStringSet object")
  nms = x
  if ( is(x, "XStringSet")) { nms <- names(x) }

  genericName = hutils::longest_suffix(nms)
  orf = stringr::str_extract( string=genericName, pattern="([Y][A-P][LR][0-9]{3}[WC](?:-[A-Z])?)|(Q[0-9]{4})|(R[0-9]{4}[WC])")
  strains = stringr::str_replace(string = nms, pattern = genericName, replacement = "")
  res = data.frame( strain_orf = orf, strain_name = strains, stringsAsFactors = F)

  what = match.arg(arg=tolower(what),choices = c('orf','strains','both'), several.ok = F)

  if(what == 'orf'){ return(unique(res$strain_orf)) }
  else if(what == 'strains'){ return(res$strain_name)}
  return(res)
}

get.REF = function(s,p,.id=NULL){
  # if no id given, extract the yeast ORF (locus name) from the sequnece names
  if(is.null(.id)){ .id = get.strain_orf(s,what='orf') }
  #ref=data.frame(orf=character(0),pos_ref=numeric(0),aa_ref=character(0),fr_ref=numeric(0),stringsAsFactors = F)
  ref = NULL
  if( count.SNP(s,p) > 0 ){
    ref = data.frame(row.names = NULL, stringsAsFactors = F,
            orf = .id,
            pos_ref = get.SNP.positions(s,p,as.ind=T),
            aa_ref = get.REF.aa(s,p),
            fr_ref = get.REF.freq(s,p)
    )
  }
  return(ref)
}

### SNP = non-WT amino acid present in at least q=(1-p) strains
get.SNP.aa = function(s,p,as.df=T){
  FR=freq.SNP(s,p,zero2na =T)
  AA=rownames(FR)
  aa.snp = !is.na(FR) & !is.ref(s,p)
  ALT = apply(aa.snp,2,function(aa){ paste0(AA[aa],collapse="|") })
  if(as.df){
    if( length(ALT) > 0){
      alt = stack( strsplit(ALT,split="\\|") )[,c('ind','values')]
      alt$ind = unfactor(alt$ind)
      ALT = setNames(alt,c('pos_alt','aa_alt'))
      #ALT = tidyr::separate_rows( tibble::enframe(ALT), value, sep='\\|') %>% dplyr::rename(pos_alt=name, aa_alt = value)
    }else{
      ALT=data.frame(pos_alt=numeric(0),aa_alt=character(0),stringsAsFactors = F)
    }
  }
  return(ALT)
}

get.SNP.freq = function(s,p,addSNP=T){
  FR=freq.SNP(s,p,zero2na =T)
  aa.snp = !is.na(FR) & !is.ref(s,p)
  fr.snp = FR[aa.snp]
  if(addSNP){
    fr.snp = colSums(replace(FR,is.ref(s,p),NA),na.rm = T)
  }
  return(fr.snp)
}

get.SNP =function(s,p,.id=NULL){
  #  library(tidyverse)
  if(is.null(.id)){ .id = get.strain_orf(s,what='orf') }
  alt=NULL
  #alt= data.frame(stringsAsFactors = F, row.names = NULL,
  #                orf=character(0),pos_alt=numeric(0),aa_alt=character(0), fr_alt=numeric(0))
  if( count.SNP(s,p) > 0 ){
    alt = data.frame(stringsAsFactors = F, row.names = NULL,
                     orf = .id,
                     get.SNP.aa(s,p),
                     fr_alt=get.SNP.freq(s,p,addSNP=F) )
  } # %>%
  # TO CALCULATE FREQUENCIES OF JOINED VARIANTS
  # as_tibble %>%
  # group_by(pos_alt) %>%
  # mutate(fr_SNP = sum(fr_alt), single_SNP = n()== 1)
  return(alt)
}

get.STOP = function(S,pmin=0.05,as.count=T){
  M = consensusMatrix(S, as.prob = T)
  STOP = which(M["*",] > pmin & M["*",] < 1)
  earlySTOP = STOP[STOP != ncol(M)]
  if(as.count) return(length(earlySTOP))
  return(earlySTOP)
}

get.SNPscore = function(S, scoremat = "BLOSUM100", masked=c()){
  mat.sub = grep( data(package="Biostrings")$results[,'Item'], pattern="(BLOSUM)|(PAM)", value=TRUE)
  scoremat = match.arg(scoremat, choices=mat.sub, several.ok = F)
  require(msa)
  ali = msa(S)
  msaConservationScore(ali,scoremat)
}

# show.positions(orf,pos,letter){
#   pos = c(19,54,79,85)
#   ipos = IRangesList(pos)
#   extract
# }

# Get strains name sharing base/residue at position pos in fasta sequences (DNA/AA)
which.strains = function(orf,pos,letter){
  if(length(pos)>1){ warning("More than one position provided, only the first will be used!") }
  if(any(pos[1]>widths(orf))){ stop(sprintf("Position %s is out of bounds!",pos[1])) }
  matched =  names(orf)[hasLetterAt(orf,letter[1],at=pos[1])]
  matched_strains = get.strain_orf(matched,what='strains')
  return(matched_strains)
}

#http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_Current_Release.tgz
# SNPDIR="/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/Matrix/RefGenome_SNP_indel/"
# yal001c=Biostrings::readDNAStringSet(file.path(SNPDIR,"YAL001C.fasta"))
# yal001c.freq = consensusMatrix(yal001c,as.prob = T)
# noSNP  = colSums(yal001c.freq==1) == 1
# SNP.freq = apply(yal001c.freq[,!noSNP],2,max)
# SNP.05 = SNP.freq < (1-0.5)
#sum(noSNP)
#sum(SNP.05)
#colMeans( letterFrequency( x = subseq(yal001c,start=which(SNP.05)[1],width=1), letter=c("A","C","G","T")) )
