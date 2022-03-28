#source("src/utils.r",local = T)
library(stringr)
library(Biostrings)

# Sequence identifiers ---------------------------------------------------------
SGD.nomenclature = function(coding=T,rna=F){
  nuclear = "[Y][A-P][LR][0-9]{3}[WC](?:-[A-Z])?"
  mito = "Q[0-9]{4}"
  plasmid = "R[0-9]{4}[WC]"
  tRNA = "t[ATCG]\\([ATCG]{3}\\)[A-P][0-9]+"
  snRNA = "snR[0-9]+[a-z]"
  rRNA = "RDN[0-9]+\\-?[0-9]"

  ORF = paste(collapse='|',sprintf("(%s)",c(nuclear,mito,plasmid)))
  RNA = paste(collapse='|',sprintf("(%s)",c(tRNA,snRNA,rRNA)))
  if( coding & rna ){ return(sprintf("%s|%s",ORF,RNA)) }
  if( coding & !rna ){ return(ORF) }
  if( !coding & rna ){ return(RNA) }
  if( !coding & !rna ){ return(nuclear) }
}

UNIPROT.nomenclature = function(){
  ACCESSION = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
  return(ACCESSION)
}

# Sequences --------------------------------------------------------------------
load.proteome = function(url,nostop=T) {
  library(Biostrings)
  p = Biostrings::readAAStringSet(filepath = url)
  if(nostop){ # Remove the trailing star from amino acid sequence (stop codon)
    star = Biostrings::subseq(p,start=-1) == '*'
    p = Biostrings::subseq(p,start=1,  end=Biostrings::width(p)-star)
  }
  return(p)
}

load.genome = function(url) {
  library(Biostrings)
  g = Biostrings::readDNAStringSet(filepath = url)
  return(g)
}

stripR = function(BS,n=1){
  has_stop = Biostrings::subseq(BS,start=-1) != '*'
  if(all(has_stop)){
    warning("None of the sequence ends with a STOP codon (symbol *)")
  }
  p =Biostrings::subseq(BS,start=1,  end=Biostrings::width(BS)-n)
  return(p)
}

stripL = function(BS,n=1){
  p =Biostrings::subseq(BS,start=n+1,  end=Biostrings::width(BS))
  return(p)
}

rm.stop = function(BS){
  subseq(BS,start = 1,end=width(BS)-1)
}

widths = function(BS){
  if(is.list(BS)){
    return( sapply(BS,width) )
  }else if( class(BS) == "AAStringSet" ){
    return(width(BS))
  }else if( class(BS) == "AAStringSetList"){
    return( sapply(nchar(BS),unique) )
  }else{
    stop('function not defined for unknown class of input')
  }
}

get.orf.filename = function(seqfile){
  return( gsub(x=basename(seqfile),pattern='\\.(fasta)$',replacement = "", ignore.case = F, perl=T) )
}

count.fasta = function(fastafile){
  return(sum(grepl("^>",readLines(fastafile))))
}

get.orf.fasta = function(fastafile){
  library(Biostrings)
  fasta = Biostrings::readAAStringSet(fastafile)
  return(names(fasta))
}

get.width = function(BS){
  df.width = data.frame(
    orf=names(BS),
    len=widths(BS),
    stringsAsFactors = F, row.names = NULL
  )
  return(df.width)
}

get.positions = function(BS){
  df.pos = data.frame(
    orf=rep(names(BS),width(BS)),
    len=rep(widths(BS),times=widths(BS)),
    wt_pos=unlist(lapply(width(BS),seq_len)),
    wt_aa = str2chr(as.character(BS)),
    stringsAsFactors = F, row.names = NULL
  )
  return(df.pos)
}

read.sequences = function(seqfiles,strip.fname=F,ncores=parallelly::availableCores(which='max')-2,type="AA"){
  library(Biostrings)
  library(tictoc)
  library(progress)
  library(parallel)
  library(parallelly)
  #library(biomartr) # NOT library
  #library(purrr)    # NOT library
  seqtype = match.arg(type, choices=c('AA','DNA'),several.ok = F)

  task="Reading fasta sequences..."
  tic(msg = task)
  message(task)
  pb = progress::progress_bar$new(total = length(seqfiles), width = 100, clear=T,
                                  format = " (:spin) :what [:bar] :percent (:current/:total # :elapsed eta: ~:eta)")

  readSequence = function(file, .pb=NULL, .pb.toprint=task,.seqtype='AA'){
    # format = "fasta", obj.type = "Biostrings"
    if(!.pb$finished){ .pb$tick(tokens=list(what=.pb.toprint)) }
    #return( biomartr::read_proteome(file,format,obj.type,...) )
    if(seqtype == 'AA'){
      return(Biostrings::readAAStringSet(file, format = 'fasta' ))
    }else if(seqtype == 'DNA'){
      return(Biostrings::readDNAStringSet(file, format = 'fasta' ))
    }
  }
  .arg.ReadSequence =  list(.pb=pb, .pb.toprint=task,.seqtype=seqtype)
  if (require('pbmcapply')) {
    message(sprintf("using 'pbmcapply' to track progress in parallel across %s cpus",ncores))
    S = pbmcapply::pbmcmapply(FUN=readSequence,  seqfiles, MoreArgs = .arg.ReadSequence,
                              mc.cores=ncores, mc.silent=F, mc.cleanup = T)
  }else if( require('parallel') ){
    message(sprintf("using 'parallel' package across %s cpus (no progress bar)",ncores))
    S = parallel::mcmapply( FUN=readSequence,  seqfiles, MoreArgs = .arg.ReadSequence,
                            mc.cores=ncores, mc.silent=F, mc.cleanup = T)
  }else{
    message("sequential computation with progress bar")
    S = mapply(FUN=readSequence,  seqfiles, MoreArgs = .arg.ReadSequence)
  }

  if(strip.fname){
    warning("filename will be used as the proteome identifier")
    names(S) = get.orf.filename(seqfiles)
  }
  toc()
  if(seqtype == 'AA'){
    return(AAStringSetList(S))
  }else if(seqtype == 'DNA'){
    return(DNAStringSetList(S))
  }
}

seq2char = function(seq){ return( unlist(strsplit(as.character(seq),"")) ) }

seq2df = function(BS){
  is_string_set = class(BS) == 'AAStringSet'
  is_string = class(BS) == 'AAString'
  if( !is_string_set & !is_string ) { stop("requires a 'AAStringSet' or 'AAString' object!") }

  if(is_string_set & length(BS)==1){
    tibble::tibble(id = names(BS), resi = 1:Biostrings::width(BS), resn = seq2char(BS) )
  }else if(is_string_set & length(BS)>1){
    lapply(BS,seq2df) %>% dplyr::bind_rows(.id='id')
  }else if(is_string){
    tibble::tibble(resi = 1:length(BS), resn = seq2char(BS) )
  }else{
    NA
  }
}