source("src/utils.r",local = T)
### Playing with SGD (Saccharomyces Genome Database) ---------------------------
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


load.sgd.CDS = function(withORF=T) {
  library(stringr)
  sgd.url = "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz"
  SGD = load.genome(sgd.url)
  regexSGD = "(S[0-9]{9})"
  if(withORF){
    # ORF identifier
    names(SGD) = str_extract(names(SGD), SGD.nomenclature() )
  }else{
    # SGD ID
    names(SGD) = str_extract(names(SGD), regexSGD)
  }
  #names(SGD) = subname(names(SGD),sep=" ",lc=F)
  return(SGD)
}

load.sgd.proteome = function(withORF=T) {
  library(stringr)
  sgd.url = "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz"
  SGD = load.proteome(sgd.url)
  regexSGD = "(S[0-9]{9})"
  if(withORF){
    # ORF identifier
    names(SGD) = str_extract(names(SGD), SGD.nomenclature() )
  }else{
    # SGD ID
    names(SGD) = str_extract(names(SGD), regexSGD)
  }
  #names(SGD) = subname(names(SGD),sep=" ",lc=F)
  return(SGD)
}

load.pombase.proteome = function(withORF=T) {
  library(stringr)
  pombase.url = "ftp://ftp.pombase.org/pombe/genome_sequence_and_features/feature_sequences/peptide.fa.gz"
  Pombase = load.proteome(pombase.url)
  regexPombaseID = "(SP[^ ]+)(?=:pep)"
  regexPombase = "(?<=:pep )(.+)(?=\\|)"
  if(withORF){
    # ORF identifier
    names(Pombase) = str_extract(names(Pombase), regexPombaseID )
  }else{
    # Pombase standard name
    names(Pombase) = str_extract(names(Pombase), regexPombase)
  }
  return(Pombase)
}

load.uniprot.proteome = function(species='yeast') {
  library(stringr)
  uniprot.url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota"
  taxon=match.arg(species, choices = c('yeast','human'), several.ok = F)
  proteomes=c(human="UP000005640_9606.fasta.gz",yeast="UP000002311_559292.fasta.gz")
  uniprot.url = sprintf("%s/%s",uniprot.url,proteomes[taxon])

  UNI = load.proteome(uniprot.url)
  regexUNIPROTAC = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
  names(UNI) = str_extract(names(UNI), regexUNIPROTAC)
  return(UNI)
}

read.proteomes = function(seqfiles,strip.fname=F){
  library(Biostrings)
  library(tictoc)
  library(progress)
  #library(biomartr) # NOT libraryD
  #library(purrr)    # NOT libraryD

  task="Reading proteomes from fasta sequences..."
  tic(msg = task)
  message(task)
  pb = progress::progress_bar$new(total = length(seqfiles), width = 100, clear=T,
                         format = " (:spin) :what [:bar] :percent (:current/:total # :elapsed eta: ~:eta)")

  readProteome = function(file, format = "fasta", obj.type = "Biostrings", .pb=NULL, ...){
    if(!.pb$finished){ .pb$tick(tokens=list(what=task)) }
    #return( biomartr::read_proteome(file,format,obj.type,...) )
    return(Biostrings::readAAStringSet(file, format = 'fasta' ))
  }
  P = mapply( seqfiles, FUN=readProteome,  MoreArgs = list(.pb=pb))
  if(strip.fname){
    warning("filename will be used as the proteome identifier")
    names(P) = get.orf.filename(seqfiles)
  }
  toc()
  return(P)
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

widths = function(BS){
  if(is.list(BS)){
    return( sapply(BS,width) )
  }else{
    width(BS)
  }
}
