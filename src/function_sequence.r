#source("src/utils.r",local = T)
library(stringr)
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

UNIPROT.nomenclature = function(){
  ACCESSION = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
  return(ACCESSION)
}

fallback = function(url,archived){
  # if current download is dead/blocked use the latest archive
  if( httr::http_error(url) ){
    warning(sprintf("falling back to last release archived (%s)...",archived),immediate. = T)
    url = archived
  }else{
    message("URL...OK")
  }
  return(url)
}

load.sgd.CDS = function(withORF=T,orf.dna="sequence/S288C_reference/orf_dna") {
  library(stringr)
  sgd.url = "http://sgd-archive.yeastgenome.org"
  cds= file.path(sgd.url,orf.dna,"orf_coding_all.fasta.gz")
  cds_archived = file.path(dirname(cds),"archive/orf_coding_all_R64-3-1_20210421.fasta.gz")
  SGD = load.genome(fallback(cds,cds_archived))
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

load.sgd.proteome = function(withORF=T,rm.stop=T, orf_protein="sequence/S288C_reference/orf_protein") {
  library(stringr)
  sgd.url = "http://sgd-archive.yeastgenome.org"
  prot= file.path(sgd.url,orf_protein,"orf_trans_all.fasta.gz")
  prot_archived = file.path(dirname(prot),"archive/orf_trans_all_R64-3-1_20210421.fasta.gz")
  SGD = load.proteome(fallback(prot,prot_archived),nostop = rm.stop)
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

load.pombase.proteome = function(withORF=T,rm.version=T) {
  library(stringr)
  pombase.url = "ftp://ftp.pombase.org/pombe/genome_sequence_and_features/feature_sequences/peptide.fa.gz"
  Pombase = load.proteome(pombase.url)
  regexPombaseID = "(SP[^ ]+)(?=:pep)"
  regexPombase = "(?<=:pep )(.+)(?=\\|)"

  orf = str_extract(names(Pombase), regexPombaseID) # ORF identifier
  if(rm.version){ orf = str_remove(orf, "(\\.[0-9]+$)") }

  gname = str_extract(names(Pombase), regexPombase) # Pombase standard name
  has_noname = is.na(gname)
  gname[has_noname] = orf[has_noname]

  if(withORF){
    names(Pombase) = orf
  }else{
    names(Pombase) = gname
  }

  return(Pombase)
}

find.uniprot_refprot = function(search,all=T){
  library(stringr)
  library(readr)
  UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
  URL_README = paste0(UNIPROT_URL,"knowledgebase/reference_proteomes/README")
  README = readr::read_lines(URL_README)
  row_header = str_subset(README, pattern = '^Proteome_ID\\tTax_ID\\t') %>%
               str_split('\t') %>% unlist()

  row_content = str_subset(README, pattern = "^UP[0-9]+\\t[0-9]+\\t")
  refprot = readr::read_tsv(file=I(row_content), col_names = row_header) %>%
            janitor::clean_names()
  if(!missing(search) & !all){
    matched = refprot %>% dplyr::filter_all(any_vars(str_detect(., search)))
    return(matched)
  }else if(!all){
    name = sprintf("(%s) %s",refprot$species_name,refprot$tax_id)
    which_prot = menu(name, graphics=interactive())
    return(refprot[which_prot,])
  }else{
    return(refprot)
  }
}

get.uniprot.proteome = function(taxid,DNA=F) {

  if(missing(taxid)){  stop("Need an uniprot taxon id") }
  UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
  SEQTYPE = ".fasta.gz"

  refprot = find.uniprot_refprot(all=T)
  found = refprot$tax_id %in% taxid
  if(!any(found)){ stop(sprintf("%s not found in the reference proteome!",taxid)) }
  TAX = str_to_title(refprot$superregnum[which(found)])
  UPID = refprot$proteome_id[which(found)]
  if(DNA){ seqtype = "_DNA.fasta.gz" }
  proteome_url = sprintf("%s/%s/%s/%s_%s",UNIPROT_URL,TAX,UPID,UPID,taxid,SEQTYPE)

  UNI = load.proteome(proteome_url)
  regexUNIPROTAC = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
  names(UNI) = str_extract(names(UNI), regexUNIPROTAC)
  return(UNI)
}

load.uniprot.proteome = function(species='yeast') {
  library(stringr)
  UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
  eukaryotes = sprintf("%s/knowledgebase/reference_proteomes/Eukaryota",UNIPROT_URL)

  ## CHANGE ON FEB 2021 - Added a subdirectory per each proteome
  taxon=match.arg(species, choices = c('yeast','human'), several.ok = F)
  proteomes=c(human="UP000005640_9606.fasta.gz",yeast="UP000002311_559292.fasta.gz")
  UP=word(proteomes[taxon],1,sep = "_")
  uniprot.url = sprintf("%s/%s/%s",eukaryotes,UP,proteomes[taxon])

  UNI = load.proteome(uniprot.url)
  regexUNIPROTAC = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
  names(UNI) = str_extract(names(UNI), regexUNIPROTAC)
  return(UNI)
}

read.proteomes = function(seqfiles,strip.fname=F,ncores=parallelly::availableCores(which='max')-2){
  library(Biostrings)
  library(tictoc)
  library(progress)
  library(parallel)
  library(parallelly)
  #library(biomartr) # NOT library
  #library(purrr)    # NOT library

  task="Reading proteomes from fasta sequences..."
  tic(msg = task)
  message(task)
  pb = progress::progress_bar$new(total = length(seqfiles), width = 100, clear=T,
                         format = " (:spin) :what [:bar] :percent (:current/:total # :elapsed eta: ~:eta)")

  readProteome = function(file, .pb=NULL, .pb.toprint=task){
    # format = "fasta", obj.type = "Biostrings"
    if(!.pb$finished){ .pb$tick(tokens=list(what=.pb.toprint)) }
    #return( biomartr::read_proteome(file,format,obj.type,...) )
    return(Biostrings::readAAStringSet(file, format = 'fasta' ))
  }

  if (require('pbmcapply')) {
    message(sprintf("using 'pbmcapply' to track progress in parallel across %s cpus",ncores))
    P = pbmcapply::pbmcmapply(FUN=readProteome,  seqfiles, MoreArgs = list(.pb=pb, .pb.toprint=task),
                          mc.cores=ncores, mc.silent=F, mc.cleanup = T)
  }else if( require('parallel') ){
    message(sprintf("using 'parallel' package across %s cpus (no progress bar)",ncores))
    P = parallel::mcmapply( FUN=readProteome,  seqfiles, MoreArgs = list(.pb=pb, .pb.toprint=task),
                            mc.cores=ncores, mc.silent=F, mc.cleanup = T)
  }else{
    message("sequential computation with progress bar")
    P = mapply(FUN=readProteome,  seqfiles, MoreArgs = list(.pb=pb, .pb.toprint=task))
  }

  if(strip.fname){
    warning("filename will be used as the proteome identifier")
    names(P) = get.orf.filename(seqfiles)
  }
  toc()
  return(AAStringSetList(P))
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
  }else if( class(BS) == "AAStringSet" ){
    return(width(BS))
  }else if( class(BS) == "AAStringSetList"){
    return( sapply(nchar(BS),unique) )
  }else{
    stop('function not defined for unknown class of input')
  }
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

rm.stop = function(BS){
  subseq(BS,start = 1,end=width(BS)-1)
}
