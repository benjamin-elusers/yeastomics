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

maskseq_ =function(BS=AAString("abcdefghijkl"), S=3, E=6){
  # Switch case at specific position of the Biostrings sequence
  tomask = subseq(x=BS, start=S,end=E) %>% toggle_case() %>% Biostrings::AAString()
  subseq(x=BS, start=S,end=E) <- tomask
  return(BS)
}

maskseq =function(BS=AAString("abcdefghijkl"), S=c(3,8), E=c(6,12)){
  # Recursive version of maskseq
  if(length(S)==0 && length(E)==0){ return(BS) }
  BS = maskseq_(BS,S[1],E[1])
  return( maskseq(BS,tail(S,-1),tail(E,-1)) )
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

# widths = function(BS){
#   if(is.list(BS)){
#     W = sapply(BS, function(x){
#       if(n_distinct(width(x))==1){ return(unique(width(x))) }
#       return(width(x))
#       })
#     return(W)
#   }else if( class(BS) %in% c("AAStringSet","DNAStringSet") ){
#     return(width(BS))
#   }else if( class(BS) %in% c("AAStringSetList","DNAStringSetList") ){
#     return( unlist( lapply(nchar(BS),unique) )  )
#   }else{
#     stop(sprintf('function not defined for unknown class of input [%s]',class(BS)))
#   }
# }

all_same <- function(x) {
  length(unique(x)) == 1
}

# Function to calculate widths for different classes
calculate_widths <- function(x) {
  if (inherits(x, c("AAStringSet", "DNAStringSet"))) {
    w = width(x)
    if(all_same(w)){ return( unique(w) ) }
    return(w)
  }else if (inherits(x, c("AAStringSetList", "DNAStringSetList")) ){
    return( unlist( lapply(nchar(x),unique) )  )
  } else if (inherits(x, c("AAMultipleAlignment", "DNAMultipleAlignment"))) {
    return(unique(ncol(x)))
  }
  return(width(x))
}

widths <- function(BS) {
  if (is.list(BS)) {
    W <- sapply(BS, calculate_widths)
    return(W)
  } else {
    return(calculate_widths(BS))
  }
}

get.orf.filename = function(seqfile){
  return( gsub(x=basename(seqfile),pattern='\\.(fasta)$',replacement = "", ignore.case = F, perl=T) )
}

count.fasta = function(fastafile){
  return(sum(grepl("^>",readLines(fastafile))))
}

get.fasta.names= function(fastafile){
  library(Biostrings)
  fasta = Biostrings::readAAStringSet(fastafile, format = 'fasta')
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
  nseq =  length(seqfiles)
  task=sprintf("Reading %s fasta sequences...",nseq)
  tic(msg = task)
  message(task)
  pb = progress::progress_bar$new(total =nseq, width = 100, clear=T,
                                  format = " (:spin) :what [:bar] :percent (:current/:total # :elapsed eta: ~:eta)")

  readSequence = function(file, .pb=NULL, .pb.toprint=task,.seqtype='AA'){
    # format = "fasta", obj.type = "Biostrings"
    if(!.pb$finished){ .pb$tick(tokens=list(what=.pb.toprint)) }
    #return( biomartr::read_proteome(file,format,obj.type,...) )
    if(seqtype == 'AA'){
      return(Biostrings::readAAStringSet(file, format = 'fasta'))
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

get_codon_table = function(){
  library(seqinr)
  codon_table = seqinr::SEQINR.UTIL$CODON.AA %>% as_tibble() %>%
                mutate(CODON=toupper(CODON),
                       AA=str_to_title(AA),
                       L=str_replace(L, "\\*", "STOP"),
                       codon_aa = paste0(CODON,"_",AA,"_",L))
  return(codon_table)
}

get_aa_score = function(string,score){
  aa_scores = get.aascales()
  scores = names(aa_scores)[-1]

  if(!missing(score)){
    scores = match.arg(score,all_scores,several.ok = T)
    aa_scores = aa_scores[,c('AA',scores)]
  }

  if(string == "" | is.na(string)){
    df_scores = set_names(rep(as.numeric(NA),length(scores)),scores) %>% as_tibble_row()
    return(df_scores)
  }

  BS = Biostrings::AAString(string)
  aa_count = alphabetFrequency(BS) %>% enframe(name = 'aa','n')


  sum_scores = left_join(aa_count,aa_scores,by=c('aa'='AA')) %>%
                  pivot_longer(cols=all_of(scores),names_to='scales',values_to='score') %>%
                  mutate(tot_score = n * score ) %>%
                  group_by(scales) %>%
                  summarize(sum_score=sum_(tot_score))

  df_scores = pivot_wider(sum_scores,names_from = 'scales',values_from = 'sum_score')

  return(df_scores)
}

normalize_sequence = function(BS){
  #data("BLOSUM62")
  BS_norm = BS %>% chartr("U","C",.) %>% chartr("O","K",.) %>% chartr("J","L",.)
  return(BS_norm)
}

count_aa = function(BS){
  seqids = BS |> names()
  col_AA   = Biostrings::AA_STANDARD
  col_noAA = setdiff(Biostrings::AA_ALPHABET,col_AA)

  freqaa = Biostrings::alphabetFrequency(BS,as.prob = F) |>
    as_tibble() |>
    mutate(ids=seqids, naa=width(BS), nseq = n_distinct(seqids)) |>
    relocate(nseq,ids,naa) |>
    group_by(ids,naa,nseq) |>
    mutate( noAA = sum_( c_across(cols=all_of(col_noAA)) ),
            AA = sum_( c_across(cols=all_of(col_AA)) ),
            f_noAA = noAA / naa,
            f_AA = AA / naa
    ) |> nest(aacount=c(col_AA,"other"))
  return(freqaa)
}


get_ambiguous_codon <- function(codon, code_map=Biostrings::IUPAC_CODE_MAP) {
  # Split the codon into its bases
  bases <- strsplit(codon, "")[[1]]
  # Get the possible replacements for each base
  replacements <- lapply(bases, function(b) strsplit(code_map[b], "")[[1]])
  # Create all combinations of the replacements
  possible_codons <- expand.grid(replacements[[1]], replacements[[2]], replacements[[3]])
  # Collapse each combination into a codon string
  apply(possible_codons, 1, paste, collapse = "")
}

get_all_codons = function(shorten=F, rm.wildcard=T, rm.negative=T){

  iupac_bases = names(Biostrings::IUPAC_CODE_MAP)
  if(rm.wildcard){ iupac_bases = setdiff( iupac_bases,"N") }
  if(rm.negative){ iupac_bases = setdiff( iupac_bases,c("V","H","D","B") ) }

  codons_long = expand_grid(b1=iupac_bases,
                       b2=iupac_bases,
                       b3=iupac_bases) |>
           rowwise() |>
           mutate(fuzzy_codon = paste0(b1,b2,b3,collapse=''),
                  n_fuzzy_bases = sum(!(c(b1,b2,b3) %in% c("U",Biostrings::DNA_BASES)))) |>
           mutate(codon = list(get_ambiguous_codon(fuzzy_codon, Biostrings::IUPAC_CODE_MAP[iupac_bases]))) |>
           unnest(col=codon) |>
           mutate(aa_codon = Biostrings::GENETIC_CODE[codon]) |>
           group_by(fuzzy_codon) |>
           mutate(n_codons = n(),
                  n_aa = n_distinct(aa_codon)) |>
           arrange(n_codons,n_aa) |>
           mutate(is_ambiguous = n_codons > 1)

    codons_short = codons_long |>
    group_by(fuzzy_codon, is_ambiguous, n_fuzzy_bases, n_codons, n_aa) %>%
    summarise(
      codon_ambiguous = paste(unique(codon), collapse="/"),
      aa_ambiguous = paste(aa_codon, collapse="/"),
      .groups = "drop"
    )
    codons = codons_long
    if(shorten){ codons = codons_short }
  return(codons)
}

