load.proteome = function(url,nostop=T) {
  library(Biostrings)
  p = Biostrings::readAAStringSet(filepath = url)
  if(nostop){ # Remove the trailing star from amino acid sequence (stop codon)
    star = Biostrings::subseq(p,start=-1) == '*'
    p = Biostrings::subseq(p,start=1,  end=Biostrings::width(p)-star)
  }
  return(p)
}

load.uniprot.proteome = function(species='yeast') {
  library(stringr)
  uniprot.url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota"
  ## CHANGE ON FEB 2021 - Added a subdirectory per each proteome
  taxon=match.arg(species, choices = c('yeast','human'), several.ok = F)
  proteomes=c(human="UP000005640_9606.fasta.gz",yeast="UP000002311_559292.fasta.gz")
  UP=word(proteomes[taxon],1,sep = "_")
  uniprot.url = sprintf("%s/%s/%s",uniprot.url,UP,proteomes[taxon])

  UNI = load.proteome(uniprot.url)
  regexUNIPROTAC = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
  names(UNI) = str_extract(names(UNI), regexUNIPROTAC)
  return(UNI)
}


# Sc=load.uniprot.proteome('yeast')
# Hs=load.uniprot.proteome('human')
#
# Biostrings::writeXStringSet(filepath = '~/Desktop/Scer-proteome.fasta', x = Sc, format = 'fasta')
# Biostrings::writeXStringSet(filepath = '~/Desktop/Hsap-proteome.fasta', x = Hs, format = 'fasta')
