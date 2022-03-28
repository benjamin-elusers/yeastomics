# FUNCTIONS ---------------------------------------------------------------
read.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con,skipNul=T)
  #closeAllConnections()
  close.connection(con)
  return(txt)
}

load.proteome = function(url,nostop=T) {
  #library(Biostrings)
  p = Biostrings::readAAStringSet(filepath = url)
  if(nostop){ # Remove the trailing star from amino acid sequence (stop codon)
    star = Biostrings::subseq(p,start=-1) == '*'
    p = Biostrings::subseq(p,start=1,  end=Biostrings::width(p)-star)
  }
  return(p)
}

find.uniprot_refprot = function(keyword,all=T,GUI=interactive()){
  #library(stringr)
  #library(readr)
  #library(janitor)
  library(magrittr) # for using pipe operator (%>%)
  UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
  URL_README = paste0(UNIPROT_URL,"knowledgebase/reference_proteomes/README")
  README = readr::read_lines(URL_README)
  row_header = stringr::str_subset(README, pattern = '^Proteome_ID\\tTax_ID\\t') %>%
    stringr::str_split('\t') %>% unlist() %>%
    stringr::str_replace(stringr::fixed('#(1)'),'n_canonical') %>%
    stringr::str_replace(stringr::fixed('#(2)'),'n_isoforms') %>%
    stringr::str_replace(stringr::fixed('#(3)'),'n_gene2acc')

  row_content = stringr::str_subset(README, pattern = "^UP[0-9]+\\t[0-9]+\\t")
  refprot = readr::read_tsv(file=I(row_content), col_names = row_header) %>%
    janitor::clean_names() %>% dplyr::arrange(tax_id)
  if(!missing(keyword)){
    matched = refprot %>% dplyr::filter(dplyr::if_any(everything(),stringr::str_detect, keyword))
    message(sprintf('%s entries matched keyword "%s"',nrow(matched),keyword))
    return(matched)
  }else if(!all){
    name = sprintf("%s (taxid %s)",refprot$species_name,refprot$tax_id)
    which_prot = menu(name, graphics=GUI,title = 'pick an organism below...(sorted by tax_id)')
    return(refprot[which_prot,])
  }else{
    return(refprot)
  }
}

get.uniprot.proteome = function(taxid,DNA=F,fulldesc=F) {
  #library(stringr)
  if(missing(taxid)){  stop("Need an uniprot taxon id") }
  UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
  SEQTYPE = ".fasta.gz"

  refprot = find.uniprot_refprot(all=T)
  found = refprot$tax_id %in% taxid
  if(!any(found)){ stop(sprintf("%s not found in the reference proteome!",taxid)) }
  TAX = stringr::str_to_title(refprot$superregnum[which(found)])
  UPID = refprot$proteome_id[which(found)]
  if(DNA){ seqtype = "_DNA.fasta.gz" }
  proteome_url = sprintf("%s/%s/%s/%s_%s%s",UNIPROT_URL,TAX,UPID,UPID,taxid,SEQTYPE)

  UNI = load.proteome(proteome_url)
  if(!fulldesc){
    regexUNIPROTAC = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
    names(UNI) = stringr::str_extract(names(UNI), regexUNIPROTAC)
  }
  return(UNI)
}


load.pfam = function(tax='559292'){
  library(magrittr) # for using pipe operator (%>%)
  #library(readr)
  #library(tidyr)
  #library(stringr)
  #library(dplyr)
  # Load domains assignments based of HMM profiles (PFAM)
  # 559292 = S.cerevisiae
  URL_PFAM = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
  .release = read.url(paste0(URL_PFAM,"/Pfam.version.gz")) %>% paste0(sep="\n")
  message("CURRENT RELEASE")
  message("---------------")
  message(.release)
  message("---------------")
  url.pfam = sprintf("%s/proteomes/%s.tsv.gz",URL_PFAM,tax)
  url.clan = sprintf("%s/Pfam-A.clans.tsv.gz",URL_PFAM)
  clans = readr::read_delim(url.clan,delim="\t",col_names=c("pfam_id",'clan_id','clan_name',"pfam_name","pfam_desc")) %>%
    dplyr::mutate( clan_id = tidyr::replace_na(clan_id,replace = 'No_clan') )

  .vers = read.url(url.pfam)[1]
  message(.vers)
  .nprot = read.url(url.pfam)[2]
  message(.nprot)
  header = read.url(url.pfam)[3]
  # parse the header (3rd commented rows)
  columns = header %>%
    stringr::str_split(pattern = "> <") %>% unlist %>%
    stringr::str_replace_all("[ \\-]","_") %>% stringr::str_remove_all("[#<>]")
  pfam = readr::read_delim(file = url.pfam, skip=2,comment = "#",delim="\t", col_names = columns,escape_double = F,guess_max = 100) %>%
    dplyr::left_join(clans,by = c('clan'="clan_id",'hmm_acc'='pfam_id','hmm_name'='pfam_name') )

  return(pfam)
}

get.superfamily.species = function(){

  URL_SUPERFAMILY = "https://supfam.org/SUPERFAMILY"
  url_gen_list = paste0(URL_SUPERFAMILY,"/cgi-bin/gen_list.cgi")
  superfamily_gen_list  = rvest::read_html(url_gen_list)
  supfam_taxlevels = superfamily_gen_list %>%
    rvest::html_elements(xpath = '//table/preceding-sibling::strong') %>%
    rvest::html_text()

  # Retrieve hyperlinks on genome names (contain abbreviaiotns in href)
  genomes_abbr = superfamily_gen_list %>%
    rvest::html_elements("a[href*='gen_list.cgi?genome=']") %>%
    rvest::html_attr('href') %>%
    stringr::str_replace(stringr::fixed("gen_list.cgi?genome="),"")

  get_genome_info_taxon_id = function(x){
    tryCatch(
      rvest::read_html(x) %>%
        rvest::html_elements("table") %>%
        .[[4]] %>%
        rvest::html_elements("td") %>%
        rvest::html_text() %>%
        .[ which(stringr::str_detect(string=.,pattern="NCBI Taxon ID:")) + 1 ]
      #finally=print(paste0("get ncbi taxon id for: ",g))
    )
  }
  url_genomes = paste0(URL_SUPERFAMILY,"/cgi-bin/info.cgi?genome=",genomes_abbr)

  tictoc::tic()
  message("Retrieving ncbi taxon id from genome information...")
  if (require('pbmcapply')) {
    ncores=parallelly::availableCores(which='max')-1
    message(sprintf("using 'pbmcapply' in parallel with %s cpus (~3mn on 10 cpus)",ncores))
    genome_ncbi_taxid = pbmcapply::pbmcmapply(FUN=get_genome_info_taxon_id,  url_genomes, mc.cores=ncores, mc.silent=F, mc.cleanup = T)
  }else{
    message("NOT PARALLEL! This may take a while (>20mn)")
    genome_ncbi_taxid = furrr::future_map_chr(url_genomes,get_genome_info_taxon_id)
  }
  tictoc::toc()

  supfam_genomes = superfamily_gen_list %>%
    rvest::html_elements("table.small_table_text") %>%
    rvest::html_table() %>%
    stats::setNames(janitor::make_clean_names(supfam_taxlevels)) %>%
    dplyr::bind_rows(.id = 'taxlevel') %>%
    janitor::clean_names() %>%
    dplyr::mutate(genome=genomes_abbr,ncbi_taxid=genome_ncbi_taxid) %>%
    dplyr::relocate(taxlevel,ncbi_taxid,genome)

  return(supfam_genomes)
}

load.superfamily = function(tax='xs'){
  # Load domain assignments from superfamily SCOP level (SUPFAM.org)
  # xs = saccharomyces cerevisiae
  #library(httr)
  #library(readr)
  #library(janitor)
  URL_SUPERFAMILY = "https://supfam.org/SUPERFAMILY"
  download.request =  httr::GET(sprintf("%s/cgi-bin/save.cgi?var=%s;type=ass",URL_SUPERFAMILY,tax))
  superfamily.txt = httr::content(download.request,as = 'text')
  #superfamily.assignment = unlist(stringi::stri_split_lines(superfamily.txt,omit_empty = T))[-c(1:2)]
  supfam = readr::read_delim(file=superfamily.txt,comment = "#",delim="\t", escape_double = F) %>% janitor::clean_names()
  return(supfam)
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
# MAIN --------------------------------------------------------------------
id = id_pfam_uni[8664]
uniprot = hs_uni[ id ]


assign_pfam_uniprot = function(taxid=9606){
  find.uniprot_refprot(keyword='homo',all=F)
  uni_seq = get.uniprot.proteome(taxid)
  pfam_data = load.pfam(taxid)
  id_pfam_uni = intersect(pfam_data$seq_id,names(uni_seq)) %>% unique
  cat(sprintf("%s"))

  df_uni = seq2df(hs_uni[uniprot])
  pfam_uni = hs_pfam %>% dplyr::filter(seq_id == id)

  pfam = tibble()
  for(d in 1:nrow(pfam_uni)){
    tmp = pfam_uni[d,]
    l=tmp$hmm_length-1
    pfam_dom =
      tibble( pos_pfam = tmp$alignment_start:tmp$alignment_end,
            pfam_id = tmp$hmm_acc,
            pfam_type = tmp$type,
            pfam_name = tmp$hmm_name,
            pfam_score = tmp$bit_score,
            pfam_Evalue = tmp$E_value,
            pfam_clan = tmp$clan,
            pfam_clan_name = tmp$clan_name,
            pfam_desc = tmp$pfam_desc
      )
    pfam = bind_rows(pfam,pfam_dom)
  }



}
