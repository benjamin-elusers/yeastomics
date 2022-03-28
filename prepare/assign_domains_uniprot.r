options(dplyr.width=Inf)

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

get.uniprot.mapping = function(taxid) {
  if(missing(taxid)){  stop("Need an uniprot taxon id") }
  UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
  EXTENSION = ".idmapping.gz"
  refprot = find.uniprot_refprot(all=T)
  found = refprot$tax_id %in% taxid
  if(!any(found)){ stop(sprintf("%s not found in the reference proteome!",taxid)) }
  TAX = stringr::str_to_title(refprot$superregnum[which(found)])
  UPID = refprot$proteome_id[which(found)]

  gene2acc_url = sprintf("%s/%s/%s/%s_%s%s",UNIPROT_URL,TAX,UPID,UPID,taxid,EXTENSION)
  mapped = readr::read_delim(gene2acc_url,delim='\t',col_names=c('uni','extdb','extid')) %>%
    dplyr::mutate(sp=taxid,upid=UPID)
  return(mapped)
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
    matched = refprot %>% dplyr::filter(dplyr::if_any(everything(),stringr::str_detect, as.character(keyword)))
    message(sprintf('%s entries matched keyword "%s"',nrow(matched),keyword))
    if(all || nrow(matched)==1){
      return(matched)
    }else{
      species = sprintf("%s (%s)",matched$species_name,matched$tax_id)
      selection = menu(species,title = 'pick a species below...')
      return(matched[selection,])
    }
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
  library(magrittr)
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
    taxonid = tryCatch(
      rvest::read_html(x) %>%
        rvest::html_elements("table") %>%
        .[[4]] %>%
        rvest::html_elements("td") %>%
        rvest::html_text() %>%
        .[ which(stringr::str_detect(string=.,pattern="NCBI Taxon ID:")) + 1 ]
      #finally=print(paste0("get ncbi taxon id for: ",g))
    )
    if(length(taxonid)){ return(taxonid) }
    return(NA)
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

  #tibble( taxid = unlist(genome_ncbi_taxid), genome=genomes_abbr )
  supfam_genomes = superfamily_gen_list %>%
    rvest::html_elements("table.small_table_text") %>%
    rvest::html_table() %>%
    stats::setNames(janitor::make_clean_names(supfam_taxlevels)) %>%
    dplyr::bind_rows(.id = 'taxlevel') %>%
    janitor::clean_names() %>%
    dplyr::mutate(genome=genomes_abbr,ncbi_taxid=unlist(genome_ncbi_taxid)) %>%
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

assign_pfam_uniprot = function(taxid=9606){
  # get uniprot proteome to dataframe
  uni_seq = get.uniprot.proteome(taxid)
  df_uni = seq2df(uni_seq)

  # get pfam data for reference proteome identifiers
  pfam_data = load.pfam(taxid)
  pfam_uni_data = pfam_data %>% dplyr::filter( seq_id %in% names(uni_seq))
  message(sprintf("filtered Pfam (%s rows) for uniprot reference proteome [%s]",nrow(pfam_uni_data),taxid))

  # extend pfam data to residue level
  pfam_res = pfam_uni_data %>%
    dplyr::group_by(seq_id,alignment_start,alignment_end) %>%
    dplyr::mutate(pfam_pos=min(alignment_start)) %>%
    tidyr::complete(pfam_pos = seq(alignment_start, alignment_end, by = 1)) %>%
    tidyr::fill(everything(),.direction = "down")

  # merge uniprot reference proteome and pfam data at residue level
  uni_pfam_assigned = dplyr::left_join(df_uni,pfam_res, by=c('id'='seq_id','resi'='pfam_pos'))
  return(uni_pfam_assigned)
}

assign_superfamily_uniprot = function(taxid="9606",supfam_sp){

  uni_tax = find.uniprot_refprot(taxid,all = F)
  uni2acc = get.uniprot.mapping(uni_tax$tax_id)
  uni_seq = get.uniprot.proteome(uni_tax$tax_id)

  # superfamily uses their own abbrevation to designate a genome (takes 5mn to be matched to NCBI taxid)
  is_number_taxid = grepl("^[0-9]+$",taxid)
  if(is_number_taxid){
    if(missing(supfam_sp)){ supfam_sp=get.superfamily.species() }
    sp_taxid = supfam_sp %>% dplyr::filter(taxid == ncbi_taxid)
    if(nrow(sp_taxid)>1){
      which_sp = sprintf('%s -> %s (%s=%s)',sp_taxid$taxlevel,sp_taxid$genome_name,sp_taxid$ncbi_taxid,sp_taxid$genome)
      selection=menu(which_sp,title = 'pick one species below...')
      tax = sp_taxid$genome[selection]
    }else if(nrow(sp_taxid)==1){
      tax = sp_taxid$genome
    }else{
      stop(sprintf('no superfamily data for taxon %s!',taxid))
    }
  }else{
    matched = supfam_sp %>% dplyr::filter(dplyr::if_any(everything(),stringr::str_detect, taxid))
    message(sprintf('%s entries matched keyword "%s"',nrow(matched),keyword))
    if( nrow(matched)>1 ){
      which_sp = sprintf('%s -> %s (%s=%s)',matched$taxlevel,matched$genome_name,matched$ncbi_taxid,matched$genome)
      selection=menu(which_sp,title = 'pick one species below...')
      tax = matched$genome[selection]
    }else if(nrow(matched)==1){
      tax = matched$genome
    }else{
      stop(sprintf('no superfamily data for keyword "%s"!',taxid))
    }
  }

  # get superfamily data for reference proteome identifiers
  # superfamily does not use uniprot (e.g. human is ensembl)
  supfam_data = load.superfamily(tax) %>%
    tidyr::separate_rows(region_of_assignment,sep = ",") %>%
    tidyr::separate(col=region_of_assignment, into = c('scop_start','scop_end'), sep = '\\-',convert = T) %>%
    tidyr::separate(col=sequence_id, into = c('db','id'), sep = '\\|',remove=F)

  if(any(uni2acc$extid %in% supfam_data$sequence_id)){
    id_mapped_uni = uni2acc %>% dplyr::filter(extid %in% supfam_data$sequence_id)
    supfam_uni_data =  dplyr::inner_join(supfam_data,id_mapped_uni,by=c('sequence_id'='extid'))

  }else{ # superfamily sequence id is not identical to uniprot mapped identifier
    id_mapped_uni = uni2acc %>% dplyr::filter(extid %in% supfam_data$id)
    supfam_uni_data =  dplyr::inner_join(supfam_data,id_mapped_uni,by=c('id'='extid'))
  }
  uniref_matched =supfam_uni_data$uni %in% names(uni_seq)
  message(sprintf("filtered Pfam (%s rows) for uniprot reference proteome [%s]",sum(uniref_matched),taxid))

  # extend superamily data to residue level
  supfam_res = supfam_uni_data %>%
    dplyr::group_by(sequence_id,uni,scop_start,scop_end) %>%
    dplyr::mutate(supfam_pos=min(scop_start)) %>%
    tidyr::complete(supfam_pos = seq(min(scop_start), max(scop_end), by = 1)) %>%
    tidyr::fill(everything(),.direction = "down")

  # get uniprot proteome to dataframe
  df_uni = seq2df(uni_seq)

  # merge uniprot reference proteome and pfam data at residue level
  uni_supfam_assigned = dplyr::left_join(df_uni,supfam_res, by=c('id'='uni','resi'='supfam_pos'))
  return(uni_supfam_assigned)
}

# MAIN --------------------------------------------------------------------
supfam_genomes= get.superfamily.species()

hs_pfam_uni = assign_pfam_uniprot(9606)
sc_pfam_uni = assign_pfam_uniprot(559292)
ec_pfam_uni = assign_pfam_uniprot(83333)


# Must validate selection
hs_supfam_uni = assign_superfamily_uniprot(9606,supfam_genomes)
sc_supfam_uni = assign_superfamily_uniprot('cerevisiae',supfam_genomes) # for S288c select 5 and then 1
ec_supfam_uni = assign_superfamily_uniprot('coli',supfam_genomes) # for K12 select 10 and then 63


