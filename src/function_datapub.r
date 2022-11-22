#source("src/utils.r",local = T)
#source("src/function_sequence.r",local = T)
#source('src/function_alignment.r',local = T)
#source('src/function_phylogenetic.r',local = T)
yeastomics_url = "https://raw.githubusercontent.com/benjamin-elusers/yeastomics"
source(file.path(yeastomics_url,"main/src/utils.r"))
.dbg = log::Logger$new("DEBUG")$
  date()$
  time()$
  hook(crayon::bgWhite)

# TO REMOVE DEPENDENCIES, THOSE FUNCTIONS WERE COPIED FROM OTHER SCRIPTS
# Utils --------------------------------------------------------------------
library(openxlsx)
# open.url <- function(file_url) {
#   con <- gzcon(url(file_url))
#   txt <- readLines(con,skipNul=T)
#   #closeAllConnections()
#   close.connection(con)
#   return(textConnection(txt))
# }
#
# read.url <- function(file_url) {
#   con <- gzcon(url(file_url))
#   txt <- readLines(con,skipNul=T)
#   #closeAllConnections()
#   close.connection(con)
#   return(txt)
# }

# strfind = function(strings, patterns){
#   # Find several patterns in set of strings
#   sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = T) })
# }


# # geometric mean and standard deviation
# geomean = function(x) {  exp(mean(log(x[x != 0 & !is.na(x)]))) }
# geosd = function(x) {  exp(sd(log(x[x != 0 & !is.na(x)]))) }

# SGD ORF regular expression
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

# Uniprot Accession regular expression
UNIPROT.nomenclature = function(){
  ACCESSION = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
  return(ACCESSION)
}

# Ensembl protein regular expression
ENSEMBL.nomenclature = function(){
  ACCESSION = "(ENSP[0-9]+)"
  return(ACCESSION)
}

clean_header = function(header){
  header %>%
    tolower %>% # CONVERT TO LOWERCASE
    gsub("[[:punct:]]","",x=.) %>% # REMOVE PUNCTUATIONS%>%
    gsub("\\s+",".",x=.) # SPACES TO UNDERSCORE
}

### END OF DEPENDENCIES ###

# Remote published proteome data -----------------------------------------------

load.byrne2005.data = function() {
  message("REF: Byrne and Wolfe, 2005, Genome Research")
  message("The Yeast Gene Order Browser: Combining curated homology and syntenic context reveals gene fate in polyploid species")
  # http://ygob.ucd.ie/
  ygob.ohno = read.csv("http://ygob.ucd.ie/browser/ohnologs.csv", # s. cerevisiae ohnologs
                       stringsAsFactors = F)[, -1]
  colnames(ygob.ohno) = c('WGD_anc', 'orf1', 'gname1', 'orf2', 'gname2', 'pid', 'rlen')
  regexORF = "^[Y][A-P][LR][0-9]{3}[WC](?:-[A-Z])?$"
  isORF1 = grepl(regexORF, ygob.ohno$orf1)
  isORF2 = grepl(regexORF, ygob.ohno$orf2)
  regexYGOB = "^(Scer_YGOB_)(Y\\w+|[^Y]\\w+)$"
  # Scer_YGOB_SDC25 = YLL016W (SDC25)
  # Scer_YGOB_YDR134C= YDR134C (CCW22)
  fake = data.frame(
    old.orf = c("Scer_YGOB_SDC25", "Scer_YGOB_YDR134C"),
    orf = c("YLL016W", "YDR134C"),
    gname = c("SDC25", "CCW22"),
    stringsAsFactors = F
  )
  ygob  = ygob.ohno %>%
    mutate(
      fake.orf = orf1 %in% fake$old.orf,
      gname1 = replace(gname1, fake.orf, fake$gname),
      orf1 = replace(orf1, fake.orf, fake$orf),
      pid = as.numeric(gsub("^(\\d+)\\%", "\\1", pid)),
      ref = 1
    ) %>%
    dplyr::rename(orf = "orf1", gname = "gname1", dup.orf = "orf2", dup.gname = "gname2" ) %>%
    dplyr::select(WGD_anc, orf, gname, ref, dup.orf, dup.gname, pid, rlen)
  return(ygob)
}

load.belle2006.data = function(){
  # Load protein half-lives
  message("REF: A. Belle, A. Tanay, L. Bitincka, R. Shamir, E.K. O’Shea, 2006, PNAS")
  message("Quantification of protein half-lives in the budding yeast proteome")
  #https://doi.org/10.1073/pnas.0605420103

  #EDIT URL NOT WORKING
  # OLD "https://www.pnas.org/highwire/filestream/591690/field_highwire_adjunct_files/0/SuppDataSet.txt"
  # NEW "https://www.pnas.org/doi/suppl/10.1073/pnas.0605420103/suppl_file/suppdataset.txt"
  # THIS IS NEW ADDRESS REDIRECTS TO ANOTHER LINK

  Sdat1.url = "https://www.pnas.org/action/downloadSupplement?doi=10.1073%2Fpnas.0605420103&file=suppdataset.txt"
  #download.file(Sdat1.url, destfile = "SuppDataSet.txt")
  halflives = read.delim(Sdat1.url,
                         header = T, skip = 10,
                         stringsAsFactors = F,
                         strip.white = T, blank.lines.skip = T
  )
  halflives$X=NULL # LAST COLUMN IS EMPTY
  colnames(halflives)=c('ORF','GENE','halflife_mn','halflife_mn_corrected')
  return(halflives)
}

load.pu2008.data = function(){
  # Load macromolecular protein complexes from yeast
  message("REF: S. Pu et al., 2008, Nucleic Acids Research")
  message("Up-to-date catalogues of yeast protein complexes")

  CYC2008.url = "http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab"
  CYC2008 = read.delim(file = CYC2008.url, sep='\t', stringsAsFactors = F,header = T, na.strings = "") %>%
    as_tibble() %>%
    group_by(Complex) %>%
    fill(PubMed_id:Jaccard_Index, .direction = "down") %>%
    add_count(name="n_members")
  return(CYC2008)
}

load.costanzo2010.data = function(){
  # load gene classification of biological functions
  message("REF: M. Costanzo, A. Baryshnikova et al., 2010, Science")
  message("The Genetic Landscape of a Cell")
  # doi: https://doi.org/10.1126/science.1180823.

  S6="https://boonelab.ccbr.utoronto.ca/supplement/costanzo2009/bioprocess_annotations_costanzo2009.xls"
  #costanzo = readxl::read_xls(path=S6, sheet = 1)
  todownload="data/bioprocess_annotations_costanzo2009.xls"
  downloaded=download.file(S6, destfile = todownload)
  if(!downloaded){
    costanzo = readxl::read_excel(path = todownload, sheet = 1, progress = T,
                                  col_names = c("ORF","GENE","bioprocess"),na = c("unknown")) %>%
      separate_rows(bioprocess,sep=";") %>%
      mutate(BIOPROCESS=str_trim(bioprocess)) %>%  dplyr::select(-bioprocess)
    #         FN=factor(FUNCTION,levels=unique(FUNCTION),labels = invert(.class_function))
    # ) %>% dplyr::select(-Function) %>%
    #   group_by(ORF) %>% mutate(FN_all=paste0(unique(FN),collapse = ""))
    file.remove(todownload)
    return(costanzo)
  }else{
    stop("Could not download the supplementary table!")
    return(NA)
  }
}

load.barton2010.data = function(by=c("aa","prot")){
  # load amino acid biosynthetic cost in yeast
  message("REF: Michael D. Barton et al., PLOS One, 2010")
  message("Evolutionary Systems Biology of Amino Acid Biosynthetic Cost in Yeast")
  #https://doi.org/10.1371/journal.pone.0011935
  S1.url = "https://doi.org/10.1371/journal.pone.0011935.s001" # AA biosynthetic cost
  S4.url = "https://doi.org/10.1371/journal.pone.0011935.s004" # Protein biosynthetic cost
  BY = match.arg(by,choices=c('aa','prot'),several.ok = F)
  # Trasncript levels are obtained from 3 dilutions and 4 sources of nutrients published in:
  # JI Castrillo et al. 2007 -> "Growth control of the eukaryote cell: a systems biology study in yeast" (Journal of Biology)
  library(seqinr)
  todownload=tempfile()
  if(BY=="aa"){
    aa.df = data.frame( aa=seqinr::a(), aaa=tolower(seqinr::aaa()))
    downloaded=download.file(S1.url, destfile = todownload)
    biosynth_cost = readr::read_csv(todownload) %>%
      left_join( aa.df, by=c('amino_acid'='aaa')) %>%
      relocate(aa)
  }else if(BY=="prot"){
    downloaded=download.file(S4.url, destfile = todownload)
    raw_cost=readr::read_csv(todownload,col_types = list(X1="_"))
    # raw_cost %>% arrange(name,environment,cost_type) %>% print(n=100)
    # trna.var = raw_cost %>% group_by(name,dilution) %>% summarise(vtrna=var(transcript_level))
    biosynth_cost = raw_cost %>%
      pivot_wider(id_cols = c(name,trna,cai,transcript_level),
                  names_from=cost_type,
                  names_glue = "{.value}_{cost_type}",
                  values_from=cost, values_fn = unique) %>%
      rename_with(.cols=contains("-"),.fn=str_replace_all,pattern="-",replacement="_")

  }
  return(biosynth_cost)
}

load.marguerat2012.data = function(raw=F){
  # load protein abundance for fission yeast (pombe)
  # !! WARNING !! paxDB has values which do not correspond to original data
  library(openxlsx)
  library(rio)
  library(hablar)
  # Load mRNA half-lives
  message("REF: S. Marguerat et al., 2012, Cell")
  message("Quantitative Analysis of Fission Yeast Transcriptomes and Proteomes in Proliferating and Quiescent Cells")
  #https://doi.org/10.1016/j.cell.2012.09.019
  SM.url = "https://ars.els-cdn.com./content/image/1-s2.0-S0092867412011269-mmc1.xlsx"

  # S4 = Copies per cell data for all RNAs and proteins, including mRNA and protein features and annotation
  S4=import(SM.url,which = 5,skip=6) %>% janitor::clean_names() %>%
    # remove columns names containing footnotes from a to e
    dplyr::rename(chromosome="chromosomea",
                  sequencability="sequencabilityb",
                  ribosome_occupancy="ribosome_occupancyc",
                  mrna_halflife="m_rna_half_lifed",
                  annotation="annotatione") %>%
    rowwise() %>%
    mutate( cpc_pro=as.numeric(mm_protein_cpc),cpc_qui=as.numeric(mn_protein_cpc)) %>% ungroup
  # S8 = Global absolute abundance estimations for all proteins identified from proliferating cells
  S8=import(SM.url,which = 9) %>% janitor::clean_names() %>% dplyr::rename_with(xxS, sx='pro', s='_')
  # S9 = Global absolute abundance estimations for all proteins identified from quiescence cells
  S9=import(SM.url,which = 10) %>% janitor::clean_names() %>% dplyr::rename_with(xxS, sx='qui', s='_')

  colexp = c("cpc_pro","cpc_qui","protein_copies_cell_pro","protein_copies_cell_qui")
  rawdata = full_join(S4 , S8 ,by=c('systematic_name'='accession_pro')) %>%
    full_join(S9, by=c('systematic_name'='accession_qui')) %>% rowwise %>%
    mutate(nval = sum.na(c_across(colexp),notNA=T) ) %>% dplyr::filter(nval>0)

  if(raw){ return(rawdata) }

  expdata = rawdata %>%
    dplyr::select(pombase="systematic_name", gname="common_name",colexp,'nval') %>%
    rowwise %>%
    mutate(cpc_avg = mean_(c_across(cpc_pro:protein_copies_cell_qui)),
           cpc_min = min_(c_across(cpc_pro:protein_copies_cell_qui)),
           cpc_max = max_(c_across(cpc_pro:protein_copies_cell_qui)),
    ) %>% ungroup() %>%
    mutate( perc_max = 100 * cpc_max / sum_(cpc_max), ppm_max = 10000 * perc_max,
            perc_min = 100 * cpc_min / sum_(cpc_max), ppm_min = 10000 * perc_min,
            perc_avg = 100 * cpc_avg / sum_(cpc_avg), ppm_avg = 10000 * perc_avg,
            perc_pro = 100 * cpc_pro / sum_(cpc_pro), ppm_pro = 10000 * perc_pro,
            perc_qui = 100 * cpc_qui / sum_(cpc_qui), ppm_qui = 10000 * perc_qui)
return(expdata)
}


load.christiano2014.data = function(sp='cerevisiae'){
  message("REF: R Christiano et al, Cell Reports, 2014")
  message("Global Proteome Turnover Analyses of the Yeasts S. cerevisiae and S. pombe")
  S1.url="https://ars.els-cdn.com/content/image/1-s2.0-S2211124714009346-mmc2.xlsx" # cerevisiae
  S2.url="https://ars.els-cdn.com/content/image/1-s2.0-S2211124714009346-mmc3.xlsx" # pombe
  if(sp == 'cerevisiae'){
    prot.turnover = rio::import(S1.url,na='n.d.') %>% as_tibble() %>% janitor::clean_names()
    #openxlsx::read.xlsx(xlsxFile = S1.url, sheet = 1, detectDates = F, skipEmptyRows = T, skipEmptyCols = T, startRow = 1)
  }else if(sp == 'pombe'){
    prot.turnover = rio::import(S2.url,na='n.d.') %>% as_tibble() %>% janitor::clean_names()
  }
  return(prot.turnover)
}

load.geisberg2014.data = function(nodesc=T){
  library(openxlsx)
  # Load mRNA half-lives
  message("REF: J.V. Geisberg et al., 2014, Cell")
  message("Global Analysis of mRNA Isoform Half-Lives Reveals Stabilizing and Destabilizing Elements in Yeast")
  #https://doi.org/10.1016/j.cell.2013.12.026
  S1.url = "https://ars.els-cdn.com/content/image/1-s2.0-S009286741301595X-mmc1.xlsx"
  #download.file(S1.url, destfile = "mmc1.xlsx")
  mrna.half = openxlsx::read.xlsx(xlsxFile = S1.url,
                                  sheet = 1, detectDates = F,
                                  skipEmptyRows = T, skipEmptyCols = T, startRow = 1)

  colnames(mrna.half)=c('chr','ORF','GNAME','type','mrna_halflife_mn','DESC')
  if(nodesc){ mrna.half$DESC = NULL }
  return(mrna.half)
}

load.dana2014.data = function(){
  # Load ribosome translation efficiency
  message("REF: A. Dana and T. Tuller, 2014, G3 (Genes|Gemomes|Genetics)")
  message("Mean of the Typical Decoding Rates: A New Translation Efficiency Index Based on the Analysis of Ribosome Profiling Data")
  message("see also -> The effect of tRNA levels on decoding times of mRNA codons (Nucleic Acids Res)")
  # https://doi.org/10.1534/g3.114.015099
  # S. cerevisiae data based on ribosome profiling data:
  # - Ingolia et al. 2009 : Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling
  # - Brar et al. 2012    : High-resolution view of the yeast meiotic program revealed by ribosome profiling

  url.MTDR = "https://www.cs.tau.ac.il/~tamirtul/MTDR/MTDR_ORF_values"
  org = "S. cerevisiae"
  conditions = c("A14201","gb15",
                 "DNA replication","Recombination","anaphase","metaphase I","metaphase II","premeiotic entry",
                 "spore packing","spores")
  studies = c("Ingolia",paste0("Brar-",conditions) )
  url.datasets = gsub(" ","%20", x=sprintf("%s/%s %s.txt",url.MTDR,org,studies))
  names(url.datasets) = conditions
  exp.names = c('ingolia2009.TE_exponential','brar2012.TE_exponential_A14201','brar2012.TE_exponential_gb15',
                'brar2012.TE_DNA_replication','brar2012.TE_recombination','brar2012.TE_anaphase','brar2012.TE_metaphase_I','brar2012.TE_metaphase_II', 'brar2012.TE_premeiotic_entry',
                'brar2012.TE_spore_packing', 'brar2012.TE_spores')
  MTDR = map(url.datasets,read.delim,header=F,col.names=c("orf","decoding_rate")) %>%
    purrr::reduce(.x=., .f = left_join, by = "orf") %>%
    rename_with(~exp.names,.cols=starts_with('decoding_rate'))
  return(MTDR)
}

load.lee2014.data = function(rawdata=F){
  # Load chemogenomic fitness signatures
  message("REF: A.Y. Lee,  R.P. St.Onge et al., 2014, Science")
  message("Mapping the Cellular Response to Small Molecules Using Chemogenomic Fitness Signatures")
  # https://doi.org/10.1126/science.1250217
  # Related to this first paper: The chemical genomic portrait of yeast: Uncovering a phenotype for all genes
  # E. Hillenmeyer, et al. Science 320,362–365 (2008).doi:10.1126/science.1150021


  url.hiphop="http://chemogenomics.pharmacy.ubc.ca/hiphop/files/supplemental"
  url.s1 = paste0(url.hiphop,"/","leesupptableS1.xlsx")
  url.s4 = paste0(url.hiphop,"/","leesupptableS4.xlsx")
  url.fitness_hom =paste0(url.hiphop,"/","fitness_defect_matrix_hom.txt")
  url.fitness_het =paste0(url.hiphop,"/","fitness_defect_matrix_het.txt")

  library(openxlsx)
  # Compound library
  compound.library =  openxlsx::read.xlsx(  xlsxFile = url.s1, sheet = 2,
                                            colNames = T, skipEmptyRows = T,
                                            skipEmptyCols = T, na.strings = c(""))

  # Fitness Defect Score Matrix, homozygous strains (right-click to download, tab-delimited text, 280 Mb)
  # [Rows=Yeast Deletion Strain, Identified by Systematic ORF name; Columns=Fitness Screens, Identified by Screen ID (SGTC_N)]
  # Fitness defect on homozygous deletion strain for 3000 small molecules screen
  library(vroom)
  if(rawdata){
    chemofit = vroom(file=url.fitness_het, col_names=T, delim="\t",
                     col_select = list(orf = 1, everything()),
                     col_types = cols(.default = col_double()))
    return(chemofit)
  }else{
    major.responses=  openxlsx::read.xlsx(  xlsxFile = url.s4, sheet = 2, na.strings = c("","NA")) %>%
                      janitor::clean_names() %>%
                      dplyr::select(response_signature,gene,median_fd)

    minor.responses = openxlsx::read.xlsx(  xlsxFile = url.s4, sheet = 3, na.strings = c("",NA)) %>%
                      janitor::clean_names() %>% hablar::convert(chr(response_signature))
    responses = bind_rows(major.responses, minor.responses) %>%
                relocate(gene) %>%
                arrange(gene,response_signature) %>%
                dplyr::rename(FD_med=median_fd)
    return(responses)
  }
}


load.filleton2015.data = function(){
  message("REF: Filleton F. et al, 2015, Epigenetics & Chromatin")
  message("The complex pattern of epigenomic variation between natural yeast strains at single-nucleosome resolution")
  # Epidiv values for all genes (intra-species variation of chromatin modifications)
  S12="https://static-content.springer.com/esm/art%3A10.1186%2Fs13072-015-0019-3/MediaObjects/13072_2015_19_MOESM12_ESM.ods"
}


load.lancaster2014.data = function(){
  message("REF: Lancaster AK et al, 2014, Bioinformatics")
  message("PLAAC: a web and command-line application to identify proteins with Prion-Like Amino Acid Composition")
  #doi:10.1093/bioinformatics/btu310
  yeast_datafile="http://plaac.wi.mit.edu/Scer-all-proteins-2014-05-17.xls"
}


load.vanleeuwen2016.data = function(single_orf=F){
  # load gene classification of biological functions (based on costanzo 2010)
  message("REF: J. Van Leeuwen et al., 2016, Science")
  message("Exploring genetic suppression interactions on a global scale")
  # doi: https://doi.org/10.1126/science.aag0839

  # Sheet 1 is the functional classes
  S7="https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aag0839&file=aag0839tables7.xlsx"
  ## EDIT 21.12.21: NEW URL THAT REPLACED THE PREVIOUS ONE IS REDIRECTING TO DOWNLOAD THE FILE
  ## Using the target url for downloading the file.
  ##  xlsx file is: https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aag0839&file=aag0839tables7.xlsx
  #S7="https://www.science.org/doi/suppl/10.1126/science.aag0839/suppl_file/aag0839tables7.xlsx"
  ## EDIT 13.12.21: OLD URL BELOW NOT VALID ANYMORE
  ## https://science.sciencemag.org/highwire/filestream/686300/field_highwire_adjunct_files/6/aag0839TableS7.xlsx"
  #temp=tempfile()
  #download.file(S7, destfile = temp)

  .class_function = c(
    'A'="Amio acid biosynth & transport",
    'B'="Autophagy",
    'C'="Cell cycle progression/meiosis",
    'D'="Cell polarity/morphogenesis",
    'E'="Chrom. seg./kinetoch./spindle/microtub.",
    'F'="Chromatin/transcription",
    'G'="DNA replication & repair/HR/cohesion",
    'H'="Drug/ion transport",
    'I'="ER-Golgi traffic",
    'J'="Golgi/endosome/vacuole sorting",
    'K'="Lipid/sterol/fatty acid biosynth & transport",
    'L'="Metabolism/mitochondria",
    'M'="Nuclear-cytoplasmic transport",
    'N'="Peroxisome",
    'O'="Protein degradation/proteosome",
    'P'="Protein folding & glycosylation/cell wall",
    'Q'="Ribosome/translation",
    'R'="RNA processing",
    'S'="Signaling/stress response",
    'X'="Highly pleiotropic"
  )

  vanleeuwen = read.xlsx( xlsxFile = S7, sheet = 2) %>%
    separate_rows(Function,sep=",") %>%
    mutate( BIOPROCESS = str_trim(Function),
            FN=factor(BIOPROCESS,levels=unique(BIOPROCESS),labels = invert(.class_function))
    ) %>% dplyr::select(-Function) %>%
    group_by(ORF) %>% mutate(FN_all=paste0(unique(FN),collapse = ""),
                             BIOPROCESS_all = paste0(unique(BIOPROCESS),collapse=" # "))
  if(single_orf){
    # Return merged functions for each orf
    vanleeuwen =  vanleeuwen %>% dplyr::select(-c(FN,BIOPROCESS)) %>% distinct()
  }
  return(vanleeuwen)
}

load.lahtvee2017.data = function(){
  # Load protein and mRNA abundance and translation efficiency
  message("REF: P-J Lahtvee et al., 2017, Cell Systems")
  message("Absolute Quantification of Protein and mRNA Abundances Demonstrate Variability in Gene-Specific Translation Efficiency in Yeast")

  # transcript level
  S2.url = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471217300881-mmc2.xlsx"
  # 10 first columns = picog of Dry Weight, 11-19 log2 foldchange, 20-29 pvalue adjusted
  measurements = c("",rep("",10),rep("fc_",9),rep("pv_",9))
  header_s2 = rio::import(S2.url,skip=1,n_max=1) %>% janitor::clean_names("lower_camel") %>% colnames %>%
    paste0(paste0(measurements,"mrna_"),.)  %>% str_replace_all("Ref[0-9]+","ref")
  transcript_level = rio::import(S2.url,skip=4,col_names=c('geneId',header_s2[-1])) %>% as_tibble %>% type_convert %>%
    dplyr::select(-starts_with("pv_"))

  # protein abundance
  S3.url = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471217300881-mmc3.xlsx"
  # 10 first columns = picog of Dry Weight, 11-19 log2 foldchange, 20-29 pvalue adjusted
  header_s3 = rio::import(S3.url,skip=1,n_max=1) %>% janitor::clean_names("lower_camel") %>% colnames %>%
    paste0(paste0(measurements,"prot_"),.)  %>% str_replace_all("Ref[0-9]+","ref")
  protein_level = rio::import(S3.url,skip=4,col_names=c('geneId',header_s3[-1])) %>% as_tibble %>% type_convert %>%
    mutate(across(starts_with("prot"), log10,.names = "log10_{.col}")) %>%
    dplyr::select(-starts_with(c("pv_","prot")))


  log10_ = function(x){ x[!is.na(x) & x>0] = log10(x[!is.na(x) & x>0]); return(x) }
  # protein-mRNA ratio (translation rate)
  S4.url = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471217300881-mmc4.xlsx"
  header_s4 = rio::import(S4.url,skip=1,n_max=1) %>% janitor::clean_names("lower_camel") %>% colnames
  translation_rate = rio::import(S4.url,skip=4,col_names=header_s4) %>% as_tibble %>% type_convert %>%
    rename_with(.fn=Pxx, px='trans',s='_',.cols=-geneId) %>%
    mutate(across(starts_with("trans_"), log10_,.names = "log10_{.col}")) %>%
    dplyr::select(-starts_with(c("trans")))

  # turnover
  S5.url = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471217300881-mmc5.xlsx"
  protein_turnover = rio::import(S5.url,skip=1) %>% janitor::clean_names() %>% as_tibble %>%
    dplyr::select(-adjusted_p_value)

  protquant = left_join(transcript_level,protein_level,by=c('geneId')) %>%
    left_join(translation_rate,by=c('geneId')) %>%
    left_join(protein_turnover,by=c('geneId'='associated_gene_id'))
  return(protquant)
}


load.villen2017.data = function(){
  library(openxlsx)

  # Load protein turnover (half-lives)
  message("REF: M. Martin-Perez and J. Villén, 2017, Cell Systems")
  message("Determinants and Regulation of Protein Turnover in Yeast")
  #"https://doi.org/10.1016/j.cels.2017.08.008"
  S1.url = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471217303411-mmc2.xlsx"
  #download.file(S1.url, destfile = "1-s2.0-S2405471217303411-mmc2.xlsx")
  turnover = openxlsx::read.xlsx(xlsxFile = S1.url,
                                 sheet = 1, detectDates = F,
                                 skipEmptyRows = T, skipEmptyCols = T, startRow = 4)

  colnames(turnover) = c('ORF','UNIPROT','GNAME','PNAME',
                         'halflife.avg','halflife.sd','halflife.cv',
                         'R1_hl.avg','R1_method','R1_reduced.x2','R1_adj.r2','R1_cv',
                         'R2_hl.avg','R2_method','R2_reduced.x2','R2_adj.r2','R2_cv')
  return(turnover)
}

load.leuenberger2017.data = function(species='S. cerevisiae',rawdata=F){
  # Load protein stability data
  library(openxlsx)
  library(tidyverse)
  message("REF: P. Leuenberger et al., 2017, Science")
  message("Cell-wide analysis of protein thermal unfolding reveals determinants of thermostability")
  species=match.arg(species, choices = c('S. cerevisiae','E. coli', 'Human HeLa Cells','T. thermophilus'), several.ok = F)
  S3.url="https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aai7825&file=aai7825_leuenberger_table-s3.xlsx"
  ## EDIT 04.01.22: OLD URL BELOW NOT VALID ANYMORE - NEW URL USES REDIRECTION TO DOWNLOAD THE FILE
  # S3.url = "https://science.sciencemag.org/highwire/filestream/690833/field_highwire_adjunct_files/2/aai7825_Leuenberger_Table-S3.xlsx"
  # REDIRECTED URL : https://www.science.org/doi/suppl/10.1126/science.aai7825/suppl_file/aai7825_leuenberger_table-s3.xlsx
  # Target url for downloading the file is:
  # https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aai7825&file=aai7825_leuenberger_table-s3.xlsx

  #download.file(S3.url, destfile = "aai7825_Leuenberger_Table-S3.xlsx" )
  LIP_MS = openxlsx::read.xlsx(xlsxFile = S3.url,
                               sheet = species, detectDates = T,
                               skipEmptyRows = T, skipEmptyCols = T) %>%
    janitor::clean_names()

  peptides = LIP_MS %>%
    separate_rows(protein_id, sep=';') %>%
    add_count(protein_id,name='npep') %>%
    distinct() %>%
    mutate( nres = round(0.01*protein_coverage*length))

  if(!rawdata){
    peptides = peptides %>%
      dplyr::select(protein_id,tm_protein,
                    protinfo,protein_coverage,length,
                    measured_domains, theoretical_number_of_domain,
                    npep )
  }
  # cols = c("Peptide.ID","Aggregation","Position","Tm.Peptide", essential,
  #         "Sequence","Pepinfo", "Secondary.Structure","Domain.logic","Domain.Name",
  #          "is.disordered", "T.90..Unfolded")

  return(peptides)
}

load.mittal2017.data=function(){
  message("REF: N. Mittal et al., 2017, Nature Communications")
  message("The Gcn4 transcription factor reduces protein synthesis capacity and extends yeast lifespan")
  S1.url="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5587724/bin/41467_2017_539_MOESM2_ESM.xlsx"
  protein_synthesis = rio::import(S1.url,skip=1) %>% janitor::clean_names() %>% type_convert()
  return(protein_synthesis)
}

load.peter2018.data =function(d){
  # Load 1011 yeast strains data
  # Sheet 1 = Strains details
  # Sheet 3 = Variable ORFs
  library(tidyverse)
  library(hablar)
  library(janitor)
  library(rio)

  message("REF: J. Peter et al., 2018, Science")
  message("Genome evolution across 1,011 Saccharomyces cerevisiae isolates")
  SM.url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0030-5/MediaObjects/41586_2018_30_MOESM3_ESM.xls"
  #"http://1002genomes.u-strasbg.fr/isolates/page8/files/1002genomes.txt"

  peter=data.frame(stringsAsFactors = F, name   = c("strains collection", "variable ORFS"), sheetnum    = c(1,3) )
  if( missing(d) || !(d %in% seq_along(peter$name)) ){
    d = menu(sprintf("sheet %s (%s)",peter$sheetnum,peter$name), graphics = FALSE, title = "Which dataset do you want to use?")
  }
  choice = peter[d,]
  with(choice,cat(sprintf("Your choice was:\n [%s] S%s - '%s'  \n-----> obtained from Supp. Tables of Peter et al., Nature (2018)\n",d,sheetnum,name)))

  if(d==1){
    collection = rio::import(file = SM.url,which = 1,progress = T,skip = 2,col_names = T,n_max = 1011,guess_max = 900) %>%
      janitor::clean_names() # Clean the column names

    strains= tibble(type.convert(collection[1:1011,])) %>%   # 1011 first rows = strains details
      hablar::convert(hablar::chr(isolate_name),
                      hablar::chr(isolation),
                      hablar::chr(geographical_origins),
                      hablar::num(number_of_singletons),
                      hablar::chr(collection_provider)
      )
    return(strains)

  }else if(d==2){
      varorf = rio::import(file = SM.url,which = 3, progress = T,skip = 2,col_names = T,n_max = 3000,guess_max = 2600) %>%
        janitor::clean_names(parsing_option=3) %>% # Clean the column names
        dplyr::select(where(~ mean(is.na(.x)) < 0.9)) %>%  #  Remove columns with 90% NA
        mutate(orf.s288c = str_extract(annotation_name,pattern = SGD.nomenclature())) %>%
        relocate(orf.s288c,origin_assignment,occurrences_confirmed_by_mapping) %>%
        dplyr::filter(!is.na( orf.s288c) )
      return(varorf)
  }
    #PARSING REFERENCES
    # refs = collection = readxl::read_excel(path = todownload, sheet = 1, progress = T,skip = 1013,col_names = F) %>%
    #   str_split(pattern="\\. ",n=2,simplify = T) %>%
    #   as_tibble %>%
    #   rename(numref =V1, complete.ref=V2) %>%
    #   mutate(year  = str_extract(complete.ref,pattern = "\\(([0-9]{4})\\)"),
    #          title = get.longest(complete.ref),
    #          ref.2 = str_remove(complete.ref,year),
    #          pages = str_extract(ref.2,"(\\d+)-(\\d+)"),
    #          ref.3 = str_remove(ref.2,title),
    #          authors = str_extract(string = ref.3, pattern = '.+\\.\\.'),
    #          ref.4 = str_remove(ref.3,authors),
    #          journal.vol = str_extract(ref.3,pattern='\\.\\..+,'),
    #          journal = str_trim(str_extract(journal.vol,pattern='.+ ')),
    #          volume = str_remove(str_extract(journal.vol,pattern='[0-9]+,'),",")
    #          #ref.2=NULL, ref.3=NULL, ref.4=NULL,journal.vol=NULL
    #   )
}

load.ho2018.data = function(noauto = T, nogfp  = F, noms   = F, notap  = F) {
  library(openxlsx)
  message("REF: B. Ho et al., 2018, Cell Systems")
  message("Unification of Protein Abundance Datasets Yields a Quantitative Saccharomyces cerevisiae Proteome")
  # https://doi.org/10.1016/j.cels.2017.12.004
  ms  = sprintf("ms.%s",
                c("LU", "PENG", "KUL", "LAW", "LAHT", "DGD", "LEE2", "THAK", "NAG", "PIC", "WEB"))
  gfp = sprintf("gfp.%s",
                c("TKA", "BRE", "DEN", "MAZ", "CHO", "YOF", "NEW", "LEE", "DAV"))
  tap = sprintf("tap.%s", "GHA")
  gro = c("ypd_mid", "min_early", "min_mid", "min_steady", "min_chem")
  cols = c( 'orf', 'gname', 'qual', 'mean.mpc', 'median.mpc', 'CV.mpc',
            sprintf("%s.%s", ms, gro[c(1, 2, 3, 5, 5, 1, 1, 3, 1, 1, 1)]),
            sprintf("%s.%s", gfp, gro[c(3, 3, 4, 3, 3, 3, 1, 1, 1)]),
            sprintf("%s.%s", tap, gro[1])
  )
  # Table S3. Protein Abundance in Molecules per Cell, before GFP Autofluorescence Filtering
  S3 =  openxlsx::read.xlsx(
    xlsxFile = "https://ars.els-cdn.com/content/image/1-s2.0-S240547121730546X-mmc4.xlsx",
    sheet = 1, rows = 3:5751,
    colNames = T, skipEmptyRows = T, skipEmptyCols = T, na.strings = c("")
  )
  colnames(S3) = cols

  # Table S4. Protein Abundance in Molecules per Cell, after GFP Autofluorescence Filtering,
  S4 =  openxlsx::read.xlsx(
    xlsxFile = "https://ars.els-cdn.com/content/image/1-s2.0-S240547121730546X-mmc5.xlsx",
    sheet = 1, rows = 3:5861,
    colNames = T, skipEmptyRows = T, skipEmptyCols = T, na.strings = c("")
  )
  colnames(S4) = cols


  col.ms  = strfind(colnames(S3),ms)
  col.gfp = strfind(colnames(S3),gfp)
  col.tap = strfind(colnames(S3),tap)
  col.exp = c(col.ms,col.gfp,col.tap)

  unified = S3
  if (noauto) { unified = S4 }
  if (nogfp)  { unified = unified[, -col.gfp] }
  if (noms)   { unified = unified[, -col.ms]  }
  if (notap)  { unified = unified[, -col.tap] }

  # Compute geometric mean and sd and number of available values from each type of experiment
  unified$EXP = rowSums(!is.na(unified[,col.exp]))
  unified$na = rowSums(is.na(unified[,col.exp]))

  unified$GFP = rowSums(!is.na(unified[,col.gfp]))
  unified$na.GFP = rowSums(is.na(unified[,col.gfp]))
  unified$GFP.avg = apply(unified[, col.gfp], 1, function(x) { geomean(x) })
  unified$GFP.sd = apply(unified[, col.gfp], 1, function(x) { geosd(x) })

  unified$MS = rowSums(!is.na(unified[,col.ms]))
  unified$na.MS = rowSums(is.na(unified[,col.ms]))
  unified$MS.avg = apply(unified[, col.ms], 1, function(x) { geomean(x) })
  unified$MS.sd = apply(unified[, col.ms], 1, function(x) { geosd(x) })

  return(unified)
}

load.meldal.2019.data = function(species='yeast'){
  # Load macromolecular protein complexes from yeast
  message("REF: B. Meldal et al., 2019, Nucleic Acids Research")
  message("Complex Portal 2018: extended content and enhanced visualization tools for macromolecular complexes")
  # https://doi.org/10.1093/nar/gky1001
  complextab.url = "ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab"

  taxon=match.arg(species, choices = c('559292'='yeast','9606'='human'), several.ok = F)
  UNI_AC = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"

  COMPLEXES = read.delim(file = sprintf("%s/%s.tsv",complextab.url,names(taxon)),
                         stringsAsFactors = F,header = T, sep = '\t', na.strings = '-',
                         col.names = c('CPLX_ID','CPLX_NAME','CPLX_ALIAS','TAX_ID',
                                       'MEMBERS_STOCHIO','CONF','EXP_EVID','GO_ANNOT',
                                       'CROSSREF','DESC',
                                       'CPLX_PROP','CPLX_ASSEMBLY',
                                       'LIG','DISEASE','AGONIST','ANTAGONIST',
                                       'COMMENT','SOURCE','EXPANDED_MEMBERS')
  ) %>% as_tibble()

  # Some members of complexes are also complexes
  # Expanded members contains only the uniprot references
  # i.e. all members are proteins i.e. complexes are flatter to their constituents)
  CX=COMPLEXES %>%
    dplyr::select("CPLX_ID","CPLX_NAME", "EXPANDED_MEMBERS",
                  "EXP_EVID","CPLX_ASSEMBLY","GO_ANNOT",
                  "LIG","DISEASE","AGONIST","ANTAGONIST") %>%
    separate_rows( EXPANDED_MEMBERS, sep='\\|') %>%
    mutate( members = str_remove(EXPANDED_MEMBERS,"\\([0-9]+\\)"),
            stochio = as.integer( str_replace(EXPANDED_MEMBERS,".+\\(([0-9]+)\\)","\\1"))) %>%
    add_count(name='n_members',CPLX_ID) %>%
    mutate(is_uniprot = str_detect(members,UNI_AC),
           is_RNA = str_detect(members,"^URS[0-9A-Z]"),
           is_complex = str_detect(members,"^CPX\\-[0-9]+"),
           is_small_molecule = str_detect(members,"^CHEBI:")) %>%
    group_by(CPLX_ID) %>% mutate( n_uniprot= sum(is_uniprot) ) %>%
    dplyr::select(-EXPANDED_MEMBERS)

  return(CX)
}

load.dubreuil2019.data = function(d){
  # Load stickiness and disorder data
  message("REF: Dubreuil, Matalon and Levy, 2019, Journal of Molecular Biology")
  message("Protein Abundance Biases the Amino Acid Composition of Disordered Regions to Minimize Non-functional Interactions")
  dubreuil=data.frame(stringsAsFactors = F,
                      name   = c("1.1 yeast.prot", "1.2 yeast.res", "1.3 yeast.res (full proteome)","1.4 yeast.prot (full proteome)",
                                 "2.1 human.prot", "2.2 human.res", "2.3 human.res (full proteome)","2.4 human.prot (full proteome)") ,

                      num    = c('16924922', '16920611', '25863816', '26844182',
                                 '16924925', '16920581', '25863810', '26844185'),
                      base_url =  rep('https://ndownloader.figshare.com/files/',8),
                      format = c('XLSX','TSV','TSV.GZ','TSV',
                                 'XLSX','TSV','TSV.GZ','TSV')
  )
  if( missing(d) || !(d %in% seq_along(dubreuil$name)) ){
    d = menu(sprintf("%s (%s)",dubreuil$name,dubreuil$format), graphics = FALSE, title = "Which dataset do you want to use?")
  }
  choice = dubreuil[d,]
  with(choice,cat(sprintf("Your choice was:\n [%s] %s\n-----> from %s%s (formatted as %s)\n",d,name,base_url,num,format)))
  data.url = sprintf('%s/%s',dubreuil$base_url[d], dubreuil$num[d])

  if( dubreuil$format[d] == 'XLSX' ){
    library(openxlsx)
    res = openxlsx::read.xlsx(data.url,colNames=T,sheet=1, startRow=1)
  } else if( dubreuil$format[d] == 'TSV'){
    res = readr::read_delim(file = data.url, col_names = T, delim = '\t',  progress = T,guess_max=50000)
  } else if( dubreuil$format[d] == 'TSV.GZ' ){
    res = readr::read_delim(file = I(read.url(data.url)), col_names = T, delim = '\t',  progress = T,guess_max=50000)
  }
  return(res)
}

load.hausser2019.data=function(show_desc=F){
  message("REF: J. Hausser et al., 2019, Nature Communications")
  message("Central dogma rates and the trade-off between precision and economy in gene expression")
  yeast_data = "https://data.mendeley.com/public-files/datasets/2vbrg3w4p3/files/9a834480-1da3-4e36-9f89-9cdd79d03382/file_downloaded"
  var_desc = rio::import(yeast_data, sheet=1) %>% as_tibble()
  if(show_desc){ print(var_desc) }
  yeast_rates = rio::import(yeast_data, sheet=2)
  colnames(yeast_rates) =c('gene','lm','lp','wRPF','wmRNA','bmEser','amEser','bm','m','bp','cv','apExp','am','isEssential','hasTATA','YEPDFit','pEst')
  return(yeast_rates %>% as_tibble %>% dplyr::select(-c("am"))) # am is constant
}

load.jarzab2020.data = function(org='S.cerevisiae'){
  library(hutils)
  message("REF: A. Jarzab et al., 2020, Nature Methods")
  message("Meltome atlas—thermal proteome stability across the tree of life")
  #https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-020-0801-4/MediaObjects/41592_2020_801_MOESM7_ESM.xlsx
  F2_url = "https://figshare.com/ndownloader/files/21653313"
  species = c("T.thermophilus", "P.torridus", "G.stearothermophilus",
              "E.coli", "B.subtilis",
              "S.cerevisiae", "M.musculus",  "H.sapiens", "C.elegans", "D.melanogaster", "D.rerio",
              "A.thaliana", "O.antarctica")
  organisms = match.arg(org,species,several.ok = T)

  F2a = rio::import(F2_url,col_names = TRUE,skip=1,sheet=1) %>%
        as_tibble %>%
        dplyr::filter(Species %in% organisms) %>%
        janitor::clean_names()

  F2b = rio::import(F2_url,col_names = TRUE,skip=1,sheet=2) %>%
    as_tibble %>%
    mutate(sp = str_replace_all(Dataset," +","") ) %>%
    dplyr::select(-Dataset) %>%
    dplyr::filter(sp %pin% organisms) %>% # Allow partial matching for H.sapiens
    janitor::clean_names() %>%
    mutate(org = str_extract(sp,organisms))

  F2=left_join(F2a,F2b, by=c('species'='org','protein_id')) %>%
     separate(col = protein_id, sep = '_',into = c('UNIPROT','GENENAME')) %>%
     dplyr::select(species,UNIPROT,GENENAME,
                   Tm_celsius=melting_point_c,
                   Tm_type=protein_classification,
                   AUC=area_under_the_melting_curve)
  return(F2)
}

load.szavitsnossan2020.data = function(){

  message("REF: J. Szavits-Nossan et al., 2020, Nucleic Acids Research")
  message("Inferring efficiency of translation initiation and elongation from ribosome profiling")
  zip.url = "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/48/17/10.1093_nar_gkaa678/1/gkaa678_supplemental_files.zip?Expires=1648233481&Signature=SPpuRD-C9Hl8cQYD5qtQnTXvPpPHpoqQCorupykk62Ix5jhnIe5JMiaVO3-rpTToxETOAe56z8BhVKVd5n6FBFvl0ovfvCoTeAsEquMZdzoBzUhSErDkllLunIXilFhF17L1Upv7wZFjPmTc3J1rK~Ckmm8EFqnnoJu5LvQBHcMN~HT8cpj3koYuBEnh1GFJ8pb17bILkNcB6qCwHnvtoW8g9e7S0fHZ2hsjYM9AlRCdwheSadSHQZgTyIJI44xHI2X51t5VC4v9LP~gGDlb8EcuYkDiRLcpQSbQkoI5Yt~BDOvymQY4AOna6UFVb1qN4Ll0j4J3Yv55ZjeXmXSm4g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"
  temp<-tempfile()
  download.file(zip.url,temp)
  data_files = unzip(temp,list = T)$Name[1:3]

  ## ribosome profiling experiments

  #Weinberg D.E. et al., 2016, Cell Reports
  #Improved ribosome-footprint and mRNA measurements provide insights into dynamics and regulation of yeast translation
  weinberg <- readr::read_delim(file = unz(temp, data_files[1]), skip=3) %>% janitor::clean_names()
  # Pop C. et al., 2014, Mol. Syst. Biol.
  # Causal signals between codon bias, mRNA structure, and the efficiency of translation and elongation
  pop <- readr::read_delim(file = unz(temp, data_files[2]), skip=3) %>% janitor::clean_names()
  # Guydosh N. et al., 2014, Cell
  # Dom34 rescues ribosomes in 3′ untranslated regions
  guydosh <- readr::read_delim(file = unz(temp, data_files[3]), skip=3) %>% janitor::clean_names()

  # TIE = Translation Initiation Efficiency
  # TEE = Translation Elongation Efficiency

  # Symbol	Meaning
  # L 	length of the mRNA (in codons, including START)
  # ℓ 	length of the ribosome (in codons)
  # α 	initiation rate [s−1]
  # ki 	elongation rate [s−1] of codon i
  # kL (or β) 	termination rate [s−1]
  # {ki} 	speed profile (elongation) of a given transcript
  # κi = ki/α 	relative (to initiation) elongation rate at codon i
  # {κi} 	relative (to initiation) elongation profile
  # ri 	experimental (normalized) density of codon i
  # {ri} 	experimental (normalized) density profile
  # r=∑Li=2ri/(L−1) 	mean density of a given gene
  # ρi 	theoretical (normalized) density of codon i
  # ρILAi 	theoretical (normalized) density of codon i in the initiation-limited approximation
  # {ρi} 	theoretical (normalized) density profile
  # ρsimi 	simulated (normalized) density of codon i
  # {ρsimi} 	simulated (normalized) density profile

  translation = weinberg %>%
                dplyr::select(orf=x1, tie, mean_tee, mean_kappa,
                              s,s_mf,s_opt,
                              mean_experimental_density, mean_th_density)
  return(translation)
}

load.vanleeuwen2020.data = function(){
  library(openxlsx)
  # Load gene dispensability inferred from bypass suppression of essential genes
  message("REF: Van Leeuwen et al., 2020, Molecular Systems Biology")
  message("Systematic analysis of bypass suppression of essential genes")
  #"https://www.embopress.org/doi/full/10.15252/msb.20209828"
  EV13.url = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209828&file=msb209828-sup-0014-DatasetEV13.xlsx"
  #download.file(EV13.url, destfile = "msb209828-sup-0014-datasetev13.xlsx" )
  dispensable = openxlsx::read.xlsx(xlsxFile = EV13.url,
                                    sheet = 2, detectDates = F,
                                    skipEmptyRows = T, skipEmptyCols = T
  )
  colnames(dispensable) = c('ORF','gene','KO_exp','disp_score','literature','disp')
  return(dispensable)
}

load.dubreuil2021.data = function(d){
  # Load evolutionary rate data for yeast
  message("REF: Dubreuil and Levy, 2021, Frontiers in Molecular Bioscences")
  message("Abundance imparts evolutionary constraints of similar magnitude on the buried, surface, and disordered regions of proteins")
  dubreuil=data.frame(stringsAsFactors = F,
                      name   = c("1.1 yeast.prot", "1.2 yeast.res", "1.3 yeast.res") ,
                      num    = c('27528953', '26404475', '26404478'),
                      base_url =  rep('https://ndownloader.figshare.com/files/',3),
                      format = c('TSV','TSV','TSV.GZ')
  )
  if( missing(d) || !(d %in% seq_along(dubreuil$name)) ){
    d = menu(sprintf("%s (%s)",dubreuil$name,dubreuil$format), graphics = FALSE, title = "Which dataset do you want to use?")
  }
  choice = dubreuil[d,]
  with(choice,cat(sprintf("Your choice was:\n [%s] %s\n-----> from %s%s (formatted as %s)\n",d,name,base_url,num,format)))
  data.url = sprintf('%s/%s',dubreuil$base_url[d], dubreuil$num[d])
  if( dubreuil$format[d] == 'TSV'){
    res = readr::read_delim(file = data.url, col_names=T, progress=T, delim='\t',guess_max=50000)
  } else if( dubreuil$format[d] == 'TSV.GZ' ){
    res = readr::read_delim(file=read.url(data.url), col_names=T, progress=T, delim='\t',guess_max=50000)
  }
  return(res)
}

# Resource databases with proteome data available ------------------------------

check.alphafold = function(uniprot,extension=c('cif','pdb'),quiet=T,as.path=T){
  URL_ALPHAFOLD = "https://alphafold.ebi.ac.uk/files"
  if(!file.exists(uniprot)){
    # Check whether a uniprot is associated to an alphafold model
    ext=match.arg(extension,extension)
    AF_uniprot=sprintf("%s/AF-%s-F1-model_v1.%s",URL_ALPHAFOLD,uniprot,ext)
    found =  RCurl::url.exists(AF_uniprot)
  }else{ # for locally predicted alphafold files, use full path instead of uniprot accession number
    AF_uniprot=uniprot # provide full path to the predicted file
    found = file.exists(AF_uniprot)
  }
  if(!quiet){ cat(uniprot," --> ",found,"\n") }
  if(found & as.path){ return(setNames(found,AF_uniprot)) }
  return(setNames(found,uniprot))
}

load.alphafold = function(uniprot,extension=c('cif','pdb'),quiet=T){
  # Load alphafold predicted structure from uniprot along with the confidence score per residue
  library(bio3d)
  has_alphafold = check.alphafold(uniprot,extension,as.path=T)
  if( has_alphafold ){
    AF = switch(extension, cif=bio3d::read.cif(names(has_alphafold),verbose=!quiet), pdb=bio3d::read.pdb(names(has_alphafold),verbose=!quiet))
    Ca = AF$atom[AF$calpha,]
    df.res = Ca %>% mutate(uni=uniprot) %>% select(uni,chain,resn=resno,resi=resid,plddt=b)
    return(df.res)
  }else{
    df.res=tibble(uni=uniprot, chain=NA, resn=NA, resi=NA,plddt=NA)
    warning(sprintf("%s was not found on the EBI/AlphaFold server!\n",uniprot))
    return(df.res)
  }
}

get.alphafold.proteome = function(id_uniprot){
  #library(doParallel)
  #library(parallel)
  library(tidyverse)
  library(tictoc)
  library(progressr)
  library(future.apply)
  plan(multisession)

  # Use a cluster to accelerate the process
  #ncpus = detectCores(logical = F)-1
  #cl <- makeCluster(ncpus, type='SOCK')
  #registerDoParallel(cl)
  #n <- 100

  n_uni = length(id_uniprot)
  # Check the existence of the alphafold model for each uniprot
  with_progress({
    p <- progressor(along = id_uniprot)
    tic("Retrieving alphafold model... (using future.apply")
    # Retrieve the alphafold predicted structure model from the input identifiers
    af =  future_lapply(1:n_uni, function(i){
      p()
      load.alphafold(uniprot=id_uniprot[i],extension='pdb')
    })
    toc()
  })

  df.af = do.call(rbind,af)

  n_af  = df.af %>% summarize(n=n_distinct(uni), na=sum(is.na(c_across(-uni))))
  cat(sprintf("--> found %s/%s predicted alphafold\n",n_af,n_uni))

  if( n_af>0 ){
    af_brk = c(0,50,70,90,100)
    af_lab = c('D','L','M','H')
    af_uni = id_uniprot[has_alphafold]
    afprot= df.af %>%
           as_tibble() %>%
           mutate(aa_af = aa321(resi),
                  plddt_bin = cut(plddt,breaks=af_brk, labels=af_lab, include.lowest=T))
    #Stop the cluster
    #stopCluster(cl)

    return(afprot)
  }else{
    stop("None of the identifiers given correspond to an alphafold model!")
  }
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

load.pfam = function(tax='559292'){
  # Load domains assignments based of HMM profiles (PFAM)
  # 559292 = S.cerevisiae
  URL_PFAM = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
  library(readr)
  .release = read.url(paste0(URL_PFAM,"/Pfam.version.gz")) %>% paste0(sep="\n")
  message("CURRENT RELEASE")
  message("---------------")
  message(.release)
  message("---------------")
  url.pfam = sprintf("%s/proteomes/%s.tsv.gz",URL_PFAM,tax)
  url.clan = sprintf("%s/Pfam-A.clans.tsv.gz",URL_PFAM)
  clans = readr::read_delim(url.clan,delim="\t",col_names=c("pfam_id",'clan_id','clan_name',"pfam_name","pfam_desc")) %>%
          mutate( clan_id = tidyr::replace_na(clan_id,replace = 'No_clan') )

  .vers = read.url(url.pfam)[1]
  message(.vers)
  .nprot = read.url(url.pfam)[2]
  message(.nprot)
  header = read.url(url.pfam)[3]
  # parse the header (3rd commented rows)
  columns = header %>%
    str_split(pattern = "> <") %>% unlist %>%
    str_replace_all("[ \\-]","_") %>% str_remove_all("[#<>]")
  pfam = readr::read_delim(file = url.pfam, skip=2,comment = "#",delim="\t", col_names = columns,escape_double = F,guess_max = 100) %>%
         left_join(clans,by = c('clan'="clan_id",'hmm_acc'='pfam_id','hmm_name'='pfam_name') )

  return(pfam)
}

load.string = function(tax="4932",phy=T,ful=T,min.score=700,vers="v11.0"){
  # Load protein links (physical PPI)
  # 4932 = S.cerevisiae
  URL_STRING = "https://stringdb-static.org/download"
  .version = vers
  .protein = ifelse(phy,"protein.physical","protein")
  .links = ifelse(ful,"links.full","links")
  STRING_dataset = sprintf("%s.%s.%s",.protein,.links,.version)
  STRING_url = sprintf("%s/%s/%s.%s.%s",URL_STRING,STRING_dataset,tax,STRING_dataset,"txt.gz")

  message("STRING VERSION: ", .version)
  message("STRING DATASET: ",STRING_dataset)
  message(sprintf("FILTERED BY:\n physical links = %s \nAND\n full evidence =  %s \nAND\n min. score >= %s", phy, ful, min.score))
  # each numeric column is an evidence score out of 1000 to asssess the existence of the link between protein
  # _transferred means information was inferred from a different organism (homology or orthologous group)

  # EDIT 04/01/22 :
  # Error in open.connection(con, "rb") :
  #    server certificate verification failed. CAfile: /etc/ssl/certs/ca-certificates.crt CRLfile: none
  # fixed by reinstalling ca-certificates in terminal:
  #  sudo apt update ; sudo apt-get install apt-transport-https ca-certificates -y ; sudo update-ca-certificates
  STRING_net = readr::read_delim(STRING_url,delim = " ") %>% filter(combined_score >= min.score)
  return(STRING_net)
}

load.pombe.orthologs = function() {
  # Load orthologs pombe-cerevisiae
  # Fused genes have parentheses to indicate which fused side they are
  bracket.before = "(?<=\\()"
  bracket.after = "(?=\\))"
  regex_fusion = paste0( bracket.before, "(FUSION-)?(N|C)?", bracket.after)
  URL_POMBASE = "https://www.pombase.org/data/" #"ftp://ftp.pombase.org/pombe/"

  url.orthologs = paste0(URL_POMBASE,"orthologs/cerevisiae-orthologs.txt")
  sp.sc = readr::read_delim(url.orthologs, delim="\t", col_names=c('PombaseID','ORFS'), comment="#", trim_ws=T) %>%
    mutate(grp = sprintf("OG_%04d",row_number())) %>%
    separate_rows(ORFS,sep = "\\|") %>%
    mutate(orf= dplyr::na_if(str_extract(ORFS,"^[^\\(]+"),"NONE"),  # get the orf,
           fused_side=str_sub(start=-1,str_extract(ORFS,regex_fusion)) )%>%  # get the parentheses content
    dplyr::filter(!is.na(orf)) %>%
    group_by(PombaseID) %>% mutate(sp1 = n() ==1 ) %>%
    group_by(orf) %>% mutate(sc1= n() == 1) %>%
    rowwise %>% mutate( is_1to1 = sp1 & sc1,
                        is_fusion = !is.na(fused_side)) %>%
    dplyr::select(-ORFS) %>% ungroup() %>% arrange(orf)

  return(sp.sc)
}

##### MobiDB #####

get.mobidb.id = function(id="A0A075B734"){

  URL_API = "https://mobidb.bio.unipd.it/api/download?"
  URL_PARAM = sprintf("acc=%s&format=tsv",id)
  URL_QUERY = paste0(URL_API,URL_PARAM)
  df_mobi = readr::read_delim(URL_QUERY) %>%
            separate(col=feature, sep='-',into=c('evidence','feature','source')) %>%
            separate_rows('start..end',sep=',') %>%
            separate('start..end',into = c('S','E'),convert = T) %>%
            mutate( feature_len = E-S+1 )
  if( nrow(df_mobi) > 0L ){
    return(df_mobi)
  }else{
    return(NULL)
  }
}

fetch.mobidb = function(ids,to.df=T){
  # Check mobidb data description at the following url:
  # https://mobidb.bio.unipd.it/about/mobidb
  library(pbmcapply)
  library(parallel)
  library(tictoc)
  cpus=detectCores()-2
  tic('get disorder from MobiDB')
  message(sprintf('retrieve disorder predictions from MobiDB for %s identifiers (uniprot) ...',n_distinct(ids)))
  mobidb = pbmcapply::pbmclapply(ids, get.mobidb.id,mc.cores=cpus)
  names(mobidb) = ids
  toc()

  if(to.df){ return( purrr::compact(mobidb) %>% bind_rows() ) }

  return(mobidb)
}

load.mobidb = function(taxon){
  # Check mobidb data description at the following url:
  # https://mobidb.bio.unipd.it/about/mobidb
  tic('get mobidb for taxon from url....')
  mobidb = readr::read_delim(sprintf("https://mobidb.bio.unipd.it/api/download?ncbi_taxon_id=%s&format=tsv",taxon)) %>%
    separate(col=feature, sep='-',into=c('evidence','feature','source')) %>%
    separate_rows('start..end',sep=',') %>%
    separate('start..end',into = c('S','E'),convert = T) %>%
    mutate( feature_len = E-S+1, ncbi_taxid = taxon ) %>%
    relocate(ncbi_taxid)
  toc()
  return(mobidb)
}

merge_mobidb = function(mobidb, gap_min=1L){
  library(GenomicRanges)
  if(missing(mobidb)){ stop("requires mobidb data...") }

  merged = mobidb %>%
    dplyr::select(seqnames=acc,start=S,end=E) %>%
    as("GRanges") %>%
    reduce(min.gapwidth=gap_min) %>% # maximum gap for merging two intervals
    as_tibble() %>%
    dplyr::rename(acc=seqnames,S=start,E=end,feature_len=width) %>%
    mutate(acc=as.character(acc)) %>%
    dplyr::select(-strand) %>%
    arrange(acc)
  return(merged)
}

##### paxDB #####
find_paxdb_downloads = function(){
  URL_PAXDB = "https://pax-db.org/downloads/latest/"
  paxdb_files = rvest::read_html(URL_PAXDB) %>%
  rvest::html_nodes("a") %>%
  rvest::html_text(trim = T) %>%
  str_subset(pattern = "/$",negate = T)
  return(paxdb_files)
}

get_paxdb_version = function(verbose=T){
  paxdb_files=find_paxdb_downloads()
  .version = str_extract(pattern="[0-9]\\.[0-9]",string = paxdb_files) %>% gtools::mixedsort(na.last = F) %>% last
  if(verbose){
    cat(sprintf("PAXDB VERSION: %s\n",.version))
  }
  return(.version)
}

load.paxdb.orthologs = function(node,show.nodes=F) {
  library(rio)
  library(RCurl)
  library(hablar)
  # get table of orthologs proteins from paxdb
  URL_PAXDB = "https://pax-db.org/downloads/latest/"
  .version= get_paxdb_version()

  paxdb_ortho=grep("paxdb-orthologs",find_paxdb_downloads(),v=T) %>% str_subset(pattern=.version) %>% paste0(URL_PAXDB,.)
  # download the archive
  temp<-tempfile()
  download.file(paxdb_ortho,temp)
  # Find the list of files in the archive (taxonomic nodes)
  paxdb_orthologs_dir = unzip(temp,list = T)$Name[1]

  paxdb_nodes = unzip(temp,list = T)$Name %>%
           str_subset("\\.txt") %>%
           basename %>%
           str_remove(".orthgroups.consistent.txt") %>%
           gtools::mixedsort()
  if(show.nodes){ return(paxdb_nodes) }
  # Check the taxonomic node selected
  node_exists= !purrr::is_empty(node)
  is_paxdb_node=F

  if(node_exists){
    numnode=grep(node,paxdb_nodes)
    if(length(numnode)==1 ){
      is_paxdb_node=T
      node_full=grep(node,paxdb_nodes,v=T)
      message("taxonomic node (",node_full,") is valid!")
    }else if(length(numnode)>1){
      warning("(",node,") is not unique! Please select a valid unique node",immediate.=T)
      is_paxdb_node=F
    }
  }

  if(!is_paxdb_node){
    warning("(",node,") is not a valid taxonomic node from paxDB orthologs!",immediate.=T)
    numnode = menu(paxdb_nodes, title = "Pick a taxonomic node from the list below",graphics = T)
  }

  node_file = paxdb_nodes[numnode]
  paxdb_node = basename(paxdb_orthologs_dir) %>% file.path(.,node_file)

  # Read orthogroups at a specific taxonomic level
  message("Retrieving orthologs for node [",node,"]...")
  node_ortho <- readr::read_delim(file = unz(temp, paxdb_node),col_names = c("NOG", "orthologs"), delim="\t",progress=T) %>%
    separate_rows(orthologs,sep=" ") %>%
    extract(orthologs,into=c('taxid','protid'),regex='(^[0-9]+)\\.(.+)') %>%
    group_by(NOG,taxid) %>% mutate(is_1to1=n()==1)
  unlink(temp)
  return(node_ortho)
}

find_paxdb_datasets = function(taxon=4932){
  URL_PAXDB = "https://pax-db.org/downloads/latest/"
  paxdb_dataset =paste0(URL_PAXDB,"datasets/")
  taxon_dir=file.path(paxdb_dataset,taxon,"/")

  .version= get_paxdb_version()

  taxon_data <- rvest::read_html(taxon_dir) %>%
    rvest::html_nodes("a") %>%
    rvest::html_attr(name='href') %>%
    stringr::str_subset(pattern = "/$",negate = T) %>%
    stringr::str_subset("\\.txt")

  Ndata = length(taxon_data)
  message(Ndata," paxDB datasets for taxon [",taxon,"]")

  taxon_url = paste0(taxon_dir,taxon_data)
  get.paxdb_header = function(urldata){
    read.url(urldata) %>%
      stringr::str_subset(pattern="^#") %>%
      as_tibble %>%
      tidyr::extract(col=value,into=c('info',"value"), regex="^#([^\\:]+)\\:(.+)$") %>%
      dplyr::mutate(info=str_trim(info),value=str_trim(value)) %>%
      dplyr::filter( !is.na(info) ) %>%
      pivot_wider(names_from='info',values_from=c(value)) %>%
      mutate(taxid=as.character(taxon))
  }

  message(sprintf('retrieving paxdb datasets information [%s]...\n',taxon))
  if(require(pbmcapply)){
    ncpu = parallelly::availableCores()-2
    cat(sprintf("Using 'pbmcapply' with %s parallel threads",ncpu))
    infos = pbmcapply::pbmcmapply(mc.cores = ncpu, taxon_url , FUN = get.paxdb_header)  %>%
            bind_rows(.id='taxon_url')
  }else{
    cat("Using only 1 cpu! please consider installing 'pbmcapply' to use parallel threads.")
    infos = map_dfr(taxon_url, get.paxdb_header)
  }
  infodata = infos %>%
              mutate( w=parse_number(weight)*0.01,ndata = n_distinct(id,filename),
                      is_integrated = integrated=='true' | ndata==1) %>%
              dplyr::select(taxid,organ,ndata,id,filename,is_integrated,
                  score,w,cov=coverage,yr=publication_year)
  return(infodata)
}

load.paxdb = function(taxon=4932,rm.zero=T){
  URL_PAXDB = "https://pax-db.org/downloads/latest/"
  .version= get_paxdb_version(F)
  mapping_uniprot=paste0(URL_PAXDB,sprintf("paxdb-uniprot-links-v%s/paxdb-uniprot-links-v%s.tsv",.version,.version))
  map2uniprot = readr::read_delim(mapping_uniprot,delim="\t",col_names = c('id_string','id_uniprot'))

  infodata <- find_paxdb_datasets(taxon)
  taxon_url = file.path(URL_PAXDB,'datasets',taxon,infodata$filename)

  # if( Ndata == 1){
  # ppm = rio::import(taxon_url) %>%
  #   rename_with(~c("paxid",'string','ppm','count')) %>%
  #    mutate(dataset = basename(taxon_url)) %>%
  #     extract(string,into=c('taxid','protid'),regex='(^[0-9]+)\\.(.+)') %>%
  #      mutate()
  # }else{
  # If only one file, read twice the url and keep the distinct rows
  message(sprintf('retrieving paxdb abundance values  [%s]...\n',taxon))
  ppm = rio::import_list(file=c(taxon_url[1],taxon_url),
                           setclass="tibble", showProgress=T,
                           rbind = TRUE,rbind_label = "dataurl",rbind_fill = T) %>%
      dplyr::rename(paxid="#internal_id", string="string_external_id", ppm="abundance") %>%
      dplyr::mutate(dataset = basename(path = dataurl)) %>% dplyr::select(-dataurl) %>%
      tidyr::extract(string,into=c('taxid','protid'),regex='(^[0-9]+)\\.(.+)',remove = F) %>%
      distinct()
  #}#

  taxon_ppm = left_join(ppm,infodata, by=c('dataset'='filename','taxid')) %>%
              left_join(map2uniprot, by=c('string'='id_string')) %>%
              arrange(protid,id_uniprot,ppm)

  if(rm.zero){
    taxon_ppm = taxon_ppm %>% mutate(ppm = ppm +min_above(ppm,above=0,na.rm=T))
  }
  return(taxon_ppm)
}

get.paxdb = function(tax=4932, abundance='integrated',rm.zero=T){
  #closeAllConnections() # Make sure to close connections
  paxdb = load.paxdb(tax,rm.zero=rm.zero)
  # if nothing selected return the integrated values
  # (if there is a single dataset, it is considered as integrated)
  targets=match.arg(abundance,c('integrated','median','mean','weighted'), T)
  message(sprintf("---> Returning (%s) abundance values...",toString(targets)))
  RES = list()

  if( "integrated" %in% targets ){
    RES$INT = paxdb %>%
      dplyr::group_by(taxid,protid,id_uniprot) %>%
      mutate(n_data = n_distinct(id)) %>%
      dplyr::group_by(taxid,organ,protid) %>%
      dplyr::filter(is_integrated) %>%
      group_by(protid) %>% mutate(n_int = n_distinct(organ)) %>%
      ungroup %>%
      dplyr::select(taxid,organ,protid,id_uniprot, ppm_int = ppm, n_data, n_int)
  }

  if( "median" %in% targets ){
    RES$MED = paxdb %>%
      dplyr::filter(!is_integrated) %>%
      group_by(taxid,organ,protid,id_uniprot) %>%
      summarise(
        # range
        ppm_max = hablar::max_(ppm),
        ppm_min = hablar::min_(ppm),
        ppm_med = hablar::median_(ppm),
        ppm_mad = mad(ppm,na.rm=T)
      )
  }

  if( "mean" %in% targets ){
    RES$AVG = paxdb %>%
      dplyr::filter(!is_integrated) %>%
      group_by(taxid,organ,protid,id_uniprot) %>%
      summarise(
        # regular average
        ppm_avg = mean_(ppm),
        ppm_sd = sd_(ppm),
        ppm_cv  = ppm_sd/ppm_avg,
      )
  }

  if( "weigthed" %in% targets ){
    RES$WT = paxdb %>%
      dplyr::filter(!is_integrated) %>%
      group_by(taxid,organ,protid,id_uniprot) %>%
      summarise(
        # weighted average/median
        wppm_sum = sum_(ppm*w),
        wppm_Wtot=sum(w),
        wppm_avg = wppm_sum/wppm_Wtot,
        wppm_sd  = sd_(ppm*w),
        wppm_med = median_(ppm*w),
        wppm_mad = mad(ppm*w,na.rm=T),
        wppm_cv  = wppm_sd/wppm_avg,
      )
  }
  return(purrr::reduce(.x=RES,.f=left_join, by = c("taxid","organ","protid","id_uniprot")) )
}

summarise_paxdb_abundance = function(taxon=4932){
  paxdb_data = load.paxdb(taxon)
  multicellular = n_distinct(na.omit(paxdb_data$organ)) > 1
  whole_org = paxdb_data %>%
              dplyr::filter(is_integrated & organ == 'WHOLE_ORGANISM') %>%
              group_by(taxid,organ,protid,id_uniprot) %>%
              summarize( ppm_wholeorg = ppm) %>%
              group_by(taxid,organ) %>%
              mutate( pc_wholeorg = 100*percent_rank(-ppm_wholeorg) )

  organ_ppm = paxdb_data %>%
              dplyr::filter(is_integrated) %>%
              group_by(taxid,organ) %>%
              mutate( pc_organ = 100*percent_rank(-ppm), nprot_organ = n_distinct(protid))%>%
              group_by(taxid,protid,id_uniprot) %>%
              summarize(
                n_organs = paste0(sprintf("%s:%s",sort(organ),nprot_organ[order(organ)]),collapse="/"),
                ppm_organs = paste0(sprintf("%s:%.2f",sort(organ),ppm[order(organ)]),collapse="/"),
                pc_organs = paste0(sprintf("%s:%.2f",sort(organ),pc_organ[order(organ)]),collapse="/")
              )

  protein_ppm = paxdb_data %>%
    group_by(taxid,protid,id_uniprot) %>%
    summarize(
      ppm_md      = median_(ppm),
      ppm_avg     = mean_(ppm),
      ppm_gmean   = geomean(ppm),

      ppm_sd      = sd_(ppm),
      ppm_se      = sd_(ppm)/n(),
      ppm_int     = median_(ppm[is_integrated]),

      n_int       = sum_(is_integrated),
      norgan      = n_distinct(organ),
      norgan_int  = n_distinct(organ[is_integrated]),
      ndata       = n_distinct(dataset),

      ppm_max     = max_(ppm),
      ppm_int_max = max_(ppm[is_integrated]),
    ) %>%
    group_by(taxid) %>%
    mutate(
      pc_md       = 100*percent_rank(-ppm_md),
      pc_avg      = 100*percent_rank(-ppm_avg),
      pc_gmean    = 100*percent_rank(-ppm_gmean),
      pc_se       = 100*percent_rank(-ppm_se),
      pc_int      = 100*percent_rank(-ppm_int),

      pc_max      = 100*percent_rank(-ppm_max),
      pc_int_max  = 100*percent_rank(-ppm_int_max),
  ) %>%
  distinct %>%
  left_join(whole_org)

  if(multicellular){
    message('Multicelllular organim -> return integrated abundance in organs')
    return( protein_ppm %>% left_join(organ_ppm) )
  }
  return(protein_ppm)
}

#test = get.paxdb(4932,abundance = 'integrated')
# test.num = test %>% ungroup() %>% dplyr::select(-c(taxid,organ,protid)) %>% as.matrix
# C=cor(test.num,use='pairwise.complete',met='spearman')
# corrplot::corrplot(C,
#                    method = 'ellipse',
#                    diag = F, type = 'upper',
#                    addCoef.col = "gold1",
#                    number.cex = 0.7,number.digits=2,number.font=1)

get.ppm.ortho = function(node="4751.fungi", raw=F, which.abundance="integrated"){
  #closeAllConnections() # Make sure to close all connections before downloading data
  # Retrive orthologs from node
  message("Selecting orthologs for [",node,"]")
  ortho = load.paxdb.orthologs(node)
  # Obtain unique taxons
  taxons = sort(unique(ortho$taxid))
  ntax = length(taxons)
  message("node ",node," has ",ntax," taxons : ",toString(taxons))
  if(raw){
    message("---> Returning raw abundance values from individual datasets...")
    ppms = map_dfr(taxons, load.paxdb)
  }else{
    ppms = map_dfr(taxons, get.paxdb,which.abundance)
  }

  ortho_ppms = inner_join(ortho,ppms, by=c('taxid','protid')) %>%
      arrange(NOG,taxid,protid,id_uniprot)

  return (ortho_ppms)
}

# show.eggnog.nodes = function(ver=5,species=4891) {
#   library(rio)
#   library(RCurl)
#   library(hablar)
#   eggnog_nodes=sprintf("http://eggnog%s.embl.de/download/latest/e%s.level_info.tar.gz",ver,ver)
#   "http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv"
#   # download the archive
#   temp<-tempfile()
#   download.file(eggnog_nodes,temp)
#
#   # Find the list of files in the archive (taxonomic nodes)
#   header_col =  readr::read_delim(temp,delim='\t',skip=0) %>% grep(pattern="^#", x = ., value = T)
#
#   eggnog_levels = readr::read_delim(file = temp)$Name %>%
#     str_subset("\\.txt") %>%
#     basename %>%
#     str_remove(".orthgroups.consistent.txt") %>%
#     gtools::mixedsort()
#
#   if(show.nodes){ return(paxdb_nodes) }
#   # Check the taxonomic node selected
#   node_exists= !purrr::is_empty(node)
#   is_paxdb_node=F
#
#   if(node_exists){
#     numnode=grep(node,paxdb_nodes)
#     if(length(numnode)==1 ){
#       is_paxdb_node=T
#       node_full=grep(node,paxdb_nodes,v=T)
#       message("taxonomic node (",node_full,") is valid!")
#     }else if(length(numnode)>1){
#       warning("(",node,") is not unique! Please select a valid unique node",immediate.=T)
#       is_paxdb_node=F
#     }
#   }
#
#   if(!is_paxdb_node){
#     warning("(",node,") is not a valid taxonomic node from paxDB orthologs!",immediate.=T)
#     numnode = menu(paxdb_nodes, title = "Pick a taxonomic node from the list below",graphics = T)
#   }
# }

##### eggNOG #####
find_eggnog_downloads = function(){
  URL_EGGNOG = "http://eggnog.embl.de/download/latest/"
  eggnog_files = rvest::read_html(URL_EGGNOG) %>%
    rvest::html_nodes("a") %>%
    rvest::html_text(trim = T) %>%
    str_subset(pattern = "/$",negate = T)
  return(eggnog_files)
}

find_eggnog_version = function(.print=T){
  eggnog_files=find_eggnog_downloads()
  .version = str_extract(pattern="^(e[0-9](\\.[0-9])?)\\.",string = eggnog_files) %>%
             gtools::mixedsort(na.last = F) %>%
             last %>% str_sub(end=-2L) #remove the '.' after the version number
  if(.print){ cat(sprintf("EggNOG VERSION: %s\n",.version)) }
  return(.version)
}

find_eggnog_taxlevels = function(.print=T){
  URL_EGGNOG = "http://eggnog.embl.de/download/latest/"
  taxlevels = rvest::read_html(paste0(URL_EGGNOG,"per_tax_level/")) %>%
    rvest::html_nodes("a") %>% # Retrieve hyperlink corresponding to taxonomic level
    rvest::html_text(trim = T) %>% # Taxonomic level id is each hyperlink text
    str_subset(pattern = "^[0-9]+/$") %>% # Retrieve the taxonomic level (directories)
    str_sub(end=-2L) # Remove the "/" of the directories
  nnodes = n_distinct(taxlevels)
  if(.print){ .info$log(sprintf("found %s eggnog taxonomic levels...",nnodes)) }
  return(taxlevels)
}

eggnog_annotations_species=function(node,species){
  URL_EGGNOG = "http://eggnog.embl.de/download/latest/"
  .ver=find_eggnog_version(.print = F)
  eggnog_node = find_eggnog_node(node,.print=T)
  #library(rotl)
  eggnog_annotation_file = sprintf("%s/%s.og_annotations.tsv",URL_EGGNOG,.ver)
  .info$log('reading eggnog annotations...')
  eggnog_annotations_node = readr::read_delim(eggnog_annotation_file,"\t",
                                              col_types = 'icfc',
                                              col_names = c('nodes','og','letter','annotation')) %>%
                            dplyr::filter(nodes == eggnog_node$id)

  .info$log('reading eggnog orthogroups...')
  members_node = get_eggnog_node(eggnog_node$id,to_long = T,.print = F) %>%
        dplyr::select( -c(algo,tree,taxon_ids,url_fasta) ) %>%
        distinct()

  .info$log('merging annotations on orthogroups...')
  annotation_sp = left_join(members_node,eggnog_annotations_node, by=c('node'='nodes','OG'='og')) %>%
                    dplyr::filter(taxid %in% species)
  return(annotation_sp)
}

find_eggnog_node=function(node,GUI=F,.print=T){

  URL_EGGNOG = "http://eggnog.embl.de/download/latest/"
  taxlevels = find_eggnog_taxlevels(.print=F)

  eggnog_tax_info = sprintf("%s/%s.taxid_info.tsv",URL_EGGNOG, find_eggnog_version(.print=F))
  egg_tax = readr::read_delim(eggnog_tax_info,delim="\t",col_types = 'ccccc',progress = F,
                              skip = 1,col_names =  c('taxid','taxon','rank','lineage_name','lineage_id')) %>%
            mutate(lineage_name = str_replace_all(lineage_name,", ", replacement = "_"))

  taxlevel = egg_tax %>% dplyr::select(-taxid,-taxon,-rank) %>% distinct() %>%
             separate_rows("lineage_id","lineage_name",sep=',',convert = T) %>%
             dplyr::rename(id=lineage_id, name=lineage_name) %>%
             filter(id %in% taxlevels) %>%
             group_by(id) %>% add_count(name='size') %>%
             distinct() %>% arrange(id)

  tax_nodes = sprintf("%-s_%-s (%-s)",taxlevel$id, taxlevel$name, taxlevel$size)

  if(missing(node)){
    chosen_node =  menu(tax_nodes, graphics = GUI, title = 'choose a taxonomic node...')
    df_node = taxlevel[chosen_node,]
  }else if(any(node == taxlevel$id)){
    df_node = taxlevel[ taxlevel$id == node, ]
  }else if(any(grepl(node,taxlevel$name,ignore.case = T))){
    df_node = taxlevel[ grep(node,taxlevel$name,ignore.case = T), ]
  }else{
    warning(sprintf("Taxonomic level (%s) is not found in EggNOG database!",node))
    return(NA)
  }

  if(.print){
    .succ$log(sprintf('taxonomic level is : %s_%s (n=%s)',df_node$id, df_node$name, df_node$size))
  }

  return(df_node)
}

get_eggnog_species = function(node,.print = T){
  URL_EGGNOG = "http://eggnog.embl.de/download/latest/"
  eggnog_node = find_eggnog_node(node,.print=F)
  eggnog_tax_info = sprintf("%s/%s.taxid_info.tsv",URL_EGGNOG,find_eggnog_version(.print=F))

  if(.print){
    .info$log(sprintf('retrieving species for taxonomic level %s_%s...',eggnog_node$id, eggnog_node$name))
  }

  sp_info = readr::read_delim(eggnog_tax_info,delim="\t",col_types='ccccc',progress = F,
                              col_names =  c('taxid','taxon','rank','lineage_name','lineage_id')) %>%
            mutate(lineage_name = str_replace_all(lineage_name,", ", replacement = "_")) %>%
            # find species with node in their lineage
            dplyr::filter(grepl(eggnog_node$name,lineage_name)) %>%
            mutate(node_id=eggnog_node$id, node_name = eggnog_node$name, node_size = eggnog_node$size) %>%
            relocate(node_id,node_name,node_size)
  return(sp_info)
}

get_eggnog_taxonomy = function(node,.print=T,only_clade=T){

  taxlevel = find_eggnog_node(node,.print = F)
  node_sp = get_eggnog_species(node,.print = F)
  SP = node_sp$taxid
  TAXLEVELS = find_eggnog_taxlevels()
  if(.print){
    .info$log(sprintf("retrieving taxonomy for %s_%s...",taxlevel$id,taxlevel$name))
  }

  clades = node_sp %>%
           separate_rows(c('lineage_id','lineage_name'), sep=',') %>%
           mutate(seen=1) %>%
           group_by(node_id,node_name,node_size, clade_id=lineage_id, clade_name=lineage_name) %>%
           summarize(clade_size=sum(seen)) %>%
           arrange(desc(clade_size),clade_id,clade_name) %>%
           ungroup() %>%
           mutate(is_eggnog = (clade_id %in% TAXLEVELS),
                  is_clade = is_eggnog & !(clade_id %in% SP),
                  is_subnode = clade_size <= node_size) %>%
           relocate(node_id,node_name,node_size) %>%

  if(only_clade){
      clades_ = clades %>%
        filter(is_clade) %>%
        rowwise() %>%
        mutate( clade_sp = list(get_eggnog_species(clade_id,.print = F) %>% pull(taxid,taxon)))
    return(clades_)
  }

  return(clades)
}

fetch_eggnog_fasta = function(og,download=F){
  URL_FASTA_EGGNOG = "http://eggnogapi5.embl.de/nog_data/text/fasta"
  url_fasta_og = sprintf("%s/%s",URL_FASTA_EGGNOG,og)

  if(url.exists(url_fasta_og)){
    fasta_og = Biostrings::readAAStringSet(url_fasta_og)
    return(fasta_og)
  }else{
    .error$log(sprintf("Orthogroup not found (%s)",og))
    return(NULL)
  }
}

get_eggnog_alignment = function(node, use_trimmed=T, max_timeout=500){

  default_timeout = getOption("timeout")
  options(timeout = max(max_timeout, default_timeout))
  path_eggnog = here::here('data','eggnog')
  dir.create(path_eggnog,showWarnings = F,recursive = T)
  ncpus = parallel::detectCores()-2

  library(archive)
  URL_EGGNOG = "http://eggnog.embl.de/download/latest/"
  eggnog_node = find_eggnog_node(node,.print=T)
  find_eggnog_version()
  url_node_info= paste0(URL_EGGNOG,"per_tax_level/",eggnog_node$id,"/")

  .info$log(sprintf('get alignment for taxonomic node %s (%s)\n',eggnog_node$id,eggnog_node$name))

  eggnog_node_files = rvest::read_html(url_node_info) %>%
    rvest::html_nodes("a") %>% # Retrieve hyperlink corresponding to taxonomic level
    rvest::html_text(trim = T) # Taxonomic level id is each hyperlink text

  ali_tar = grep("_raw_algs.tar",eggnog_node_files,v=T)
  if(use_trimmed){
    .info$log('[using trimmed alignments]')
    ali_tar = grep("_trimmed_algs.tar",eggnog_node_files,v=T)
  }
  url_ali_node = paste0(url_node_info,ali_tar)

  eggnog_node_ali = file.path(path_eggnog,ali_tar)
  if( !file.exists(eggnog_node_ali) ){ # download cause it is large archive
    .info$log('download alignment archive for taxonomic node...')
    download.file(url = url_ali_node, eggnog_node_ali, cacheOK = T )
  }
  .info$log('extract alignments for taxonomic node...')
  ali_files = archive::archive(eggnog_node_ali) |>
              mutate( og = str_replace(string = basename(path), '\\..+gz$',''))
  archive::archive_extract(eggnog_node_ali, dir=path_eggnog)

  .info$log(sprintf('reading alignments in parallel (%s cores)...\n',ncpus))
  node_ali = pbmcapply::pbmclapply(X = file.path(path_eggnog,ali_files$path),
                                   FUN=Biostrings::readAAMultipleAlignment,
                                    mc.cores = ncpus)

  names(node_ali) = ali_files$og

  options(timeout = default_timeout)
  return(node_ali)
}

get_eggnog_node = function(node,to_long=F,.print=T){

  URL_EGGNOG = "http://eggnog.embl.de/download/latest/"
  URL_FASTA_EGGNOG = "http://eggnogapi5.embl.de/nog_data/text/fasta"

  find_eggnog_version(.print = .print)
  eggnog_node = find_eggnog_node(node,.print=.print) %>% rename_with(~paste0("node_",.))
  #library(rotl)
  url_node_info = paste0(URL_EGGNOG,"per_tax_level/",eggnog_node$node_id,"/")
  eggnog_node_files = rvest::read_html(url_node_info) %>%
                      rvest::html_nodes("a") %>% # Retrieve hyperlink corresponding to taxonomic level
                      rvest::html_text(trim = T) # Taxonomic level id is each hyperlink text

  members_file = paste0(url_node_info,grep("_members.tsv.gz",eggnog_node_files,v=T))
  trees_file = paste0(url_node_info,grep("_trees.tsv.gz",eggnog_node_files,v=T))
  library(ape)
  node_trees = readr::read_delim(file = trees_file , delim='\t',
                                 col_names=c('node','OG','algo','tree'),
                                 col_types='icfc',progress = F)

  node_members = readr::read_delim(members_file,delim = '\t',
                    col_names = c('node','OG','og_np','og_ns','string_ids','taxon_ids'),
                    col_types = 'iciicc',progress = F) %>%
                 left_join(eggnog_node, by=c('node'='node_id')) %>%
                 relocate(starts_with('node')) %>%
                 mutate(og_one2one = (og_np == og_ns) ) %>%
                 mutate(url_fasta = sprintf("%s/%s",URL_FASTA_EGGNOG,OG)) %>%
                 left_join(node_trees, by=c('node','OG'))

  if(to_long){
    if(.print){
      .info$log("returning long format *rows= unique(taxon + string_id) per orthogroup*")
    }
    node_members = node_members %>%
                   separate_rows(string_ids,sep=',',convert = T) %>%
                   separate('string_ids', c('taxid','string'), sep = '\\.', extra='merge',fill = 'right',remove = F)
  }

  return(node_members)
}

find_eggnog_subnode=function(node_clade,subnode){

  if(missing(node_clade)){
    .error$log("requires species information at a taxonomic level... (use get_eggnog_taxonomy(taxid))")
  }else if( length(node_clade)==1 && is_number(node_clade) ){
    node_clade = get_eggnog_taxonomy(node_clade)
  }

  node_subnodes = node_clade %>%
                  filter(is_clade &  is_subnode & clade_size > 1) %>%
                  mutate(clade_desc = sprintf("%s_%s (n=%s)",clade_id,clade_name,clade_size))

  nsub=length(subnode)
  df_subnode = node_subnodes %>% filter(clade_id %in% subnode )
  .info$log(sprintf('checking the (%s) following subnodes: %s',nsub,paste(subnode,collapse=" ")))

  subnode_found = n_distinct(df_subnode$clade_id)

  if( subnode_found < nsub ){
    na_subnode = setdiff(subnode,node_subnodes$clade_id)
    .warn$log(sprintf('node %s (%s) does not contain any of those %s subnodes: %s ',taxlevel$id, taxlevel$name,n_distinct(na_subnode),paste0(na_subnode,collapse=" ")))

    if(subnode_found == 0){
      clade_choices = node_subnodes %>% mutate( rownum = row_number() ) %>% pull(clade_desc,rownum)
      .error$log(sprintf("None of the input clades were found within the taxonomic level %s_%s !",taxlevel$id, taxlevel$name))
      choice = select.list(choices=clade_choices, multiple=T, graphics = F, title='pick a valid clade...')
      subnode_found = n_distinct(choice)
      df_subnode = node_subnodes[names(choice),]
    }
  }

  nid=unique(df_subnode$node_id)
  nname=unique(df_subnode$node_name)
  nfound = paste(df_subnode$clade_id,collapse=' ')
  .succ$log(sprintf('node %s (%s) contain those %s subnodes: %s ',nid,nname,subnode_found,nfound))

  return(df_subnode)
}

count_taxons_eggnog_node = function(node, subnode=1){
  # e.g.
  # node=4751
  # subnode=c(4890,5204,451866,4891,147541,147545,147550,147548)

  taxlevel = find_eggnog_node(node)
  df_subnode = find_eggnog_subnode(taxlevel$id,subnode)

  node_clades = get_eggnog_taxonomy(taxlevel$id)
  node_taxons  = get_eggnog_species(taxlevel$id)
  node_species = node_taxons %>% pull(taxid,taxon)

  .info$log('count number of orthologs/species in orthogroups...')
  #tictoc::tic('count number of orthologs/species in orthogroups...')
  node_members = get_eggnog_node(taxlevel$id,.print = F,to_long = T) %>%
                 dplyr::rename(node_id=node) %>%
                 left_join(node_taxons, by='taxid')


  node_orthologs = node_members %>%
    dplyr::select(starts_with('node'),OG,taxid,string) %>%
    group_by(OG) %>%
    mutate( taxid = factor(taxid,node_species),
            id = paste0(taxid,".",string),
            node_orthogroup = list(id)) %>%
    nest( node_orthologs = c(taxid,string,id) )
  #tictoc::toc()

  .info$log('compute orthogroups statistic for the taxonomic level...')
  og_info =  node_members %>%
             dplyr::select( starts_with('node'), OG, algo, url_fasta, tree) %>%
             distinct()

  og_stats = node_members %>%
             dplyr::select( starts_with(c('node','og'),ignore.case=T), taxid,string ) %>%
             distinct() %>%
             group_by(OG,taxid) %>%
               add_count(name='og_northo') %>%
             group_by(OG) %>%
               mutate(og_n1to1 = sum(og_northo == 1)) %>%
             dplyr::select(-taxid,-string,-og_northo) %>%
             distinct()

  df_og = left_join(og_info,og_stats,by = c("node_id", "node_name", "node_size", "OG")) %>%
          left_join(node_orthologs,by = c("node_id", "node_name", "node_size", "OG"))

  clade_list = list()
  for(i in 1:nrow(df_subnode)){

    df_clade = df_subnode[i,]
    .info$log(sprintf('%2d/ count species/orthologs for clade  %s_%s (n=%s)...',i,df_clade$clade_id, df_clade$clade_name, df_clade$clade_size))
    subnode_species = get_eggnog_species(df_clade$clade_id) %>% pull(taxid,taxon)

    #tic('clade ortholog')
    clade_orthologs = left_join(df_clade,node_members,by = c("node_id", "node_name", "node_size")) %>%
                      dplyr::select(starts_with(c('node','clade'),ignore.case = T),OG,taxid,string) %>%
                      filter( taxid %in% subnode_species ) %>%
                      group_by(OG,taxid) %>%
                      add_count(name='clade_northo') %>%
                      group_by(OG) %>%
                      mutate( taxid = factor(taxid,node_species),
                              id = paste0(taxid,".",string),
                              clade_ns = n_distinct(taxid),
                              clade_np = n_distinct(string),
                              clade_one2one = clade_ns==clade_np,
                              clade_f  = clade_ns / clade_size,
                              clade_n1to1  = sum(clade_northo == 1),
                              clade_f1to1 = clade_n1to1 / clade_ns
                      ) %>%
                      dplyr::select(-clade_northo) %>%
                      distinct() %>%
                      mutate( clade_orthogroup = list(id) ) %>%
                      nest( clade_orthologs = c(taxid,string,id) )
    #toc()

    clade_list[[df_clade$clade_desc]] = left_join(df_og,clade_orthologs, by = c("node_id", "node_name", "node_size", "OG"))
  }
  if( nrow(df_subnode) == 1 ){ return( unlist(clade_list) )}
  return(clade_list)
}

find.common.ancestor= function(lineage){
  L = strsplit(lineage,',')
  LCA = Reduce(intersect, L)
  MRCA = sapply(L,function(x){ setdiff(unlist(but.last(x)), LCA) })
  return(MRCA)
}

##### ELM (Eukaryotic Linear Motifs) #####

get_elm_motifs = function(){
  library(readr)
  URL_ELM = "http://elm.eu.org/"
  url_motifs = sprintf("%s/elms/elms_index.tsv",URL_ELM) # MOTIF REGEX
  elm_motifs = read_delim(url_motifs,delim = "\t",comment = "#")
  elm_cols = c('elm_acc','elm_id','elm_name','elm_desc','elm_regex','elm_prob','n_instances','n_pdb_instances')

  motifs = elm_motifs %>% dplyr::rename( set_names( names(elm_motifs), elm_cols ) ) %>%
           mutate(elm_func=str_split_fixed(elm_id,pattern = '_',n = 2)[,1])
  return(motifs)
}

get_elm_instances = function(){
  library(readr)
  URL_ELM = "http://elm.eu.org/"
  url_instances = sprintf("%s/instances.tsv?q=*&taxon=&instance_logic=",URL_ELM) # MOTIF OCCURRENCES
  elm_instances = read_delim(url_instances,delim = "\t",comment = "#")
  elm_cols = c('elmi_acc','elm_func','elm_id','uni_name','uni','uni_acc','elm_start','elm_end','references','method','elm_conf','pdb','org')

  instances = elm_instances %>%
              dplyr::rename( set_names( names(elm_instances), elm_cols ) ) %>%
              mutate(org=str_wrap(org,20)) %>%
              relocate('org','uni','uni_name','uni_acc',
                       'elmi_acc','elm_id','elm_func','elm_start','elm_end')
  return(instances)
}

get_elm_interactions = function(){
  library(readr)
  URL_ELM = "http://elm.eu.org/"
  url_interactions = sprintf("%s/interactions/as_tsv",URL_ELM) # MOTIF INTERACTIONS
  elm_interactions = read_delim(url_interactions, col_types='cccciiiinnccc_', delim = "\t") # fix column format and skip the last column (14th column)
  elm_cols = c("elm_id","domain_id","elm_interactor","domain_interactor",'elm_start','elm_end',
               'domain_start','domain_end',"Kd_min","Kd_max","PMID","elm_tax",'domain_tax')

  interactions = elm_interactions %>%
                 dplyr::rename( set_names( names(elm_interactions), elm_cols ) ) %>%
                 relocate(Kd_min,Kd_max,
                          elm_tax,elm_id,elm_interactor,elm_start,elm_end,
                          domain_tax,domain_id,domain_interactor,domain_start,domain_end) %>%
                 separate(col=domain_tax, into=c("domain_taxid",'domain_species',NA), sep="[\\(\\)]") %>%
                 separate(col=elm_tax, into=c("elm_taxid",'elm_species',NA), sep="[\\(\\)]")

  return(interactions)
}

get_elm = function(){

  ELM = left_join(get_elm_instances(),get_elm_motifs(),by=c('elm_id')) %>%
        left_join(interactions, by=c("elm_id","uni"="elm_interactor",'elm_start','elm_end')) %>%
        dplyr::select(-method,-pdb,-PMID,-n_pdb_instances,-references)

  return(ELM)
}


# Reference sequences ----------------------------------------------------------
##### SGD #####
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

##### Uniprot #####

get_uniprot_id = function(accession){
  UNIPROT_URL = sprintf("https://rest.uniprot.org/uniprotkb/%s.tsv",accession)
  #OX   NCBI_TaxID=7955 {ECO:0000312|Proteomes:UP000000437};
  if( httr::http_error(UNIPROT_URL) ){ return(NULL) }
  res = readr::read_delim(UNIPROT_URL,show_col_types = FALSE, progress = F)
  return(res)
}

get_uniprot_ids = function(accessions){
  msg = sprintf("Retrieve %s uniprot ids from URL...",n_distinct(accessions))
  if(require(pbmcapply)){
    ncpus=parallel::detectCores()-1
    msg = paste0(msg, sprintf("\n using %s CPUs in parallel",ncpus))
    message(msg)
    df_uniprot = pbmcapply::pbmclapply(X=accessions, get_uniprot_id, mc.cores = ncpus) %>%
      purrr::compact() %>%
      bind_rows()
  }else{
    message(msg)
    warning('NOT IN PARALLEL (might be long depending on number of IDS)')
    df_uniprot = lapply(accessions,get_uniprot_id) %>% purrr::compact() %>% bind_rows()
  }
  return(df_uniprot)
}

parse_uniprot_fasta_header = function(fasta_header){
#  uni_desc = get.uniprot.proteome(9606,DNA = F,fulldesc = T) %>% names

  #IDS = str_extract(uni_desc,pattern="^[^ ]+")
  df_ids = str_split_fixed(fasta_header,pattern = '[\\| ]', n=4) %>%
            set_colnames(c('DB','AC','ID','DESC')) %>%
            as_tibble()

  df_name = str_split_fixed(df_ids$DESC,' OS=',n = 2) %>%
             set_colnames(c('NAME','DESC')) %>%
             as_tibble() %>%
             mutate(DESC = paste0("OS=",DESC))

  df_info = tibble(
    OS = str_extract("(?<=OS\\=)([^\\=]+)(?= OX\\=)", string = df_name$DESC),
    OX = str_extract("(?<=OX\\=)([^\\=]+)(?= (GN|PE)\\=)", string = df_name$DESC),
    GN = str_extract("(?<=GN\\=)([^\\=]+)(?= PE\\=)", string = df_name$DESC),
    PE = str_extract("(?<=PE\\=)([^\\=]+)(?= SV\\=)", string = df_name$DESC),
    SV = str_extract("(?<=SV\\=)([^\\=]+)$", string = df_name$DESC)
  )

  df_uni_desc = bind_cols(df_ids,df_name,df_info) %>%
                dplyr::select(-starts_with("DESC")) %>% type_convert() %>%
                relocate(DB,OX,OS,AC,ID,GN,NAME,PE,SV)
  return(df_uni_desc)
}

get_uniprot_reference = function(taxon=9606){
  cdna_name = get.uniprot.proteome(taxon,DNA=T,fulldesc=T) %>% names
  prot_name = get.uniprot.proteome(taxon,DNA=F,fulldesc=T) %>% names

  id_dna = str_split_fixed(string=cdna_name,pattern = '\\|',n = 3) %>%
    set_colnames(c('db','uni','id_cdna')) %>% as_tibble %>%
    mutate(ensp = str_extract(id_cdna,ENSEMBL.nomenclature()))
  id_prot = prot_name %>% parse_uniprot_fasta_header()

  # PROTEOME first
  id_reference =left_join(id_prot,id_dna,by=c('DB'='db','AC'='uni'))
  return(id_reference)
}

get.uniprot.mapping = function(taxid, targetdb=NULL) {
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
  if(!is.null(targetdb)){
    all_db = sort(unique(mapped$extdb))
    dbs = match.arg(targetdb,choices = all_db, several.ok = T)
     mapped = mapped %>% dplyr::filter( extdb %in% dbs )
  }
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
    if(length(keyword)==1){
      matched = get_rows_by_keyword(word = keyword, df = refprot)
      message(sprintf('%s entries matched keyword "%s"',nrow(matched),keyword))
    }else{
      matched = find_keywords(refprot,keyword,strict=T)
      message(sprintf('%s entries matched keywords "%s"',nrow(matched),toString(unique(keyword))))
    }
    if(all){
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

  if(missing(taxid)){  stop("Need an uniprot taxon id") }
  UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
  SEQTYPE = ".fasta.gz"

  refprot = find.uniprot_refprot(all=T)
  found = refprot$tax_id %in% taxid
  if(!any(found)){ stop(sprintf("%s not found in the reference proteome!",taxid)) }
  TAX = str_to_title(refprot$superregnum[which(found)])
  UPID = refprot$proteome_id[which(found)]
  if(DNA){
    SEQTYPE = "_DNA.fasta.gz"
    genome_url = sprintf("%s/%s/%s/%s_%s%s",UNIPROT_URL,TAX,UPID,UPID,taxid,SEQTYPE)
    UNI = load.genome(genome_url)
  }else{
    proteome_url = sprintf("%s/%s/%s/%s_%s%s",UNIPROT_URL,TAX,UPID,UPID,taxid,SEQTYPE)
    UNI = load.proteome(proteome_url)
  }
  if(!fulldesc){
    regexUNIPROTAC = "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
    names(UNI) = str_extract(names(UNI), regexUNIPROTAC)
  }
  return(UNI)
}

load.uniprot.proteome = function(species='yeast') { # Older version of get.uniprot.proteome
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

query_uniprot_subloc = function(uniprot, taxon, MAX_QUERY=200, todf=T){
  #.org=''
  UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/search?query="

  .accession=''
  if( missing(uniprot) & missing(taxon) ){
    stop('Query requires  a taxon id (also accepts an optional list of uniprot)...')
  }else if(missing(uniprot) & !missing(taxon)){
    uniprot = get.uniprot.proteome(taxon,DNA=F) %>% names %>% unique
    #.org = paste0('organism_id:',taxon,'&')
    message(sprintf('Using the uniprot acession from taxon %s \n',taxon))
  }

  n_uni = n_distinct(uniprot)
  n_queries = ceiling(n_uni / MAX_QUERY)

  queries = character(length = n_queries)
  for(i in 1:n_queries){
    i0 = (i-1)*MAX_QUERY+1
    num = seq(i0,len=MAX_QUERY)
    uni_query = uniprot[num] %>% na.omit() %>% as.vector
    .accession = paste0("%28",str_c('accession:',uni_query,collapse='+OR+'),"%29&")
    url_query = sprintf('%s%sfields=accession,cc_subcellular_location&size=%s&format=tsv',UNIPROT_REST,.accession,MAX_QUERY)
    queries[i] =url_query
  }

  run_uniprot_tsvquery = function(tsv_query){
    #print(tsv_query)
    query_res = read_delim(tsv_query,delim = '\t')
    return(query_res)
  }


  tictoc::tic('retrieve subcellular locations from uniprotKB...')
  message(sprintf("Retrieving subcellular locations from a list of uniprot (n=%s)...",n_uni))
  if (require('pbmcapply')) {
    ncores=parallelly::availableCores(which='max')-1
    message(sprintf("using 'pbmcapply' in parallel with %s cpus",ncores))
    uni_subloc = pbmcapply::pbmclapply(X=queries,FUN=run_uniprot_tsvquery, mc.cores=ncores, mc.silent=F, mc.cleanup = T)
  }else{
    message("NOT PARALLEL! This may take a while (>20mn)")
    uni_subloc = lapply(queries, run_uniprot_tsvquery)
  }
  tictoc::toc()

  df_subloc = uni_subloc %>%
    bind_rows() %>%
    dplyr::rename('UNIPROT'='Entry','SUBCELLULAR_LOCATION'='Subcellular location [CC]') %>%
    dplyr::filter(!is.na(SUBCELLULAR_LOCATION))

  subloc = df_subloc %>% deframe()
  uni.loc = parse.uniprot.subcellular_locations(subloc)


  if( todf ){
    uni.isloc  = uni.loc %>%  mutate(seen=1) %>%
      pivot_wider( id_cols = c('id','HAS_FOCI','HAS_ISOFORM'),
                   names_from = 'loc', values_from = 'seen',
                   values_fn=list(seen = sum),values_fill = list(seen=0))
    return(uni.isloc)
  }
  return(uni.loc)
}

##### Ensembl #####
ENS_MIRROR='asia'

get_ensembl_version = function(latest=T,withURL=T,URL_FTP_ENSEMBL="http://ftp.ensembl.org/pub/"){
  ensembl_releases = rvest::read_html(URL_FTP_ENSEMBL) %>%
    rvest::html_nodes("a") %>%
    rvest::html_text(trim = T) %>%
    str_subset('release-') %>%
    str_sub(end = -2) # remove backslash at the end

  if(withURL){ ensembl_releases = paste0(URL_FTP_ENSEMBL,ensembl_releases) }
  latest_release = gtools::mixedsort(ensembl_releases)[1]
  if(latest){ return(latest_release) }
  return(ensembl_releases)
}

get_ensembl_species= function(){
  URL_ENSEMBL = "https://www.ensembl.org/info/data/ftp/index.html"
  library(rvest)
  ens_ftp <- URL_ENSEMBL %>% read_html()
  # single species data table

  link2data = ens_ftp %>% html_element(css='.data_table') %>% html_elements("a") %>% html_attr("href")
  length(link2data)
  url_ftp_pub='http://ftp.ensembl.org/pub/'

  ss_table <- ens_ftp %>% html_element(css='.data_table') %>% html_table() %>%
    janitor::clean_names() %>%
    mutate(org = str_replace(species,'^(.+[a-z])[A-Z].+$','\\1'),
           spname  = str_replace(species,'^.+[a-z]([A-Z].+$)','\\1')) %>%
    dplyr::select(-x,-species) %>% relocate(org,spname)

  # COULD NOT EXTRACT URL LINKS AS A TABLE (so manually reconstructed links from ftp url and species name)
  ss_table_links = ss_table %>%
       mutate( species = tolower(ss_table$spname) %>% str_replace_all(' ','_')  )%>%
       mutate( dna_fasta = sprintf('%s/%s/%s/%s',url_ftp_pub,'current_fasta',species,'dna'),
               c_dna_fasta = sprintf('%s/%s/%s/%s',url_ftp_pub,'current_fasta',species,'cdna'),
               cds_fasta = sprintf('%s/%s/%s/%s',url_ftp_pub,'current_fasta',species,'cds'),
               nc_rna_fasta = sprintf('%s/%s/%s/%s',url_ftp_pub,'current_fasta',species,'ncrna'),
               protein_sequence_fasta = sprintf('%s/%s/%s/%s',url_ftp_pub,'current_fasta',species,'pep'),
               annotated_sequence_embl = sprintf('%s/%s/%s/',url_ftp_pub,'current_embl',species),
               annotated_sequence_gen_bank = sprintf('%s/%s/%s/',url_ftp_pub,'current_genbank',species),
               gene_sets_gtf = sprintf('%s/%s/%s/',url_ftp_pub,'current_gtf/',species),
               gene_sets_gff3 =sprintf('%s/%s/%s/',url_ftp_pub,'current_gff3/',species),
               other_annotations_tsv = sprintf('%s/%s/%s/',url_ftp_pub,'current_tsv',species),
               other_annotations_rdf = sprintf('%s/%s/%s/',url_ftp_pub,'current_rdf',species),
               other_annotations_json = sprintf('%s/%s/%s/',url_ftp_pub,'current_json',species),
               whole_databases=sprintf('%s/%s/%s/',url_ftp_pub,'current_mysql',species),
               variation_gvf=sprintf('%s/%s/%s/%s',url_ftp_pub,'current_variation',species,'/gvf'),
               variation_vcf=sprintf('%s/%s/%s/%s',url_ftp_pub,'current_variation',species,'/vcf'),
               variation_vep=sprintf('%s/%s/%s/%s',url_ftp_pub,'current_variation',species,'/vep'),
               regulation_gff=sprintf('%s/%s/%s/',url_ftp_pub,'current_regulation',species),
               data_files=sprintf('%s/%s/%s/',url_ftp_pub,'current_data_files',species),
               bam_big_wig=sprintf('%s/%s/%s/',url_ftp_pub,'current_bamcov',species)) %>%
    dplyr::select(-other_annotations,-gene_sets)

  return(ss_table_links)
}

find_ensembl_sptree = function(treename="", URL_SPTREE="http://ftp.ensembl.org/pub/current_compara/species_trees/"){

  #URL_FTP_ENSEMBL = host
  #"http://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi/current/compara/species_trees/fungi_protein-trees_default.nh"
  #URL_SPTREE = paste0(URL_FTP_ENSEMBL,)
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  httr::set_config(httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"))

  sptrees = rvest::read_html(URL_SPTREE) %>%
    rvest::html_elements("a") %>%
    rvest::html_text(trim = T) %>%
    str_subset(pattern = "\\.nh$")

  file_sptrees = basename(sptrees)
  if( is.null(treename) ){ treename = "" }
  url_tree = str_subset(file_sptrees, pattern = treename)

  if( length(url_tree) != 1 ){
    sptree_choice =menu(title = 'choose an Ensembl species tree file...',graphics = F,choices = file_sptrees)
    url_tree = sptrees[sptree_choice]
  }

  return(paste0(URL_SPTREE,url_tree))
}

get_ensembl_sptree = function(treename=NULL,URL_SPTREE="http://ftp.ensembl.org/pub/current_compara/species_trees/"){

  url_tree = find_ensembl_sptree(treename,URL_SPTREE) %>% utils::URLencode()
  treefile = basename(url_tree) %>% fs::path_ext_set('.nh')
  yeastomics_tree = here::here('data','ensembl',treefile)

  if( !file.exists(yeastomics_tree) ){
    dir.create(dirname(yeastomics_tree),showWarnings = F,recursive = T)
    download.file(url_tree,yeastomics_tree)
  }

  sptree=ape::read.tree(yeastomics_tree)

  return(sptree)
}

get_ensembl_19mammals = function(){

  mammals = tribble(~taxid, ~ens_pre, ~spname, ~org,
        9606,'ENSG', 'Homo sapiens', 'human',
        9913,'ENSBTAG', 'Bos taurus', 'cow',
        9615,'ENSCAFG', 'Canis lupus familiaris', 'dog',
        9685,'ENSFCAG', 'Felis catus', 'cat',
        9785,'ENSLAFG', 'Loxodonta africana', 'elephant',
        9598,'ENSPTRG', 'Pan troglodytes', 'chimpanzee',
        9739,'ENSTTRG','Tursiops truncatus','dolphin',
        9796,'ENSECAG','Equus caballus','horse',
        9595,'ENSGGOG','Gorilla gorilla gorilla','gorilla',
        30608,'ENSMICG','Microcebus murinus','mouse_lemur',
        9669,'ENSMPUG','Mustela putorius furo','ferret',
        10090,'ENSMUSG','Mus musculus','mouse',
        132908,'ENSPVAG','Pteropus vampyrus','bat',
        13616,'ENSMODG','Monodelphis domestica','opossum',
        9986,'ENSOCUG','Oryctolagus cuniculus','rabbit',
        9601,'ENSPPYG','Pongo abelii','orang_outan',
        10116,'ENSRNOG','Rattus norvegicus','rat',
        9823,'ENSSSCG','Sus scrofa','pig',
        43179,'ENSSTO','Ictidomys tridecemlineatus','squirrel'
  ) %>%
    mutate(sp=str_replace(spname,'^([A-Z]).+ ([a-z]+)','\\1\\2') %>% str_to_lower())
  return(mammals)
}

get_ensembl_sp = function(full_name, sep='_'){
  library(stringr)
  library(purrr)
  nsep = str_count(full_name,pattern = sep)
  spname = str_split(full_name, sep)

  long  = map_chr(spname, ~last(.x) %>% str_to_lower())
  short = map_chr(spname, ~head(.x,-1) %>% str_sub(start = 1,end=1) %>% str_to_lower() %>% concat())

  return(paste0(short,long))
}

get_ensembl_vertebrates=function(){
  library(readr)
  library(tidyverse)
  url_ens_vertebrates=paste0(get_ensembl_version(),"/species_EnsemblVertebrates.txt")
  url_uni_vertebrates=paste0(get_ensembl_version(),"/uniprot_report_EnsemblVertebrates.txt")

  ens_vertebrates=read_delim(url_ens_vertebrates) %>% mutate(species_id=str_replace_all(species_id,'\t',''))
  uni_vertebrates=read_delim(url_uni_vertebrates) %>% mutate(uniprotCoverage=str_replace_all(uniprotCoverage,'\t',''))

  vertebrates = inner_join(ens_vertebrates,uni_vertebrates,
                           by = c("#name", "species", "division", "taxonomy_id", "genebuild",
                                  'assembly'='assembly_name','assembly_accession'='assembly_id') ) %>%
                readr::type_convert()

  colnames(vertebrates) = c('organism','species','division','tax_id',
                            'assembly_version','assembly_accession','genebuild',
                            'variation','microarray','pan_compara','peptide_compara','genome_alignments',
                            'other_alignments', 'cored_db','species_id','n_protein_coding','n_swissprot','n_trembl','coverage')

  vertebrates$sp = get_ensembl_sp(vertebrates$species)

  return(vertebrates)
}

get_ensembl_fungi=function(){
  library(readr)
  library(tidyverse)
  URL_FTP_FUNGI = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi/current/"
  url_ens_fungi=paste0(get_ensembl_version(URL_FTP_ENSEMBL = URL_FTP_FUNGI),"/species_EnsemblFungi.txt")
  url_uni_fungi=paste0(get_ensembl_version(URL_FTP_ENSEMBL = URL_FTP_FUNGI),"/uniprot_report_EnsemblFungi.txt")

  ens_fungi=read_delim(url_ens_fungi) %>% mutate(species_id=str_replace_all(species_id,'\t',''))
  uni_fungi=read_delim(url_uni_fungi) %>% mutate(uniprotCoverage=str_replace_all(uniprotCoverage,'\t',''))

  fungi = inner_join(ens_fungi,uni_fungi,
                     by = c("#name", "species", "division", "taxonomy_id", "genebuild",
                            'assembly'='assembly_name','assembly_accession'='assembly_id') ) %>%
          readr::type_convert()

  colnames(fungi) = c('organism','species','division','tax_id',
                      'assembly_version','assembly_accession','genebuild',
                      'variation','microarray','pan_compara','peptide_compara','genome_alignments',
                      'other_alignments', 'cored_db','species_id','n_protein_coding','n_swissprot','n_trembl','coverage')

  fungi$sp = get_ensembl_sp(fungi$species)

  return(fungi)
}


###### biomart ensembl #####
get_ensembl_biomart = function(mart=NULL,quiet=T){
  library(biomaRt)
  library(dplyr)
  genomes = listEnsemblGenomes(includeHosts = T) %>% add_column(type='EnsemblGenomes')
  ensembl = listMarts(includeHosts = T,verbose = !quiet) %>% as_tibble %>% add_column(type='Ensembl')
  ensmarts = bind_rows(genomes,ensembl) %>%
                mutate(vnum=str_extract(database,'[0-9]+'),
                       host = paste0('https://',host),port = 443) # enforce  https
  if(is.null(mart)){
    m=menu(ensmarts$version,title = 'pick an ensembl biomart...')
    mart = ensmarts$biomart[m]
  }
  MART = ensmarts %>% filter( biomart == mart )

  BIOMART = useMart(biomart = MART$biomart, host = MART$host, path = MART$path, port = MART$port)
  return(BIOMART)
}

find_ensembl_datasets = function(BIOMART,quiet=T){
  library(biomaRt)
  if(missing(BIOMART)){ BIOMART = get_ensembl_biomart() }

  # host='https://www.ensembl.org/',mart='ensembl'
  # ens <- useEnsembl(biomart = mart, host = host, mirror=ENS_MIRROR)
  ens_dataset = listDatasets(BIOMART,verbose=!quiet) %>% dplyr::as_tibble() %>%
    dplyr::mutate(
      sp=str_split_fixed(dataset,'_',n=3)[,1],
      org = str_split_fixed(description,' genes ',n=2)[,1] %>% str_trim
    )
  return(ens_dataset)
}

get_ensembl_dataset = function(mart,organism){
  if(missing(mart)){ mart = get_ensembl_biomart() }
  mart = get_ensembl_biomart(mart)
  BM_org = find_ensembl_datasets(mart)

  if(missing(organism) ){
    organism = select.list(choices= BM_org$description,
                  title = 'pick an ensembl dataset from selected biomart...')
  }

  has_dataset = str_which(BM_org$dataset,pattern=fixed(organism,ignore_case = T))
  dataset_name = BM_org$dataset[ has_dataset ]
  no_dataset = length(has_dataset) == 0L
  if(no_dataset){
    has_organism = str_which(BM_org$description,pattern=fixed(organism,ignore_case = T))
    norg = length(has_organism)
    if(norg == 0L){
      .error$log(paste0(organism," dataset: not found in the selected mart!"))
      return(NA)
    }else if(norg > 1L){
      .warn$log(paste0("More than one organism matched :",BM_org$description[has_organism]))
      has_organism=has_organism[1]
      .warn$log(paste0("using the first one :",BM_org$description[has_organism]))
    }
    dataset_name = BM_org$dataset[ has_organism ]
  }
  .succ$log(paste0(organism," dataset: ",dataset_name))
  DATASET = useDataset(dataset_name, mart)
  #message(DATASET)
  return(DATASET)
}

get_ensembl_prot = function(BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human'),
                            verbose=T){
  library(biomaRt)
  att_prot = c('ensembl_peptide_id','uniprotswissprot')
  att_struct = c('gene_biotype')

  # Get representative human proteome with UniProt/SwissProt identifiers
  hs_ensp = getBM(mart = BIOMART,
                  attributes = c(att_prot,att_struct),
                  filters=c('biotype','transcript_biotype','with_uniprotswissprot'),
                  values=list('protein_coding','protein_coding',T),
                  uniqueRows = T, bmHeader = F) %>% as_tibble()
  nr = nrow(hs_ensp)
  np = n_distinct(hs_ensp$ensembl_peptide_id)
  nu = n_distinct(hs_ensp$uniprotswissprot)
  if(verbose)
    message(sprintf('rows = %7s | proteins = %6s | uniprot = %6s',nr,np,nu))
  return(hs_ensp)
}

get_ensembl_tx = function(verbose=T,
                          BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human'),
                          ENSG,ENSP){
  library(biomaRt)
  #ens <- biomaRt::useEnsembl(biomart = mart, host = host, dataset = dataset, mirror=ENS_MIRROR)
  #ens=biomaRt::useMart(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
  att_gene = c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id')
  att_struct = c('start_position','end_position','cds_length',
                 'transcript_length','transcript_start','transcript_end',
                 'ensembl_exon_id','rank','exon_chrom_start','exon_chrom_end','is_constitutive')
  filters = c('ensembl_gene_id'=list(ENSG) )

  #listAttributes(BIOMART)
  ens_canonical =  getBM(mart = BIOMART,
                     attributes = c(att_gene,'transcript_is_canonical'),
                     filters=c('ensembl_gene_id'),values=list(ENSG),
                     uniqueRows = T, bmHeader = F) %>% as_tibble() %>%
                     dplyr::rename(ensg=ensembl_gene_id, enst=ensembl_transcript_id, ensp=ensembl_peptide_id,
                                   is_canonical=transcript_is_canonical) %>%
                     mutate(is_ensgref = ensg %in% ENSG, is_enspref = ensp %in% ENSP,
                            is_canonical=replace_na(is_canonical,0),
                            has_ensp = (ensp!=""),
                            canonical = factor(1*(is_ensgref & is_enspref) + 2*(is_canonical>0),ordered = T )) %>%
                     group_by(ensg) %>% mutate(has_enspref = factor(sum(is_enspref)) ) %>%
                     ungroup() %>%
                     arrange(ensg,ensp) %>%
                     filter(is_enspref | canonical>=2 & has_enspref==0 )

  # Get representative of human transcriptome
  ensg_trans = getBM(mart = BIOMART,
                  attributes = c(att_gene,att_struct),
                  filters='ensembl_gene_id',values=list(ENSG),
                  uniqueRows = T, bmHeader = F) %>% as_tibble() %>%
                dplyr::rename(ensg=ensembl_gene_id, enst=ensembl_transcript_id, ensp=ensembl_peptide_id) %>%
                ungroup()

  ensp_trans = getBM(mart = BIOMART,
                     attributes = c(att_gene,att_struct),
                     filters='ensembl_peptide_id',values=list(ENSP),
                     uniqueRows = T, bmHeader = F) %>% as_tibble() %>%
                dplyr::rename(ensg=ensembl_gene_id, enst=ensembl_transcript_id, ensp=ensembl_peptide_id) %>%
                ungroup()

  ens_trans = left_join(ens_canonical,ensg_trans) %>% left_join(ensp_trans) %>%
              ungroup() %>% distinct() %>%
              mutate( gene_len = end_position-start_position+1,
                      exon_len = exon_chrom_end-exon_chrom_start+1 ) %>%
              dplyr::rename(cds_len=cds_length, transcript_len=transcript_length) %>%
              group_by(ensg,enst,ensp) %>%
              mutate(n_exons = n_distinct(ensembl_exon_id), n_exons_mini = sum(is_constitutive),
                     has_introns = n_exons_mini>1,
                     tot_exon_len = sum(exon_len) ) %>%
              group_by(ensg) %>%
              mutate(n_transcripts = n_distinct(enst),
                     n_proteins = n_distinct(ensp)) %>%
              dplyr::select(-start_position,-end_position,-transcript_start,-transcript_end,
                            -ensembl_exon_id,-exon_chrom_start,-exon_chrom_end,-tot_exon_len,-exon_len,-rank,-is_constitutive,
                            -is_ensgref,-has_enspref) %>%
              ungroup() %>% distinct()

  nr = nrow(ens_trans)
  ng = n_distinct(ens_trans$ensg)
  nt = n_distinct(ens_trans$enst)
  np = n_distinct(ens_trans$ensp)
  if(verbose)
    message(sprintf('rows = %7s | genes = %6s | transcripts = %6s | proteins = %6s',nr,ng,nt,np))
  return(ens_trans)
}

get_ensembl_gc = function(ENSG,BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human')){
  library(biomaRt)
  #ens=biomaRt::useMart(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')

  filters = c('ensembl_gene_id'=list(ENSG))
  hs_gc_gene=getBM(attributes=c("ensembl_gene_id",'percentage_gene_gc_content'), mart=BIOMART,
                 filters=names(filters), values=as.list(filters),
                 uniqueRows = T, bmHeader = F) %>% as_tibble() %>%
             dplyr::rename(ensg=ensembl_gene_id) %>% distinct()
  return(hs_gc_gene)
}

get_ensembl_chr = function(as.df=T,remove_patches=T,ENSG,
                           BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human')){
  library(biomaRt)

  chromosomes = c(1:22,'MT','X','Y')
  #ens=biomaRt::useMart(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
  att_gene = c('ensembl_gene_id','chromosome_name')
  filters = c('ensembl_gene_id'=list(ENSG))

  hs_chr=getBM(attributes=att_gene, mart=BIOMART,
               filters=names(filters), values=as.list(filters),
               uniqueRows = T, bmHeader = F) %>% as_tibble() %>%
    mutate(is_patched = !(chromosome_name %in% chromosomes)) %>%
    dplyr::rename(ensg=ensembl_gene_id)  %>% distinct()

  if(remove_patches){
    hs_chr = hs_chr %>% dplyr::filter(!is_patched)  %>% distinct()
  }

  patched_chr = setdiff(unique(hs_chr$chromosome_name),chromosomes)

  if(as.df){
    hs_chr= hs_chr %>%
      mutate( chr_val=T,
              chromosome_name = fct_other(chromosome_name, keep=chromosomes, other_level = "other")) %>%
      pivot_wider(id_cols=c('ensg'),
                  names_from=chromosome_name, names_prefix='chr_',
                  values_from = chr_val, values_fill = F )  %>% distinct()

  }
  return(hs_chr)
}


get_hs_GC = function(with_uniprot=T,
                     BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human')){
  library(biomaRt)
  att_gene = c('ensembl_gene_id','ensembl_peptide_id','uniprotswissprot','percentage_gene_gc_content')
  #ens=biomaRt::useMart(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
  filters = list('biotype'='protein_coding','transcript_biotype'='protein_coding')
  if(with_uniprot){ filters = c(filters,'with_uniprotswissprot'=T)  }

  hs_gc_gene=getBM(attributes=att_gene, mart=BIOMART,
                   filters=names(filters), values=as.list(filters),
                   uniqueRows = T, bmHeader = F) %>% as_tibble()
  return(hs_gc_gene)
}

get_hs_chr = function(BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human'),
                      as.df=T,remove_patches=T,with_uniprot=T){
  library(biomaRt)
  chromosomes = c(1:22,'MT','X','Y')
  att_gene = c('ensembl_gene_id','ensembl_peptide_id','uniprotswissprot','chromosome_name')
  #ens=biomaRt::useMart(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
  filters = c('biotype'='protein_coding','transcript_biotype'='protein_coding')
  if(with_uniprot){ filters = c(filters,'with_uniprotswissprot'=T)  }

  hs_chr=getBM(attributes=att_gene, mart=BIOMART,
               filters=names(filters), values=as.list(filters),
               uniqueRows = T, bmHeader = F) %>% as_tibble() %>%
         mutate(is_patched = !(chromosome_name %in% chromosomes))

  if(remove_patches){
    hs_chr = hs_chr %>% dplyr::filter(is_patched)
  }

  patched_chr = setdiff(unique(hs_chr$chromosome_name),chromosomes)

  if(as.df){
    hs_chr= hs_chr %>%
            mutate( chr_val=T,
                    chromosome_name = fct_other(chromosome_name, keep=chromosomes, other_level = "other")) %>%
      pivot_wider(id_cols=c('ensembl_gene_id','ensembl_peptide_id','uniprotswissprot'),
                  names_from=chromosome_name, names_prefix='chr_',
                  values_from = chr_val, values_fill = F )

  }
  return(hs_chr)
}

get_hs_transcript = function(BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human'),
                             verbose=T,longest_transcript=F,with_uniprot=T){
  library(biomaRt)
  att_gene = c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id')
  att_pos = c('chromosome_name','start_position','end_position')
  att_struct = c('cds_length','transcript_length','transcript_start','transcript_end','ensembl_exon_id','rank','exon_chrom_start','exon_chrom_end','is_constitutive')
  att_uni = c('uniprotswissprot')
  att_type = c('gene_biotype','transcript_biotype')

  filters = c('biotype'='protein_coding','transcript_biotype'='protein_coding')

  if(with_uniprot){ filters = c(filters,'with_uniprotswissprot'=T)  }

  #hs_ens = useEnsembl('ensembl','hsapiens_gene_ensembl',mirror=ENS_MIRROR)
  # Get representative human proteome with UniProt/SwissProt identifiers
  hs_ensg = getBM(mart = BIOMART,
                  attributes = c(att_gene,att_pos,att_struct,att_uni,att_type),
                  filters=names(filters), values=as.list(filters),
                  uniqueRows = T, bmHeader = F) %>% as_tibble() %>%
    mutate( gene_length = end_position-start_position+1,
            exon_length = exon_chrom_end-exon_chrom_start+1 ) %>%
    group_by(ensembl_gene_id,ensembl_transcript_id,ensembl_peptide_id,uniprotswissprot) %>%
    mutate(n_exons = n_distinct(ensembl_exon_id), n_exons_mini = sum(is_constitutive),
           has_introns = n_exons_mini>1,
           tot_exon_len = sum(exon_length) ) %>%
    group_by(ensembl_gene_id) %>%
    mutate(n_transcripts = n_distinct(ensembl_transcript_id),
           n_proteins = n_distinct(ensembl_peptide_id),
           n_uniprot = n_distinct(uniprotswissprot))

  if(longest_transcript){
    hs_ensg = hs_ensg %>%
      group_by(ensembl_gene_id,uniprotswissprot) %>%
      dplyr::filter(transcript_length ==  max(transcript_length) ) #%>%
    #dplyr::select(-c(start_position,end_position,transcript_start,transcript_end,rank,
    #               ensembl_exon_id,is_constitutive,exon_chrom_start,exon_chrom_end,exon_length)) %>%
    #distinct()
  }

  nr = nrow(hs_ensg)
  ng = n_distinct(hs_ensg$ensembl_gene_id)
  nt = n_distinct(hs_ensg$ensembl_transcript_id)
  np = n_distinct(hs_ensg$ensembl_peptide_id)
  nu = n_distinct(hs_ensg$uniprotswissprot)
  if(verbose)
    message(sprintf('rows = %7s | genes = %6s | transcripts = %6s | proteins = %6s | uniprot = %6s',nr,ng,nt,np,nu))
  #dplyr::select(-start_position,-end_position,-transcript_start,-transcript_end,-exon_chrom_start,-exon_chrom_end)
  return(hs_ensg)
}

get_ens_filter_ortho = function(BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human')){
  library(biomaRt)
 # hs_ens = useEnsembl(host = host, biomart = mart, dataset = dat, mirror=ENS_MIRROR)
  filter_ortho = searchFilters(BIOMART,'homolog') %>%
                 as_tibble %>%
                 mutate(sp=str_split_fixed(name,'_',n=3)[,2]) %>%
                 mutate(Org = str_replace(description,pattern = "Orthologous (.+) Genes", replacement = "\\1")) %>%
                 dplyr::select(-description) %>%
                 mutate(org = tolower(Org))
  return(filter_ortho)
}

query_ens_ortho <- function(sp_ortho,COUNTER=1,
                            BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human'),
                            species='hsapiens',
                            dataset_suffix="gene_ensembl") {
  tictoc::tic('query Ensembl orthologs')
  ortholog = sprintf('with_%s_homolog',sp_ortho)
  #ens=useEnsembl(biomart = mart, dataset = dataset, host = host, mirror=ENS_MIRROR)

  BM_att = listAttributes(BIOMART,what='name')

  t0 = proc.time()
  out <- tryCatch({
    att_gene = c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id')
    # att_species = c('ensembl_gene','associated_gene_name','ensembl_peptide','canonical_transcript_protein','subtype',
    #                 'perc_id','perc_id_r1','goc_score','wga_coverage','orthology_confidence')
    #att_ortho = intersect(sprintf("%s_homolog_%s",sp_ortho,att_species), )

    att_ortho = grep(x=listAttributes(BIOMART,page='homologs',what='name'), paste0(sp_ortho,"_") ,value=T)
    .info$log(sprintf("Trying to fetch orthologs '%s' VS. '%s'...\n",species,sp_ortho))

    if(length(att_ortho)==0){
      .error$log(sprintf('No attributes for current species: %s!\n',sp_ortho))
      return(NULL)
    }

    Q=getBM(mart=BIOMART,
            attributes=c(att_gene,att_ortho),
            filters="", values="",
            uniqueRows = T, bmHeader = F)

    ortho_ensp = grep("ensembl_peptide",att_ortho,value=T)
    Q$no_ortholog = is.na(Q[[ortho_ensp]]) | (Q[[ortho_ensp]] == "")

    return(Q %>% as_tibble())
  },
  error=function(cond) {
    .error$log(sprintf('Failed to retrieve orthologs for current species: %s!\n',species))
    return(NULL)
  },
  finally={
    #.error$log(sprintf('done [%s]\n',COUNTER))
    elapsed=(proc.time()-t0)['elapsed']
    if(elapsed>5){ tictoc::toc() }
  }
  )
  return(out)
}

query_ens_txlen <- function(Fi,Va,ORG,COUNTER=1,verbose=T, debug=F,
                            BIOMART=get_ensembl_dataset('ENSEMBL_MART_ENSEMBL','human')) {
  if(verbose){ tictoc::tic('query Ensembl transcript length') }
  t0 = proc.time()
  out <- tryCatch({
    if(verbose){
      cat(sprintf("Trying to fetch structure of genes from '%s' [%s]...",BIOMART@dataset,ORG))
    }

    if(debug){ .dbg$log(c("Dataset to be used for query:",BIOMART@dataset)) }

    BM_att = listAttributes(BIOMART,what='name')

    att_gene = c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id')
    att_struct = c('cds_length','transcript_length')
    att_valid = intersect(c(att_gene,att_struct),BM_att)
    no_len = any(!att_struct %in% att_valid)

    if(debug){ .dbg$log(c("Valid attributes for query:",att_valid)) }

    if(no_len){
      warning('gene/transcript length unavailable!')
      pos_name = c("%s_position","transcript_%s","exon_chrom_%s")
      att_pos = c('ensembl_exon_id', 'is_constitutive',sprintf(pos_name,'start'), sprintf(pos_name,'end')) %>% sort
      att_valid = intersect(c(att_gene,att_pos),BM_att)
      if(any(!att_pos %in% att_valid)){
        warning('No information about gene length available')
        return(NULL)
      }
    }

    if(!missing(Fi) && !missing(Va)){
      Q=getBM(attributes=att_valid, mart=BIOMART, filters = Fi,  values = Va, uniqueRows = T, bmHeader = F) %>%
        group_by(ensembl_gene_id) #%>% dplyr::filter(transcript_length == max(transcript_length))
      return(Q %>% as_tibble())
    }else{
      Q=getBM(attributes=att_valid, mart=BIOMART, uniqueRows = T, filters="", values="", bmHeader = F) %>%
        group_by(ensembl_gene_id) #%>% dplyr::filter(transcript_length == max(transcript_length))

      if(no_len){
        if(debug){ .dbg$log("Computing gene/cds/transcript length...") }
        Qtx = Q %>%
          filter(ensembl_peptide_id != "" & !is.na(ensembl_peptide_id) & is_constitutive==1) %>%
          mutate(gene_length = end_position-start_position+1,
                 transcript_length = transcript_end-transcript_start+1,
                 exon_length = exon_chrom_end-exon_chrom_start+1 ) %>%
          dplyr::select(-contains(c('start','end'))) %>% ungroup() %>% distinct() %>%
          group_by(ensembl_gene_id,ensembl_transcript_id,ensembl_peptide_id) %>%
          mutate(cds_length = sum_(exon_length)) %>%
          dplyr::select(-contains('exon')) %>% ungroup() %>% distinct() %>%
          dplyr::select(-is_constitutive,-gene_length)
        return(Qtx %>% as_tibble())
      }
      return(Q %>% as_tibble())
    }
  },
  error=function(cond) {
    cat(sprintf('Failed to retrieve transcript length for current species: %s!\n',ORG))
    message(cond)
    return(NULL)
  },
  finally={
    if(verbose){
      cat(sprintf('done [%s]\n',COUNTER))
      elapsed=(proc.time()-t0)['elapsed']
      if(elapsed>5){ tictoc::toc() }
    }
  }
  )
  return(out)
}

##### NCBI Taxonomy #####
find_ncbi_lineage = function(){
  url_ncbi_tax = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
  file_ncbi_tax = basename(ncbi_tax)
  if( !file.exists(file_ncbi_tax) ){ download.file(ncbi_tax,file_ncbi_tax) }

  rk_lineage = 'rankedlineage.dmp'
  lineage = 'taxidlineage.dmp'

  untar(file_ncbi_tax,rk_lineage)
  hutils::replace_pattern_in(file_contents = "\t\\|",replace='', file_pattern=rk_lineage)
  ncbi_lineage_rk <- readr::read_delim(rk_lineage, delim='\t',
                                       col_names = c('tax_id','tax_name','species','genus','family','order','class','phylum','kingdom','superkingdom'))
  unlink(rk_lineage)

  # untar(file_ncbi_tax,lineage)
  # hutils::replace_pattern_in(file_contents = "\t\\|",replace='', file_pattern=lineage)
  # ncbi_lineage <- readr::read_delim(lineage, delim='\t',col_names = c('tax_id','node_ids'))
  # unlink(lineage)
  unlink(file_ncbi_tax)
  return(ncbi_lineage_rk)
}

find_ncbi_taxid = function(spnames,dbfile='data/ncbi/accessionTaxa.sql',verbose=F){
  library(taxonomizr)
  if( file.exists(dbfile) ){
    taxaId<-getId(taxa = spnames ,sqlFile = dbfile)
    if(verbose){ print(taxaId) }
    return(taxaId)
  }else{
    stop("NCBI database not found...Prepare the NCBI database with:\n prepareDatabase(dbfile)")
    #prepareDatabase(dbfile)
    return(NA)
  }
}

load.hgnc = function(with_protein=T, all_fields=F){

  hgnc = read_delim("http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt",delim = '\t') %>%
          dplyr::rename(uni=uniprot_ids,ensg=ensembl_gene_id)

  minified = c('uni','ensg','symbol','name','location','gene_group','locus_group','locus_type')

  if(with_protein){
    hgnc = hgnc %>% filter(locus_group == 'protein-coding gene' & locus_type=='gene with protein product') %>%
           dplyr::select(-locus_group,-locus_type)
  }

  if(!all_fields){ hgnc = dplyr::select(hgnc,any_of(minified)) }
  return(hgnc)
}
# Shortcut for loading datasets ------------------------------------------------

load.abundance = function(){
  # UNIFIED ABUNDANCE OF S.CEREVISIAE PROTEOME
  abundance = load.ho2018.data() %>%
    dplyr::select(orf, MPC = mean.mpc, MDPC = median.mpc, gfp=GFP.avg, ms=MS.avg) %>% # select the average expression
    group_by(orf) %>% dplyr::summarise(across(where(is.numeric),log10)) %>% # apply log10 to expression
    rowwise() %>% mutate( rowNA = all(is.na(c_across(where(is.numeric)))) ) %>% # check how many expression values
    dplyr::filter(rowNA || !is.na(MPC) ) %>% dplyr::select(-rowNA) # remove no values

  return(abundance)
}

load.annotation = function(only_ids=F){
  # Preloaded uniprot data can be generated in 5min with:
  #   uni = load.uniprot.features(tax="559292",refdb="UNIPROTKB")
  #   sgd = load.sgd.features()

  uni_feat = preload(here::here('data','uniprot-features.rds'),load.uniprot.features(tax="559292",refdb="UNIPROTKB")) %>%
    dplyr::select(-c(REVIEWED,COMMENTS,SUBLOC))
  sgd_desc =  preload(here::here('data','uniprot-sgd-annotation.rds'),load.sgd.features())
  biofunc = load.vanleeuwen2016.data(single_orf=T)
  enog_annot = eggnog_annotations_species(4891,4932)

  annotation = full_join(sgd_desc,uni_feat,by=c("SGD","UNIPROT"='UNIPROTKB')) %>%
    full_join(biofunc,by='ORF')  %>%
    full_join(enog_annot, by =c('ORF'='string')) %>%
    filter(!is.na(ORF) | is.na(UNIPROT)) %>%
    mutate(GENENAME = ifelse(is.na(GENENAME),ORF,GENENAME)) %>%
    mutate(letter = fct_drop(letter)) %>%
    relocate(SGD,GENENAME,ORF,UNIPROT,PNAME,
             L,FAMILIES,FUNCTION,ROLE,BIOPROCESS_all,enog_annot,
             LOC,COMPLEX,ORTHO,OTHER,KEYWORDS,
             EXISTENCE,SCORE)


  if(only_ids){
    identifiers = annotation %>% dplyr::select(UNIPROT,ORF,GENENAME,SGD,OG) %>%
      dplyr::filter(!duplicated(ORF) & !duplicated(UNIPROT) & !duplicated(SGD) & !duplicated(GENENAME))
    return(identifiers)
  }

  return(annotation)
}

load.network = function(net=c('string','intact')){
  ref_net = match.arg(net,choices = net, several.ok = F)
  if(ref_net=='string'){
    network = load.string(tax="4932",phy=F, ful=T, min.score = 700) %>%
      mutate(ORF1 = str_extract(protein1,SGD.nomenclature()),
             ORF2 = str_extract(protein2,SGD.nomenclature())
      ) %>% relocate(ORF1,ORF2) %>% dplyr::select(-c(protein1,protein2))
  }else if(ref_net=='intact'){
    network = load.intact.yeast(min.intact.score = 0.4) %>%
      dplyr::rename(ORF1=protA,ORF2=protB)
    #network = load.intact()
  }
  return(network)
}

load.clade = function(clade1='schizo',clade2='sacch.wgd'){
  # BRANCH LENGTH IN FUNGI CLADES
  # Normalized branch length (Kc):
  #  Kc = SUM(branch lengths in clade subtree Tc) / SUM(branch lengths in species tree Ts)
  enog.cols = c('Tax','NOG','nog.1to1','RAXML','FunCat','STRING','orf','sp1','sp2','uni1','uni2','ppm1','ppm2')
  clade.cols = c('ppm1.log10','ppm2.log10','clade1','clade2','XX','YY')
  clade = get_clade_data(g1=clade1,g2=clade2,rate = 'ratio') %>%
    dplyr::select(all_of(enog.cols), all_of(clade.cols), starts_with(clade1), starts_with(clade2)) %>%
    dplyr::rename(Kc1.log10=XX,Kc2.log10=YY)
  return(clade)
}

load.fungi.evo = function(){
  ### FUNGI CONTAIN ALL DATA ABOUT THE FUNGI LINEAGE EVOLUTIONARY RATE
  fungi = load.dubreuil2021.data(1) %>%
    dplyr::select('ORF','UNIPROT','PPM',c(starts_with('EVO.'))) %>%
    filter(!(is.dup(UNIPROT) & is.na(EVO.FULL) & is.na(PPM))) %>%
    #dplyr::mutate(across(starts_with("EVO."),function(x){ x/mean_(x) },.names="norm.{.col}")) %>% # scale and center all evo rate
    dplyr::mutate(across(starts_with("log10.EVO."),log10,.names = "log10.{.col}")) # # apply log10 to rate4site
  return(fungi)
}
