#source("src/utils.r",local = T)
#source("src/function_sequence.r",local = T)
#source('src/function_alignment.r',local = T)
#source('src/function_phylogenetic.r',local = T)

# TO REMOVE DEPENDENCIES, THOSE FUNCTIONS WERE COPIED FROM OTHER SCRIPTS
# Utils --------------------------------------------------------------------

library(openxlsx)
open.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con,skipNul=T)
  #closeAllConnections()
  close.connection(con)
  return(textConnection(txt))
}

read.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con,skipNul=T)
  #closeAllConnections()
  close.connection(con)
  return(txt)
}


strfind = function(strings, patterns){
  # Find several patterns in set of strings
  sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = T) })
}

get.longest = function(S, s='\\.'){
  # get the longest string in a list of splitted string
  library(stringr)
  L = str_split(string = S, pattern = s)
  long=sapply(L,function(x){ nc=nchar(x); which.max(nc)})
  sapply(1:length(L),function(i){ L[[i]][long[i]] })
}

# geometric mean and standard deviation
geomean = function(x) {  exp(mean(log(x[x != 0 & !is.na(x)]))) }
geosd = function(x) {  exp(sd(log(x[x != 0 & !is.na(x)]))) }

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
  Sdat1.url = "https://www.pnas.org/highwire/filestream/591690/field_highwire_adjunct_files/0/SuppDataSet.txt"
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

load.pu.2008.data = function(){
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
                                  col_names = c("ORF","GENE","FUNCTIONS"),na = c("unknown")) %>%
      separate_rows(FUNCTIONS,sep=";") %>%
      mutate(FUNCTION=str_trim(FUNCTIONS))
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
      pivot_wider(names_from=cost_type,
                  names_glue = "{.value}_{cost_type}",
                  values_from=cost) %>%
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
  message("Mean of the Typical Decoding Rates: A New Translation Efficiency Index Based on the Analysis of Ribosome Profiling Data ")
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

load.vanleeuwen2016.data = function(){
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
    mutate( FUNCTION = str_trim(Function),
            FN=factor(FUNCTION,levels=unique(FUNCTION),labels = invert(.class_function))
    ) %>% dplyr::select(-Function) %>%
    group_by(ORF) %>% mutate(FN_all=paste0(unique(FN),collapse = ""))

  return(vanleeuwen)
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
  S3.url = "https://science.sciencemag.org/highwire/filestream/690833/field_highwire_adjunct_files/2/aai7825_Leuenberger_Table-S3.xlsx"
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
                    essential,npep )
  }
  # cols = c("Peptide.ID","Aggregation","Position","Tm.Peptide",
  #         "Sequence","Pepinfo", "Secondary.Structure","Domain.logic","Domain.Name",
  #          "is.disordered", "T.90..Unfolded")

  return(peptides)
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
    res = readr::read_delim(file = read.url(data.url), col_names = T, delim = '\t',  progress = T,guess_max=50000)
  }
  return(res)
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
  if(!file.exists(uniprot)){
    # Check whether a uniprot is associated to an alphafold model
    url_alphafold="https://alphafold.ebi.ac.uk/files"
    ext=match.arg(extension,extension)
    AF_uniprot=sprintf("%s/AF-%s-F1-model_v1.%s",url_alphafold,uniprot,ext)
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

load.superfamily = function(tax='xs'){
  # Load domain assignments from superfamily SCOP level (SUPFAM.org)
  # xs = saccharomyces cerevisiae
  library(httr)
  library(readr)
  download.request =  httr::GET(sprintf("http://supfam.org/SUPERFAMILY/cgi-bin/save.cgi?var=%s;type=ass",tax))
  superfamily.txt = httr::content(download.request,as = 'text')
  superfamily.assignment = unlist(stringi::stri_split_lines(superfamily.txt,omit_empty = T))[-c(1:2)]
  supfam = readr::read_delim(superfamily.assignment,comment = "#",delim="\t", escape_double = F) %>% janitor::clean_names()
  return(supfam)
}

load.pfam = function(tax='559292'){
  # Load domains assignments based of HMM profiles (PFAM)
  # 559292 = S.cerevisiae
  library(readr)
  .release = read.url("http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz") %>% paste0(sep="\n")
  message("CURRENT RELEASE")
  message("---------------")
  message(.release)
  message("---------------")

  url.pfam = sprintf("http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes/%s.tsv.gz",tax)
  url.clan = sprintf("http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz")
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

load.string = function(tax="4932",phy=T,ful=T,min.score=700){
  # Load protein links (physical PPI)
  # 4932 = S.cerevisiae
  STRING_baseurl = "https://stringdb-static.org/download"
  .version = "v11.0"
  .protein = ifelse(phy,"protein.physical","protein")
  .links = ifelse(ful,"links.full","links")
  STRING_dataset = sprintf("%s.%s.%s",.protein,.links,.version)
  STRING_url = sprintf("%s/%s/%s.%s.%s",STRING_baseurl,STRING_dataset,tax,STRING_dataset,"txt.gz")

  message("STRING VERSION: ", .version)
  message("STRING DATASET: ",STRING_dataset)
  message(sprintf("FILTERED BY:\n physical links = %s \nAND\n full evidence =  %s \nAND\n min. score >= %s", phy, ful, min.score))
  # each numeric column is an evidence score out of 1000 to asssess the existence of the link between protein
  # _transferred means information was inferred from a different organism (homology or orthologous group)
  STRING_net = readr::read_delim(STRING_url,delim = " ") %>% filter(combined_score >= min.score)
  return(STRING_net)
}

load.pombe.orthologs = function() {
  # Load orthologs pombe-cerevisiae
  # Fused genes have parentheses to indicate which fused side they are
  bracket.before = "(?<=\\()"
  bracket.after = "(?=\\))"
  regex_fusion = paste0( bracket.before, "(FUSION-)?(N|C)?", bracket.after)

  url.orthologs = "ftp://ftp.pombase.org/pombe/orthologs/cerevisiae-orthologs.txt"
  sp.sc = readr::read_delim(url.orthologs, delim="\t", col_names=c('PombaseID','ORFS'), comment="#", trim_ws=T) %>%
    mutate(grp = sprintf("OG_%04d",row_number())) %>%
    separate_rows(ORFS,sep = "\\|") %>%
    mutate(orf= dplyr::na_if(str_extract(ORFS,"^[^\\(]+"),"NONE"),  # get the orf,
           fused_side=str_sub(start=-1,str_extract(ORFS,regex_fusion)) )%>%  # get the parentheses content
    filter(!is.na(orf)) %>%
    group_by(PombaseID) %>% mutate(sp1 = n() ==1 ) %>%
    group_by(orf) %>% mutate(sc1= n() == 1) %>%
    rowwise %>% mutate( is_1to1 = sp1 & sc1,
                        is_fusion = !is.na(fused_side)) %>%
    dplyr::select(-ORFS) %>% ungroup() %>% arrange(orf)

  return(sp.sc)
}

load.paxdb.orthologs = function(node,show.nodes=F) {
  # get table of orthologs proteins from paxdb
  paxdb_ortho="https://pax-db.org/downloads/latest/paxdb-orthologs.zip"
  library(rio)
  library(RCurl)
  library(hablar)
  # download the archive
  temp<-tempfile()
  download.file(paxdb_ortho,temp)
  # Find the list of files in the archive (taxonomic nodes)
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

  node = paxdb_nodes[numnode]
  paxdb_node = sprintf("paxdb-orthologs-4.1/%s.orthgroups.consistent.txt",node)
  # Read orthogroups at a specific taxonomic level
  message("Retrieving orthologs for node [",node,"]...")
  node_ortho <- readr::read_delim(unz(temp, paxdb_node),col_names = c("NOG", "orthologs"), delim="\t",progress=T) %>%
    separate_rows(orthologs,sep=" ") %>%
    extract(orthologs,into=c('taxid','protid'),regex='(^[0-9]+)\\.(.+)') %>%
    group_by(NOG,taxid) %>% mutate(is_1to1=n()==1)

  return(node_ortho)
}

load.paxdb = function(taxon=4932){

  paxdb_dataset ="https://pax-db.org/downloads/latest/datasets/"
  taxon_dir=file.path(paxdb_dataset,taxon,"/")
  taxon_data <- RCurl::getURL(taxon_dir,dirlistonly = TRUE) %>%
    stringr::str_extract_all('(?<=\\<a href\\=\\")(.+\\.txt)(?=\\">)') %>%
    unlist

  Ndata = length(taxon_data)
  message(Ndata," paxDB datasets for taxon [",taxon,"]")

  taxon_url = paste0(taxon_dir,taxon_data)
  get.paxdb_header = function(urldata){
    read.url(urldata) %>%
      str_subset(pattern="^#") %>%
      as_tibble %>%
      extract(col=value,into=c('info',"value"), regex="^#([^\\:]+)\\:(.+)$") %>%
      mutate(info=str_trim(info),value=str_trim(value)) %>%
      filter( !is.na(info) ) %>%
      pivot_wider(names_from='info',values_from=c(value)) %>%
      mutate(taxid=as.character(taxon))
  }

  infodata = map_dfr(taxon_url,get.paxdb_header) %>%
    mutate( w=parse_number(weight)*0.01,
            ndata = n_distinct(id,filename),
            is_integrated = integrated=='true' | ndata==1) %>%
    dplyr::select(taxid,organ,ndata,id,filename,is_integrated,score,w,cov=coverage,yr=publication_year)

  # if( Ndata == 1){
  # ppm = rio::import(taxon_url) %>%
  #   rename_with(~c("paxid",'string','ppm','count')) %>%
  #    mutate(dataset = basename(taxon_url)) %>%
  #     extract(string,into=c('taxid','protid'),regex='(^[0-9]+)\\.(.+)') %>%
  #      mutate()
  # }else{
  # If only one file, read twice the url and keep the distinct rows
    ppm = rio::import_list(file=c(taxon_url[1],taxon_url),
                           setclass="tibble",
                           rbind = TRUE,rbind_label = "dataurl",rbind_fill = T) %>%
      dplyr::rename(paxid="#internal_id", string="string_external_id", ppm="abundance") %>%
      mutate(dataset = basename(path = dataurl)) %>% dplyr::select(-dataurl) %>%
      extract(string,into=c('taxid','protid'),regex='(^[0-9]+)\\.(.+)') %>%
      distinct()
  #}#

  taxon_ppm = left_join(ppm,infodata, by=c('dataset'='filename','taxid')) %>%
    arrange(protid,ppm)
  return(taxon_ppm)
}

get.paxdb = function(tax=4932, abundance='integrated'){
  closeAllConnections() # Make sure to closse connections
  paxdb = load.paxdb(tax)
  # if nothing selected return the integrated values
  # (if there is a single dataset, it is considered as integrated)
  targets=match.arg(abundance,c('integrated','median','mean','weighted'), T)
  message(sprintf("---> Returning (%s) abundance values...",toString(targets)))
  RES = list()

  if( "integrated" %in% targets ){
    RES$INT = paxdb %>%
      dplyr::group_by(taxid,organ,protid) %>%
      mutate(ppm_n = ndata-sum(is_integrated)) %>%
      filter(is_integrated) %>%
      dplyr::select(taxid,organ,protid, ppm_int = ppm,ppm_n)
  }

  if( "median" %in% targets ){
    RES$MED = paxdb %>%
      dplyr::filter(!is_integrated) %>%
      group_by(taxid,organ,protid) %>%
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
      group_by(taxid,organ,protid) %>%
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
      group_by(taxid,organ,protid) %>%
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
  return(purrr::reduce(.x=RES,.f=left_join, by = c("taxid","organ","protid")) )
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
  closeAllConnections() # Make sure to close all connections before downloading data

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
      arrange(NOG,taxid,protid)

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

get.eggnogg.node=function(node=4890){
  url_eggnog = "http://eggnog5.embl.de/download/latest/per_tax_level/"
  #paste0(url_eggnog,node)
}
