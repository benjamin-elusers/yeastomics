source("src/utils.r",local = T)
source("src/function_sequence.r",local = T)
source('src/function_alignment.r',local = T)
source('src/function_phylogenetic.r',local = T)
library(openxlsx)
open.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con)
  closeAllConnections()
  return(textConnection(txt))
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
  colnames(halflives)=c('ORF','GENE','half-life.mn','half-life.mn.corrected')
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
  costanzo = readxl::read_xls(path=, sheet = 1)
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
  message("REF: Michael D. Barton et al., PLOS One, 2010")
  message("Evolutionary Systems Biology of Amino Acid Biosynthetic Cost in Yeast")
  #https://doi.org/10.1371/journal.pone.0011935
  S1.url = "https://doi.org/10.1371/journal.pone.0011935.s001" # AA biosynthetic cost
  S4.url = "https://doi.org/10.1371/journal.pone.0011935.s004" # Protein biosynthetic cost
  BY = match.arg(by,choices=c('aa','prot'),several.ok = F)
  # Trasncript levels are obtained from 3 dilutions and 4 sources of nutrients published in:
  # JI Castrillo et al. 2007 -> "Growth control of the eukaryote cell: a systems biology study in yeast" (Journal of Biology)

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
                                names_glue = "{.value}.{cost_type}",
                                values_from=cost) %>%
                    rename_with(.cols=contains("-"),.fn=str_replace_all,pattern="-",replacement="_")

  }
  return(biosynth_cost)
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

  colnames(mrna.half)=c('chr','ORF','GNAME','type','mrna.half-life.mn','DESC')
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
  studies = c("Ingolia",paste0("Brar-",conditions), )
  url.datasets = gsub(" ","%20", x=sprintf("%s/%s %s.txt",url.MTDR,org,studies))
  names(url.datasets) = conditions
  exp.names = c('ingolia.exponential','brar.exponential_A14201','brar.exponential_gb15',
                'brar.DNA_replication', 'brar.recombination', 'brar.anaphase', 'brar.metaphase_I', 'brar.metaphase_II', 'brar.premeiotic_entry',
                'brar.spore_packing', 'brar.spores')
  MTDR = map(url.datasets,read.delim,header=F,col.names=c("orf","decoding_rate")) %>%
                  purrr::reduce(.x=., .f = left_join, by = "orf") %>%
                  rename_with(~exp.names,.cols=starts_with('decoding_rate'))
  return(MTDR)
}

load.lee2014.data = function(){
  # Load chemogenomic fitness signatures
  message("REF: A.Y. Lee,  R.P. St.Onge et al., 2014, Science")
  message("Mapping the Cellular Response to Small Molecules Using Chemogenomic Fitness Signatures")
  # https://doi.org/10.1126/science.1250217
  # Related to this first paper: The chemical genomic portrait of yeast: Uncovering a phenotype for all genes
  # E. Hillenmeyer, et al. Science 320,362–365 (2008).doi:10.1126/science.1150021

  library(openxlsx)
  # Compound library
  compound.library =  openxlsx::read.xlsx(
    xlsxFile = "http://chemogenomics.pharmacy.ubc.ca/hiphop/files/supplemental/leesupptableS1.xlsx",
    sheet = 2, colNames = T, skipEmptyRows = T, skipEmptyCols = T, na.strings = c(""))

  # Fitness Defect Score Matrix, homozygous strains (right-click to download, tab-delimited text, 280 Mb)
  # [Rows=Yeast Deletion Strain, Identified by Systematic ORF name; Columns=Fitness Screens, Identified by Screen ID (SGTC_N)]
  # Fitness defect on homozygous deletion strain for 3000 small molecules screen
  chemofit = read.delim(file="http://chemogenomics.pharmacy.ubc.ca/hiphop/files/supplemental/fitness_defect_matrix_hom.txt",
                        header=T, sep="\t", stringsAsFactors = F)
  return(chemofit)
}

load.vanleeuwen2016.data = function(){
  # load gene classification of biological functions (based on costanzo 2010)
  message("REF: J. Van Leeuwen et al., 2016, Science")
  message("Exploring genetic suppression interactions on a global scale")
  # doi: https://doi.org/10.1126/science.aag0839

  # Sheet 1 is the functional classes
  S7="https://science.sciencemag.org/highwire/filestream/686300/field_highwire_adjunct_files/6/aag0839TableS7.xlsx"
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
                         'half-life.avg','half-life.sd','half-life.cv',
                         'R1_hl.avg','R1_method','R1_reduced.x2','R1_adj.r2','R1_cv',
                         'R2_hl.avg','R2_method','R2_reduced.x2','R2_adj.r2','R2_cv')
  return(turnover)
}

load.leunberger2017.data = function(species='S. cerevisiae',rawdata=F){
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
                               skipEmptyRows = T, skipEmptyCols = T)

  peptides = LIP_MS %>%
    add_count(Protein_ID,name='npep') %>%
    distinct() %>%
    mutate( nres = round(0.01*Protein.Coverage*Length))

  if(!rawdata){
    peptides = peptides %>%
      dplyr::select(Protein_ID,Tm.Protein,
                    Protinfo,Protein.Coverage,Length,
                    Measured.Domains, Theoretical.Number.of.Domain,
                    Essential,npep )
  }
  # cols = c("Peptide.ID","Aggregation","Position","Tm.Peptide", "Tm.Protein",
  #          "Length","Sequence","Protein_ID","Pepinfo","Protinfo",
  #          "Secondary.Structure","Domain.logic","Domain.Name",
  #          "Theoretical.Number.of.Domain", "Measured.Domains","is.disordered",
  #          "Protein.Coverage","Protein.Abundance","Essential","T.90..Unfolded")
  #
  # "ID"           "ErrorMessage" "Position"     "Length"       "T0Abundance"  "lL"           "SSE"          "RMSE"         "R2"
  # [10] "dSTm"         "h"            "mf"           "mu"           "sf"           "su"           "tm"           "sig"          "h_cil"
  # [19] "h_ciu"        "mf_cil"       "mf_ciu"       "mu_cil"       "mu_ciu"       "sf_cil"       "sf_ciu"       "su_cil"       "su_ciu"
  # [28] "tm_cil"       "tm_ciu"       "sig_cil"      "sig_ciu"      "ProteinID"
  #
  return(peptides)
}

load.peter2018.data =function(){
  # Load 1011 yeast strains data
  # Sheet 1 = Strains details
  library(tidyverse)
  library(hablar)

  message("REF: J. Peter et al., 2018, Science")
  message("Genome evolution across 1,011 Saccharomyces cerevisiae isolates")
  SM.url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0030-5/MediaObjects/41586_2018_30_MOESM3_ESM.xls"
  todownload="41586_2018_30_MOESM3_ESM.xls"
  downloaded=download.file(SM.url, destfile = todownload)
  if(!downloaded){
    # collection = gdata::read.xls(xls = todownload, sheet = 1, method='tab',verbose = F,
    #                         blank.lines.skip=T, skip = 2, quote='',
    #                         stringsAsFactor=F) %>%
    #              mutate(across(everything(), ~ str_remove_all(., '\\\"'))) %>%
    #              mutate(across(everything(), ~ str_trim(.))) %>%
    #              mutate(X.Total.number.of.SNPs.=readr::parse_number(X.Total.number.of.SNPs.),
    #              X.Number.of.singletons.=readr::parse_number(X.Number.of.singletons.))
    collection = readxl::read_excel(path = todownload, sheet = 1, progress = T,skip = 2,col_names = T,n_max = 1011,guess_max = 900)
    #%>%
    clean_header = colnames(collection) %>%
      tolower %>% # CONVERT TO LOWERCASE
      gsub("[[:punct:]]","",x=.) %>% # REMOVE PUNCTUATIONS%>%
      gsub("\\s+",".",x=.) # SPACES TO UNDERSCORE

    #  str_replace_all('\\.+', rep= '\\.') %>% # REMOVE CONSECUTIVE DOTS
    #  stringr::str_sub(start = 3, end = -2) %>% # REMOVE starting 'X.' and ending '.'


    names(collection)=clean_header
    strains= tibble(type.convert(collection[1:1011,])) %>%   # 1011 first rows = strains details
      hablar::convert(hablar::chr(isolate.name),
                      hablar::chr(isolation),
                      hablar::chr(geographical.origins),
                      hablar::num(number.of.singletons),
                      hablar::chr(collection.provider)
      )
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
    file.remove(todownload)
    return(strains)
  }else{
    stop("Could not download the supplementary table!")
  }
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
    res = read.delim(file = data.url, header=T, sep='\t', stringsAsFactors = F)
  } else if( dubreuil$format[d] == 'TSV.GZ' ){
    res = read.delim(file=open.url(data.url), header=T, sep='\t',stringsAsFactors = F)
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
                      num    = c('26404466', '26404475', '26404478'),
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
    res = read.delim(file = data.url, header=T, sep='\t', stringsAsFactors = F)
  } else if( dubreuil$format[d] == 'TSV.GZ' ){
    res = read.delim(file=open.url(data.url), header=T, sep='\t',stringsAsFactors = F)
  }
  return(res)
}


# Resource databases with proteome data available ------------------------------

load.superfamily = function(tax='xs'){
  # Load domain assignments from superfamily SCOP level (SUPFAM.org)
  # xs = saccharomyces cerevisiae
  library(httr)
  library(readr)
  download.request =  httr::GET(sprintf("http://supfam.org/SUPERFAMILY/cgi-bin/save.cgi?var=%s;type=ass",tax))
  superfamily.txt = httr::content(download.request,as = 'text')
  superfamily.assignment = unlist(stringi::stri_split_lines(superfamily.txt,omit_empty = T))[-c(1:2)]
  supfam = readr::read_delim(superfamily.assignment,comment = "#",delim="\t", escape_double = F)
  return(supfam)
}

load.pfam = function(tax='559292'){
  # Load domains assignments based of HMM profiles (PFAM)
  # 559292 = S.cerevisiae
  library(readr)
  url.pfam = sprintf("http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes/%s.tsv.gz",tax)
  header = readLines(open.url(url.pfam))[3]
  # parse the header (3rd commented rows)
  columns = header %>%
            str_split(pattern = "> <") %>% unlist %>%
            str_replace_all("[ \\-]","_") %>% str_remove_all("[#<>]")
  pfam = readr::read_delim(file = url.pfam, skip=2,comment = "#",delim="\t", col_names = columns,escape_double = F,guess_max = 100)
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

  # each numeric column is an evidence score out of 1000 to asssess the existence of the link between protein
  # _transferred means information was inferred from a different organism (homology or orthologous group)
  STRING_net = readr::read_delim(STRING_url,delim = " ") %>% filter(combined_score >= min.score)
  return(STRING_net)
}


