source("src/utils.r",local = T)
source("src/function_sequence.r",local = T)
source('src/function_alignment.r',local = T)
# Remote published proteome data -----------------------------------------------
load.dubreuil2019.data = function(d){
  # Load stickiness and disorder data
  message("REF: Dubreuil, Matalon and Levy ")
  message("Protein Abundance Biases the Amino Acid Composition of Disordered Regions to Minimize Non-functional Interactions")
  dubreuil=data.frame(stringsAsFactors = F,
                      name   = c("1.1 yeast.prot", "1.2 yeast.res", "1.3 yeast.res (full proteome)", "2.1 human.prot", "2.2 human.res") ,
                      num    = c('16924922', '16920611', '23164307', '16924925', '16920581'),
                      base_url =  rep('https://ndownloader.figshare.com/files/',5),
                      format = c('XLSX','TSV','TSV.GZ','XLSX','TSV')
  )
  if( missing(d) || !(d %in% seq_along(dubreuil$name)) ){
    d = menu(sprintf("%s (%s)",dubreuil$name,dubreuil$format), graphics = FALSE, title = "Which dataset do you want to use?")
  }
  choice = dubreuil[d,]
  with(choice,cat(sprintf("Your choice was:\n [%s] %s\n-----> from %s%s (formatted as %s)\n",d,name,base_url,num,format)))
  data.url = sprintf('%s/%s',dubreuil$base_url[d], dubreuil$num[d])

  if( dubreuil$format[d] == 'XLSX' ){
    require(openxlsx)
    res = openxlsx::read.xlsx(data.url,colNames=T,sheet=1, startRow=1)
  } else if( dubreuil$format[d] == 'TSV'){
    res = read.delim(file = data.url, header=T, sep='\t', stringsAsFactors = F)
  } else if( dubreuil$format[d] == 'TSV.GZ' ){
    res = read.delim(file=open.url(data.url), header=T, sep='\t',stringsAsFactors = F)
  }
  return(res)
}

load.leunberger2017.data = function(species='S. cerevisiae',rawdata=F){
  # Load protein stability data
  require(openxlsx)
  require(tidyverse)
  message("REF: P. Leuenberger et al., 2017, Science")
  message("Cell-wide analysis of protein thermal unfolding reveals determinants of thermostability")
  match.arg(species, choices = c('S. cerevisiae','E. coli', 'Human HeLa Cells','T. thermophilus'), several.ok = F)
  S3.url = "https://science.sciencemag.org/highwire/filestream/690833/field_highwire_adjunct_files/2/aai7825_Leuenberger_Table-S3.xlsx"
  #download.file(S3.url, destfile = "aai7825_Leuenberger_Table-S3.xlsx" )
  LIP_MS = read.xlsx(xlsxFile = S3.url,
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

load.jackson2018.data =function(){
  # Load 1011 yeast strains data
  # Sheet 1 = Strains details
  require(tidyverse)
  require(hablar)
  get.longest = function(S, s='\\.'){
    require(stringr)
    L = str_split(string = S, pattern = s)
    long=sapply(L,function(x){ nc=nchar(x); which.max(nc)})
    sapply(1:length(L),function(i){ L[[i]][long[i]] })
  }

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

load.vanleeuwen2020.data = function(){
  # Load gene dispensability inferred from bypass suppression of essential genes
  message("REF: Van Leeuwen et al., 2020, Molecular Systems Biology")
  message("Systematic analysis of bypass suppression of essential genes")
  #"https://www.embopress.org/doi/full/10.15252/msb.20209828"
  EV13.url = "https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20209828&file=msb209828-sup-0014-DatasetEV13.xlsx"
  #download.file(EV13.url, destfile = "msb209828-sup-0014-datasetev13.xlsx" )
  dispensable = read.xlsx(xlsxFile = EV13.url,
                          sheet = 2, detectDates = F,
                          skipEmptyRows = T, skipEmptyCols = T
  )
  colnames(dispensable) = c('ORF','gene','KO_exp','disp_score','literature','disp')
  return(dispensable)
}

load.villen2017.data = function(){
  # Load protein turnover (half-lives)
  message("REF: M. Martin-Perez and J. Villén, 2017, Cell Systems")
  message("Determinants and Regulation of Protein Turnover in Yeast")
  #"https://doi.org/10.1016/j.cels.2017.08.008"
  S1.url = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471217303411-mmc2.xlsx"
  #download.file(S1.url, destfile = "1-s2.0-S2405471217303411-mmc2.xlsx")
  turnover = read.xlsx(xlsxFile = S1.url,
                       sheet = 1, detectDates = F,
                       skipEmptyRows = T, skipEmptyCols = T, startRow = 4)

  colnames(turnover) = c('ORF','UNIPROT','GNAME','PNAME',
                         'half-life.avg','half-life.sd','half-life.cv',
                         'R1_hl.avg','R1_method','R1_reduced.x2','R1_adj.r2','R1_cv',
                         'R2_hl.avg','R2_method','R2_reduced.x2','R2_adj.r2','R2_cv')
  return(turnover)
}

load.belle2006.data = function(){
  # Load protein half-lives
  message("REF: A. Belle, A. Tanay, L. Bitincka, R. Shamir, E.K. O’Shea, 2006, PNAS")
  message("Quantification of protein half-lives in the budding yeast proteome")
  #https://doi.org./10.1073/pnas.0605420103
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

load.geisberg2014.data = function(nodesc=T){
  # Load mRNA half-lives
  message("REF: J.V. Geisberg et al., 2014, Cell")
  message("Global Analysis of mRNA Isoform Half-Lives Reveals Stabilizing and Destabilizing Elements in Yeast")
  #https://doi.org/10.1016/j.cell.2013.12.026
  S1.url = "https://ars.els-cdn.com/content/image/1-s2.0-S009286741301595X-mmc1.xlsx"
  #download.file(S1.url, destfile = "mmc1.xlsx")
  mrna.half = read.xlsx(xlsxFile = S1.url,
                        sheet = 1, detectDates = F,
                        skipEmptyRows = T, skipEmptyCols = T, startRow = 1)

  colnames(mrna.half)=c('chr','ORF','GNAME','type','mrna.half-life.mn','DESC')
  if(nodesc){ mrna.half$DESC = NULL }
  return(mrna.half)
}

load.ho2018.data = function(noauto = T,
                                 nogfp  = F,
                                 noms   = F,
                                 notap  = F) {
  require(openxlsx)
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
  S3 =  read.xlsx(
    xlsxFile = "https://ars.els-cdn.com/content/image/1-s2.0-S240547121730546X-mmc4.xlsx",
    sheet = 1, rows = 3:5751,
    colNames = T, skipEmptyRows = T, skipEmptyCols = T, na.strings = c("")
  )
  colnames(S3) = cols

  # Table S4. Protein Abundance in Molecules per Cell, after GFP Autofluorescence Filtering,
  S4 =  read.xlsx(
    xlsxFile = "https://ars.els-cdn.com/content/image/1-s2.0-S240547121730546X-mmc5.xlsx",
    sheet = 1, rows = 3:5861,
    colNames = T, skipEmptyRows = T, skipEmptyCols = T, na.strings = c("")
  )
  colnames(S4) = cols

  # get the correct columns matching each experiment
  strfind = function(strings, patterns){
    sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = T) })
  }

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
  geomean = function(x) {  exp(mean(log(x[x != 0 & !is.na(x)]))) }
  geosd = function(x) {  exp(sd(log(x[x != 0 & !is.na(x)]))) }
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

get.ygob.pair = function(ygob = load.ygob.ohnologs() ){
  return(ygob[, c('orf', 'dup.orf')])
}

get.sc.ohno = function(myseq) {
  ygob = load.byrne2005.data()
  bogy = ygob %>%
    dplyr::rename( orf = dup.orf, gname = dup.gname, dup.orf = orf, dup.gname = gname ) %>%
    mutate(ref = 2) %>%
    dplyr::select(WGD_anc, orf, gname, ref, dup.orf, dup.gname, pid, rlen)
  ohno = ygob %>% bind_rows(bogy)

  if(missing(myseq)){ myseq = load.sgd.proteome() }
  ohno.seq = get.pair.prot(prot=myseq, pair = get.ygob.pair(ohno))
  ohno.ali = align.pair.prot(p1=ohno.seq$s1, p2=ohno.seq$s2, mat='BLOSUM62')
  aafreq = alphabetFrequency(alignedSubject(ohno.ali))
  gaps = rowSums(aafreq[,c("-","+")])

  pair = get.ygob.pair(ohno) %>%
    mutate( L1 = width(ohno.seq$s1),
            L2 = width(ohno.seq$s2),
            RLEN = round(pmin(L1,L2) / pmax(L1,L2), 2),
            PID1 = round(pid( ohno.ali, "PID1" ),1),
            PID2 = round(pid( ohno.ali, "PID2" ),1),
            PID3 = round(pid( ohno.ali, "PID3" ),1),
            PID4 = round(pid( ohno.ali, "PID4" ),1),
            SCORE.B100 = score( ohno.ali ),
            S = nmatch( ohno.ali ),
            N = nmismatch( ohno.ali ),
            G = gaps,
            SGDID = AnnotationDbi::select(org.Sc.sgd.db,keys = ohno$orf,columns ="SGD",keytype = 'ORF')[,2],
            SGDID.dup = AnnotationDbi::select(org.Sc.sgd.db,keys = ohno$dup.orf,columns ="SGD",keytype = 'ORF')[,2]
    ) %>% right_join(ohno)

  return(pair)
}


