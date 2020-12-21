source("utils.r")
# Local proteome data ----------------------------------------------------------

# Precomputed data #
#==================#
# raw.aligned.dataset = file.path(data.path, "202-aligned-residues-datasets.rds")
#if( !file.exists(dataset) ){ stop(sprintf("Cannot find the dataset at : %s",dataset)) }
# PROTEOME=readRDS('./data/PROTEIN-EVORATE.rds')
# ALIGNED.DATA = readRDS('./data/RESIDUE-EVORATE.rds')
#==================#

load.emmanuel.data = function(toolbox="/data/elevy/70_R_Data/bin/RToolBox_yeast_general.R"){
  source(toolbox,local = T)
  library(AnnotationDbi)
  library(org.Sc.sgd.db)
  library(GO.db)
  library(dplyr)
  SC = get.proteome.table() %>%
    mutate(
      has_len = !is.na(len),
      rel_diso1 = diso1/len, rel_diso2 = diso2/len, diso05 = diso05/len,
      has_diso05 = !is.na(diso05), has_diso1 = !is.na(diso1), has_diso2 = !is.na(diso2),
      has_tox = !is.na(over.tox),
      has_loc_ymd = !is.na(loc.ymd), has_loc_dtt = !is.na(loc.dtt),
      has_loc_h2o2 = !is.na(loc.h2o2), has_loc_starv = !is.na(loc.starv),
      has_viability = !is.na(viable)
    )
  return(SC)
}

load.aligned.data = function(data.path='../data/'){
  message("
          Merged dataset with residue-level informations from:
          (1) SGD S288C proteome sequence
          (2) Uniprot reference yeast proteome (aligned to SGD)
          (3) Dubreuil et al. (2019) for Abundance, Disorder (IUP+D2P2) & stickiness (+aaindex)
          (4) Rate4Site based on fungi lineage Wapinsky et al. (2007) of 14 yeast species
          (5) 3Dcomplex for yeast quaternary structures (based on X-Ray from PDB)"
  )
  return(readRDS(dataset))
}

# Evolutionary rate per regions on quaternary structures
get.evo3d =  function(path.data='../data/'){
  ALIGNED.RES = get.aligned.data(path.data) %>%
    mutate(iup=IUP20!=0, noiup = (IUP20==0),
           d2p2=d2p2_diso>=7, nod2p2= d2p2_diso<3,
           dom = domain!=0 & noiup, nodom = antidomain!=0) %>%
    dplyr::select( ends_with('.id'),
                   starts_with('l_'),starts_with('pos_'),
                   starts_with('aa_'), starts_with('gap_'),
                   starts_with('has_'),
                   starts_with('R4S_',ignore.case = F),
                   ends_with('.3d'),
                   c(iup,noiup,d2p2,nod2p2,dom,nodom,l_evo, l_pdb, gap_pdb),
                   -c(code.id, chain_name.id, org_ref.id,accession.id,uniprotAcc.id,
                      R4S_nmsa, R4S_orf, R4S_orf_sd, R4S_uni, R4S_uni_sd, R4S_std,
                      aa_sgd, gap_sgd,
                      aa_uniprot, gap_uniprot,
                      pos_dubreuil, aa_dubreuil, l_dubreuil, gap_dubreuil,
                      noseq.3d, npdb.3d,
                      homo.3d, best_BU.3d,hydro_kyte.3d, stickiness_ec.3d ) ) %>%
    filter(has_R4S)
  EVO3D= ALIGNED.RES %>% filter(has_R4S & has_pdb)
  return(EVO3D)
}

# Evolutionary rate per regions on quaternary structure with ASA cutoff
get.evo3d.byasa = function(asa,aligned.data=NULL){
  message(sprintf('Surface is ASA above %s%%.',asa))
  if(is.null(aligned.data)){  aligned.data = get.evo3d() }
  byasa = aligned.data %>%
    group_by(UNIPROTKB.id,SGD.id,ORF.id,code.chain.id) %>%
    mutate( CUTOFF_LEN=20, CUTOFF_ASA = asa,
            L_evo = sum_(has_R4S), L_pdb= sum_(has_pdb),
            ASA_under = asa_rel_in_BU.3d <= 25 & asa_rel_alone.3d <=25,  ASA_over = asa_rel_in_BU.3d > CUTOFF_ASA,
            L_buried=sum_(ASA_under), L_surface=sum_(ASA_over),
            L_iup = sum_(iup),  L_d2p2 = sum_(d2p2),
            L_noiup = sum_(noiup),  L_nod2p2 = sum_(nod2p2),
            L_dom = sum_(dom), L_nodom = sum_(nodom),
            R_full = ifelse(L_evo>=CUTOFF_LEN, mean_(R4S_norm), NA),
            R_pdb = ifelse(L_pdb>=CUTOFF_LEN, mean_(R4S_norm[!gap_pdb]), NA),
            R_surface   = ifelse(L_surface>=CUTOFF_LEN,   mean_(R4S_norm[ASA_over]),   NA),
            R_buried    = ifelse(L_buried>=CUTOFF_LEN,    mean_(R4S_norm[ASA_under]),    NA),
            R_iup   = ifelse(L_iup>=CUTOFF_LEN,   mean_(R4S_norm[iup]),   NA),
            R_d2p2   = ifelse(L_d2p2>=CUTOFF_LEN,   mean_(R4S_norm[d2p2]),   NA),
            R_noiup   = ifelse(L_noiup>=CUTOFF_LEN,   mean_(R4S_norm[noiup]),   NA),
            R_nod2p2   = ifelse(L_nod2p2>=CUTOFF_LEN,   mean_(R4S_norm[nod2p2]),   NA),
            R_dom   = ifelse(L_dom>=CUTOFF_LEN,   mean_(R4S_norm[dom]),   NA),
            R_nodom   = ifelse(L_nodom>=CUTOFF_LEN,   mean_(R4S_norm[nodom]),   NA)
    )    %>%
    dplyr::select( UNIPROTKB.id, SGD.id, code.chain.id, ORF.id,
                   has_iup, has_pdb, has_R4S,
                   starts_with(c('R_','L_','CUTOFF_'))
    ) %>%
    filter(!is.na(UNIPROTKB.id)) %>%
    ungroup() %>% distinct() %>%
    rename_with(function(x){ evo.cols[x]},.cols=names(evo.cols))
}

get.sc.ohno = function(myseq) {
  ygob = load.ygob.ohnologs()
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

load.leunberger2017.data = function(species='S. cerevisiae'){
  # Load protein stability data
  library(openxlsx)
  message("REF: P. Leuenberger et al., 2017, Science")
  message("Cell-wide analysis of protein thermal unfolding reveals determinants of thermostability")
  match.arg(species, choices = c('S. cerevisiae','E. coli', 'Human HeLa Cells','T. thermophilus'), several.ok = F)
  S3.url = "https://science.sciencemag.org/highwire/filestream/690833/field_highwire_adjunct_files/2/aai7825_Leuenberger_Table-S3.xlsx"
  #download.file(S3.url, destfile = "aai7825_Leuenberger_Table-S3.xlsx" )
  LIP_MS = read.xlsx(xlsxFile = S3.url,
                     sheet = species, detectDates = T,
                     skipEmptyRows = T, skipEmptyCols = T)
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
  return(LIP_MS)
}

load.jackson2018.data =function(){
  # Load 1011 yeast strains data
  # Sheet 1 = Strains details
  require(tidyverse)

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
      str_replace_all('\\.+', rep= '\\.') %>% # REMOVE CONSECUTIVE DOTS
      stringr::str_sub(start = 3, end = -2) %>% # REMOVE starting 'X.' and ending '.'
      tolower # CONVERT TO LOWERCASE

    names(collection)=clean_header

    strains= tibble(type.convert(collection[1:1011,])) %>%   # 1011 first rows = strains details
      hablar::convert(chr(isolate.name),
                      chr(isolation),
                      chr(geographical.origins),
                      num(number.of.singletons),
                      chr(collection.provider)
      )
    #PARSING REFERENCES
    refs = grep("^[0-9]+\\. ",collection$isolate.name,v=T) %>%
      str_split(pattern="\\. ",n=2,simplify = T) %>%
      as_tibble %>%
      rename(numref =V1, complete.ref=V2) %>%
      mutate(year  =  str_extract(complete.ref,pattern = "\\(([0-9]{4})\\)"),
             title =  get.longest(complete.ref),
             ref.2 =  str_remove(complete.ref,year),
             pages = str_extract(ref.2,"(\\d+)-(\\d+)"),
             ref.3 =  str_remove(ref.2,title),
             authors = str_extract(string = ref.3, pattern = '.+\\.\\.'),
             ref.4 = str_remove(ref.3,authors),
             journal.vol = str_extract(ref.3,pattern='\\.\\..+,'),
             journal = str_trim(str_extract(journal.vol,pattern='.+ ')),
             volume = str_remove(str_extract(journal.vol,pattern='[0-9]+,'),",")
             #ref.2=NULL, ref.3=NULL, ref.4=NULL,journal.vol=NULL
      )
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

load.merged.abundance = function(noauto = T,
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
