source(here::here("src","__setup_yeastomics__.r"))
library(log)
.info  = infoLog()
.error =  errorLog()
.warn  = warningLog()
.succ  = successLog()

#library(extrafont)
#library(extrafontdb)
#library(remotes)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
#extrafont::font_import()
#loadfonts()
# library(hrbrthemes)
# library(showtext)
# showtext::showtext_auto()
# font_add_google("Roboto", "roboto")
# font_add_google("Roboto Mono", "roboto-mono")
# font_add_google("Roboto Condensed", "roboto-condensed")
# font_add_google("Open Sans", "open-sans")
# font_add_google("Montserrat", "montserrat")
# font_add_google("Poppins", "poppins")
# font_add_google("Oswald", "oswald")
# font_add_google("Varela Round", "varela")
# font_add_google("Pacifico", "pacifico")
# font_add_google("Bangers", "bangers")
# font_add_google("Abel", "abel")
# font_add_google("Lobster", "lobster")
# font_add_google("Quicksand", "quicksand")
# font_add_google("Indie Flower", "indie-flower")
# font_add_google("Dancing Script", "dancing-script")
# font_add_google("Titillium Web", "titillium-web")
# font_add_google("Permanent Marker", "permanent-marker")
# font_add_google("Amatic SC", "amatic-sc")
# font_add_google("Patrick Hand", "patrick-hand")
# font_add_google("Special Elite", "special-elite")
# font_add_google("Poiret One", "poiret-one")
# font_add_google("Luckiest Guy", "luckiest-guy")
# font_add_google("Bungee", "bungee")

#hrbrthemes::import_roboto_condensed()
#hrbrthemes::import_titillium_web()
AXIS_TITLE_SIZE=16
AXIS_TEXT_SIZE=12
TEXT_SIZE=4
#AXIS_TEXT_FONT = 'Roboto'
# th_txt_size = hrbrthemes::theme_ipsum_rc(
#   axis_title_just = 'm', axis = 'xy', axis_col = 'black',
#   grid = F, base_family = 'titillium-web',
#   axis_text_size = AXIS_TEXT_SIZE,
#   axis_title_size = AXIS_TITLE_SIZE)

th_txt_size = hrbrthemes::theme_ipsum(
  axis_title_just = 'm', axis = 'xy', axis_col = 'black',
  grid = F, base_family = 'Helvetica',
  axis_text_size = AXIS_TEXT_SIZE,
  axis_title_size = AXIS_TITLE_SIZE)



# 1. LOAD DATASETS -------------------------------------------------------------
load.abundance = function(){
  # UNIFIED ABUNDANCE OF S.CEREVISIAE PROTEOME
  abundance = load.ho2018.data() %>%
    dplyr::select(orf, MPC = mean.mpc, MDPC = median.mpc, gfp=GFP.avg, ms=MS.avg) %>% # select the average expression
    group_by(orf) %>% dplyr::summarise(across(where(is.numeric),log10)) %>% # apply log10 to expression
    rowwise() %>% mutate( rowNA = all(is.na(c_across(where(is.numeric)))) ) %>% # check how many expression values
    dplyr::filter(rowNA || !is.na(MPC) ) %>% dplyr::select(-rowNA) # remove no values

  return(abundance)
}

load.annotation = function(){
  # Preloaded uniprot data can be generated in 5min with:
  #   uni = load.uniprot.features(tax="559292",refdb="UNIPROTKB")
  #   sgd = load.sgd.features()

  uni_feat = read_rds(here('data','uniprot-features.rds')) %>%
      dplyr::select(-c(REVIEWED,COMMENTS,SUBLOC))
  sgd_desc = read_rds(here('data','uniprot-sgd-annotation.rds'))
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


  return(annotation)
}

load.network = function(net=c('string','intact'),taxon=4932){
  ref_net = match.arg(net,choices = net, several.ok = F)
  if(ref_net=='string'){
    network = load.string(tax=taxon,phy=F, ful=T, min.score = 700)
    if(taxon==4932){
      network = network %>% relocate(protein1,protein2) %>%
      mutate(protein1 = str_extract(protein1,SGD.nomenclature()),
             protein2 = str_extract(protein2,SGD.nomenclature())
      )
    }else if(taxon==9606){
      network = network %>%
        mutate(protein1 = str_extract(protein1,ENSEMBL.nomenclature()),
               protein2 = str_extract(protein2,ENSEMBL.nomenclature())
        )
    }
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

load.strains.evo = function(){
  ### STRAINS CONTAIN ALL DATA ABOUT THE FUNGI LINEAGE EVOLUTIONARY RATE AND THE S.CEREVISIAE POPULATION
  strains = readRDS(here("data","PROTEIN-EVO-FUNGI-SNP.rds")) %>%
    dplyr::select(c(starts_with(c('EVO.','SNP.')),'PPM','ORF','UNIPROT','IS_FUNGI','IS_STRAINS')) %>%
    filter(!(is.dup(UNIPROT) & is.na(EVO.FULL_R4S) & is.na(PPM))) %>%
    #dplyr::mutate(across(starts_with(c("EVO.","SNP.")),function(x){ x/mean_(x) },.names="norm.{.col}")) %>% # scale and center all EVO/SNP rate
    dplyr::mutate(across(starts_with(c("EVO.","SNP.")),log10,.names = "log10.{.col}")) %>% # apply log10 to SNP rate4site
    relocate(ORF,UNIPROT,PPM,IS_FUNGI,IS_STRAINS) %>%
  return(strains)
}

load.properties=function(tolong=F){
  prop=readRDS(get.last.file(here("output"),"proteome-properties")) %>%
    dplyr::select(-contains("chemsig_X")) %>% ungroup()
  if(tolong){
    long_prop = prop %>%
      pivot_longer(cols = starts_with('cat_'),
                   names_to = c('categories','source','property'),
                   names_pattern="cat_(.+)\\.(.+)\\.(.+)",
                   values_to = "has_prop") %>%
      mutate(col_prop = paste0(categories,'.',source,'.',property)) %>%
      dplyr::filter(has_prop)
    return(long_prop)
  }
  return(prop)
}

load.features=function(tolong=F){
  feat = readRDS(get.last.file(here("output"),"proteome-features")) %>%
          ungroup() %>%
          distinct()

  if(tolong){
    long_feat = feat %>%
    pivot_longer(cols = starts_with('cat_'),
                 names_to = c('categories','source','feature'),
                 names_pattern="cat_(.+)\\.(.+)\\.(.+)",
                 values_to = "value") %>%
    mutate(col_feat = paste0(categories,'.',source,'.',feature)) %>%
    dplyr::filter(!is.na(value)) %>% add_count(col_feat,name='size')
    return(long_feat)
  }
  return(feat)
}

normalize_features=function(feat,automated=F){

  # FEATURES NORMALIZATION
  feat_norm = feat

  # a. Normalize codons counts (0-100) => 64 COLUMNS ---------------------------
  A=col_codons = grep("cat_transcriptomics.sgd.[ATCG]{3}",colnames(feat),v=T)
  feat_norm[,col_codons] =  100 * feat[,col_codons] / (feat$cat_genomics.sgd.len+1)

  # b. Normalize amino acid frequencies (0-100) => 31 COLUMNS  -----------------
  # 04.01.22 changed to use sgd sequences instead of uniprot
  B=col_f_aa = grep("cat_biophysics.sgd.f_",colnames(feat),v=T)
  feat_norm[,col_f_aa] =  100 * feat[,col_f_aa]

  # c. Normalize all other fraction (0-100) => 4 COLUMNS  ----------------------
  C=col_frac = c("cat_genomics.sgd.pGC","cat_genomics.byrne2005.RLEN","cat_transcriptomics.coRdon.CU_fop",
               "cat_biophysics.d2p2.f")
  #"cat_biophysics.dubreuil2019.IUP20_f","cat_biophysics.dubreuil2019.IUP30_f","cat_biophysics.dubreuil2019.IUP40_f")
  feat_norm[,col_frac] =  100 * feat[,col_frac]

  # d. Normalize centrality measures variables => 44 COLUMNS ---------------------------------
  #install.packages('bestNormalize')
  library('bestNormalize')
  D=col_cent = feat %>% dplyr::select(contains('.cent_')) %>% colnames

  message("Find best transformation using 'bestNormalize' on network centrality measures...")
  icent=1
  for(cent in col_cent){
    cat(sprintf("%2s/%2s - [%30s]...\n",icent,length(col_cent),cent))
    CENT = as.numeric(feat %>% pull(cent))
    CENT_norm = predict(bestNormalize(CENT,quiet = T))
    feat_norm[,col_cent] = CENT_norm
    icent=icent+1
  }

  # cent.range =  feat %>% ungroup() %>%
  #               summarise( across(all_of(col_cent),.fns = function(x){ max_(x) / min_(x[x>0]) }) )
  #
  # log10_= function(x){ x[!is.na(x) & x>0] = log10(x[!is.na(x) & x>0]) ; return(x) }
  #
  # D=col_cent_log10 = names(cent.range)[cent.range>50]
  # feat_norm[,col_cent_log10] = apply(feat[,col_cent_log10],2, log10_)

  # e. Normalize halflife proteins => 2 COLUMNS ---------------------------------
  E=col_halflife_prot= c("cat_biophysics.belle2006.HL_prot","cat_biophysics.christiano2014.HL_prot","cat_transcriptomics.geisberg2014.HL_mrna", "cat_biophysics.villen2017.HL_prot")
  feat_norm[,col_halflife_prot] =  feat[,col_halflife_prot]/3600

  feat_norm[,"cat_biophysics.pepstats.mw"]  = feat[,"cat_biophysics.pepstats.mw"] / 1000
  # f. Normalize count variables => 14 COLUMNS ---------------------------------
  FF=col_doublings = c("cat_genomics.sgd.len","cat_transcriptomics.sgd.prot_size", "cat_transcriptomics.sgd.len_cds",
                    "cat_genomics.byrne2005.S","cat_genomics.byrne2005.N","cat_genomics.byrne2005.G",
                    "cat_biophysics.d2p2.nseg","cat_biophysics.d2p2.L","cat_biophysics.d2p2.Lsegmax",
                    "cat_genomics.peter2018.snp_count",
                    "cat_biophysics.dubreuil2019.IUP20_L","cat_biophysics.dubreuil2019.IUP30_L",
                    "cat_biophysics.dubreuil2019.IUP40_L")
  feat_norm[,col_doublings] =  log2(feat[,col_doublings]+1)

  n=length(feat)
  a=length(A); b=length(B); c=length(C); d=length(D); e=length(E); f=length(FF)

  z = setdiff(colnames(feat), c(A,B,C,D,E,FF))
  #sapply(feat[,z],max_)

  if(automated==T){
    message("Find best transformation using 'bestNormalize' for all other variables...")
    iother=1
    for(other in z){
      cat(sprintf("%2s/%2s - [%30s]...\n",iother,length(z),other))
      OTHER = as.numeric(feat %>% pull(other))
      OTHER_norm = predict(bestNormalize(OTHER,quiet = T))
      feat_norm[,other] = OTHER_norm
      iother=iother+1
    }
  }
  return(feat_norm)
}


# 2. FIX MISSING VALUES --------------------------------------------------------
check_missing_var = function(input){
  library(skimr)
  # Filter columns that contains NA
  input_na = input %>% dplyr::select(where(~ any(is.na(.x))))
  missing_var = input_na %>% skimr::skim(.) %>% filter(complete_rate<1)
  return(missing_var)
}

#### A. IN COLUMNS ####
get_binary_col = function(df,only.names=F){
  # Binary variables (only 2 outcomes)
  df_bin = df %>% dplyr::select(where(~is.binary(.x)))
  if(only.names){ return( colnames(df_bin) ) }
  return(df_bin)
}

remove_rare_vars = function(df,min_obs=2){

  # Find binary variables with rare observations (preferably 0's and singletons)
  binary_vars = get_binary_col(df)
  rare_vars = binary_vars[ colSums(binary_vars,na.rm = T) < min_obs ]
  n_rare = length(rare_vars)
  .succ$log(sprintf("Excluding %s/%s predictors with less than %s observations\n",n_rare,ncol(df),min_obs))

  # Find variables with near-0 variance
  num_vars = df %>% dplyr::select(where(is.numeric)) %>% apply(., 2,var_) %>% na.omit()
  constant_vars = num_vars[ num_vars < .Machine$double.eps ]
  n_constant = length(constant_vars)
  .succ$log(sprintf("Excluding %s/%s predictors with low varianace (<%.1e) \n",n_constant,ncol(df), .Machine$double.eps))

  .succ$log(sprintf("Excluded %s/%s useless (rare/low variance) predictors\n",n_constant+n_rare,ncol(df)))
  df_fixed = df %>% dplyr::select(-all_of(c(names(rare_vars),names(constant_vars))))
  return(df_fixed)
}


#### B. IN ROWS ####
#### 1. codons counts ####
get_codons_col = function(df,col_prefix='cat_transcriptomics.sgd.'){
  regex_codons=paste0(col_prefix,get.codons4tai())
  res = df %>% dplyr::select(matches(regex_codons,ignore.case = F))
  return(res)
}

retrieve_missing_codons = function(orf_missing){
  .info$log(sprintf("Retrieving missing codons counts for %s ORF...\n",n_distinct(orf_missing)))
  # Retrieve CDS and count codons for orf with missing values
  cod = paste0(get.codons4tai(),"$")
  cds = load.sgd.CDS()
  codon_counts = coRdon::codonCounts(coRdon::codonTable(cds[orf_missing]))
  codon_missing = tibble(ORF=orf_missing, as_tibble(codon_counts))
  return(codon_missing)
}

fix_missing_codons = function(df,col_prefix='cat_transcriptomics.sgd.'){
  # Replace orf with missing values with retrieved codons counts from CDS
  orf_missing = df %>% column_to_rownames('ORF') %>% get_codons_col(col_prefix) %>% find_na_rows() %>% rownames()
  if( length(orf_missing) == 0){
    .warn$log("No missing values found for any codons count columns!")
    return(df)
  }

  col_codons=get_codons_col(df,col_prefix) %>% colnames()

  codon_table = seqinr::SEQINR.UTIL$CODON.AA %>% as_tibble() %>%
                mutate(CODON=toupper(CODON), AA=str_to_title(AA),
                       L=str_replace(string = L, "\\*", "STOP"),
                       codon_aa = paste0(CODON,"_",AA,"_",L))
  .warn$log("Replace columns of codons counts with missing values...\n")
  df_na_codon = retrieve_missing_codons(orf_missing)
  codon_aa=codon_table$codon_aa[codon_table$CODON %in% colnames(df_na_codon)]
  df_na_codon = df_na_codon %>%  dplyr::rename_with(.cols=-ORF,~paste0(col_prefix,codon_aa))
  df_fixed = coalesce_join(x = df, y=df_na_codon, by = "ORF")
  return(df_fixed)
}

#### 2. protein length / Average Molecular Weight ####
fix_missing_peptide_stats = function(df,id='ORF', taxon=4932,
                                     col_len='cat_transcriptomics.sgd.prot_size',
                                     col_mw='cat_biophysics.pepstats.mw',
                                     col_mw_avg='cat_transcriptomics.pepstats.mean_MW',
                                     col_charge='cat_biophysics.pepstats.netcharge',
                                     col_pi='cat_biophysics.pepstats.pI'
                                     ){
  col_pep = c(col_len,col_mw,col_mw_avg,col_charge,col_pi)

  df$has_unique_id = !is.na(df[[id]]) & !duplicated(df[[id]])
  df_unique = df %>% drop_na(id) %>% rowwise %>% filter(has_unique_id) %>% column_to_rownames(id)
  orf_missing =  df_unique %>% dplyr::select(all_of(col_pep)) %>% find_na_rows() %>% rownames()
  if(length(orf_missing)==0){
    .warn$log("No missing values found for any of the peptide stats columns!")
    return(df)
  }
  .warn$log("Replace columns with missing values for protein length / average molecular weight...\n")
  if(taxon==4932){
    prot_missing = load.sgd.proteome()[orf_missing] %>% as.character
  }else if(taxon==9606){
    prot_missing = get.uniprot.proteome(taxid = 9606,DNA = F,fulldesc = F)[orf_missing] %>% as.character
  }

  library(Peptides)
  df_na_pep = tibble(id=orf_missing,
                    missing_len = Peptides::lengthpep(prot_missing),
                    missing_mw  = Peptides::mw(prot_missing),
                    missing_mw_avg = missing_mw / missing_len,
                    missing_charge = Peptides::charge(prot_missing),
                    missing_pi = Peptides::pI(prot_missing),
                    PEPSTATS.AA_costly = missing_mw_avg > 118,
                    PEPSTATS.AA_cheap = missing_mw_avg <= 105) %>%
    dplyr::rename(!!col_len:=missing_len,
                  !!col_mw:=missing_mw,
                  !!col_mw_avg:=missing_mw_avg,
                  !!col_charge:=missing_charge,
                  !!col_pi:=missing_pi
                  )
  if(taxon == 9606){
    df_na_pep = df_na_pep %>% dplyr::rename(uniprot=id)
  }
  df_fixed = coalesce_join(x = df, y=df_na_pep, by = c('uniprot'))
  return(df_fixed)
}

#### 3. network centrality ####
get_centrality_col = function(df,col_prefix="cat_interactions.string."){
  regex_centrality = paste0("^",col_prefix,"cent_")
  res = df %>% dplyr::select(matches(regex_centrality,ignore.case = F))
  return(res)
}

retrieve_missing_centrality = function(orf_missing,type='string',taxon=4932){
  interactions=load.network(type,taxon)
  centrality_low.rds = here::here('data',paste0(taxon,'-',type,'-low_stringency_centrality.rds'))
  centrality_low = preload(saved.file = centrality_low.rds,
                            loading.call = {
                              cent=interactions %>%
                              #filter(ORF1 %in% orf_missing | ORF2 %in% orf_missing) %>%
                              dplyr::select(protein1,protein2) %>%
                              network.centrality(fromTo = ., namenet = toupper(type)) %>%
                              filter(ids %in% orf_missing)
                            },
                            doing = 'retrieving centrality with min score 700')

  return(centrality_low)
}

fix_missing_centrality = function(df,id='ORF',col_prefix='cat_interactions.string.',taxon=4932){
  # Replace orf with missing values for centrality with 0's
  net_type = str_extract(tolower(col_prefix),'(string|intact)')
  df$has_unique_id = !is.na(df[[id]]) & !duplicated(df[[id]])
  df_unique = df %>% drop_na(id) %>% rowwise %>% filter(has_unique_id) %>% column_to_rownames(id)

  if(net_type=='string'){
    orf_missing =  df_unique %>% get_centrality_col(col_prefix) %>% find_na_rows() %>% rownames()
  }else if(net_type=='intact'){
    orf_missing = df_unique %>% get_centrality_col(col_prefix) %>% find_na_rows() %>% rownames()
  }

  if( length(orf_missing) == 0 ){
    .warn$log("No missing values found for any of the network centrality measures used!")
    return(df)
  }

  .warn$log("Replace columns with missing values for network centrality measures...\n")
  df_missing = tibble(ids = orf_missing)
  df_na_centrality = retrieve_missing_centrality(orf_missing,net_type,taxon) %>%
                     right_join(df_missing,by='ids') %>%
                     dplyr::rename_with(.cols=starts_with('cent_'), Pxx, px=col_prefix, s='')

  numVar = df_na_centrality %>% dplyr::select(where(is.numeric)) %>% as.matrix
  .info$log("Remaining missing values of centralities are replaced by imputation with 'missForest'!")
  library(missForest)
  df_na_centrality = missForest::missForest(numVar)$ximp %>% as_tibble() %>% mutate(ids=orf_missing)
  if(net_type=='intact'){
    df_na_centrality = df_na_centrality %>% dplyr::rename(UNIPROTKB=ORF)
    df_fixed = coalesce_join(x = df, y=df_na_centrality, by = 'UNIPROTKB')
  }else{
    if( taxon==9606){
      df = df %>% mutate(across(where(is.character) & starts_with('string.cent_'), as.numeric))
      df_na_centrality = df_na_centrality %>% rename(ensp=ids)
    }
    df_fixed = coalesce_join(x = df, y=df_na_centrality, by = c('ensp'))
  }

  return(df_fixed)
}


fix_missing_paxdb = function(df){
  .warn$log("Replace columns with missing paxdb ortholog abundance values...\n")

  paxdb = df %>% dplyr::select(contains('paxdb')) %>%
    mutate(cat_transcriptomics.paxdb.ortho_ppm_n=replace_na(cat_transcriptomics.paxdb.ortho_ppm_n,0)) %>%
    as.matrix()
  df_na_paxdb = missForest::missForest(paxdb)$ximp %>% as_tibble() %>% mutate(ORF=df$ORF)
  df_impute = coalesce_join(x = df, y=df_na_paxdb, by = 'ORF')
  return(df_impute)
}
### WORKFLOW PROCESSING MISSING VALUES
PROCESS_MISSING_VALUES = function(MAT, IDS, taxon=4932){

  miss.0=check_missing_var(MAT)
  # Predictors with missing values must be corrected
  #  I. Fix missing observations for certain genes (rows):
  #     A) missing codons counts
  PREDICTORS.1 = fix_missing_codons(MAT,col_prefix="cat_transcriptomics.sgd.")
  miss.1=check_missing_var(PREDICTORS.1)

  #     B) missing protein length/MW
  PREDICTORS.2 = fix_missing_peptide_stats(PREDICTORS.1,
                                           col_len='cat_transcriptomics.sgd.prot_size',
                                           col_mw='cat_biophysics.pepstats.mw',
                                           col_mw_avg='cat_transcriptomics.pepstats.mean_MW',
                                           col_charge='cat_biophysics.pepstats.netcharge',
                                           col_pi='cat_biophysics.pepstats.pI')
  miss.2=check_missing_var(PREDICTORS.2)

  #     C) missing centrality values (STRING and INTACT network are treated individually)
  PREDICTORS.3 = fix_missing_centrality(df=PREDICTORS.2,col_prefix="cat_interactions.string.",taxon)
  PREDICTORS.4 = fix_missing_centrality(PREDICTORS.3,col_prefix="cat_interactions.intact.",taxon) # only available for yeast
  #     D) missing protein disorder (D2P2/IUP)
  #     E) missing transcriptomics in barton 2010

  #  II.  Remove unnecessary variables (columns):
  #     A) rarely observed data (less than 2 observations)
  PREDICTORS.5 = PREDICTORS.4 %>% filter(ORF %in% IDS) %>% remove_rare_vars()
  miss.5 = check_missing_var(PREDICTORS.5)

  PREDICTORS.6 = fix_missing_paxdb(PREDICTORS.5)
  miss.6 = check_missing_var(PREDICTORS.6)

  PREDICTORS.7 = PREDICTORS.6
  .warn$log(sprintf("Replace NAs by fixed values for ohnologs, pangenome, limited-proteolysis data...\n"))
  PREDICTORS.7$cat_genomics.byrne2005.RLEN = replace_na(PREDICTORS.6$cat_genomics.byrne2005.RLEN,1)
  PREDICTORS.7$cat_genomics.byrne2005.PID1 = replace_na(PREDICTORS.6$cat_genomics.byrne2005.RLEN,100)
  PREDICTORS.7$cat_genomics.byrne2005.pid = replace_na(PREDICTORS.6$cat_genomics.byrne2005.pid,100)
  PREDICTORS.7$cat_genomics.byrne2005.score_B100 =  replace_na(PREDICTORS.6$cat_genomics.byrne2005.score_B100, median_(PREDICTORS.6$cat_genomics.byrne2005.score_B100))
  PREDICTORS.7$cat_genomics.byrne2005.S[is.na(PREDICTORS.6$cat_genomics.byrne2005.S)] = PREDICTORS.6$cat_genomics.sgd.len[is.na(PREDICTORS.6$cat_genomics.byrne2005.S)]
  PREDICTORS.7$cat_genomics.byrne2005.N = replace_na(PREDICTORS.6$cat_genomics.byrne2005.N,0)
  PREDICTORS.7$cat_genomics.byrne2005.G = replace_na(PREDICTORS.6$cat_genomics.byrne2005.G,0)
  PREDICTORS.7$cat_genomics.peter2018.f_gene_loss = replace_na(PREDICTORS.6$cat_genomics.peter2018.f_gene_loss,0)
  PREDICTORS.7$cat_genomics.peter2018.strains_presence = replace_na(PREDICTORS.6$cat_genomics.peter2018.strains_presence,1011)
  PREDICTORS.7$cat_biophysics.leuenberger2017.LIP_nres = replace_na(PREDICTORS.6$cat_biophysics.leuenberger2017.LIP_nres,0)
  PREDICTORS.7$cat_biophysics.leuenberger2017.LIP_npep = replace_na(PREDICTORS.6$cat_biophysics.leuenberger2017.LIP_npep,0)
  PREDICTORS.7$cat_biophysics.leuenberger2017.LIP_protein_coverage = replace_na(PREDICTORS.6$cat_biophysics.leuenberger2017.LIP_protein_coverage,0)

  # Replace missing values by the median
  median_aa_dom = PREDICTORS.6 %>% dplyr::select(ends_with("_dom")) %>% summarise(across(everything(),.fns=median_))
  median_aa_iup20 = PREDICTORS.6 %>% dplyr::select(ends_with("_iup20")) %>% summarise(across(everything(),.fns=median_))
  median_snp = PREDICTORS.6 %>% dplyr::select(starts_with("cat_genomics.peter2018.snp")) %>% summarise(across(everything(),.fns=median_))
  var_to_median = c(median_aa_dom,median_aa_iup20,median_snp)
  .warn$log(sprintf("Replace NAs by median for amino-acid average propensities and snp evolutionary rate (%s variables)...\n",length(var_to_median)))
  PREDICTORS.7 = replace_na(data = PREDICTORS.7, var_to_median )

  # HANDLING REMAINING MISSING VALUES (AUTOMATICALLY)
  var_has_na = PREDICTORS.7 %>% dplyr::select_if(~any(is.na(.)) )
  .warn$log(sprintf("Imputation with 'missForest' to replace missing values in the remaining %s variables...\n",length(var_has_na)))
  freq_na = apply(var_has_na,2,function(x){ mean_(is.na(x)) })
  var_with_na = PREDICTORS.7 %>% dplyr::select(names(var_has_na)) %>% as.matrix
  library(missForest)
  var_without_na = missForest::missForest(var_with_na)$ximp %>% as_tibble()
  PREDICTORS.7[,names(var_without_na)] = var_without_na

  var_na = PREDICTORS.7 %>% dplyr::select_if(~any(is.na(.)) )
  freq_na = apply(var_na,2,function(x){ mean_(is.na(x)) })
  if(length(freq_na) == 0){
    .succ$log("All missing values were successfully replaced!")
  }else{
    .warn$log(sprintf("Missing values remain in %s variables!",length(freq_na)))
    print(freq_na)
  }
  return(PREDICTORS.7)
}

# 3. PLOTS FIT ----------------------------------------------------------------
make_plot_1A = function(dat=EVOLUTION, X='PPM', Y="log10.EVO.FULL",
                        XLAB='Protein Abundance (log10)', YLAB='Evolutionary Rate Orthologs (log10)',
                        ANNOT=ANNOTATION, id=c('ORF','UNIPROT'),
                        col_lm = 'black', col_stat2d='black',
                        centerY=T,
                        fill_points='#CCCCCC', col_points='#CCCCCC', shape_points=19, stroke_points=1, size_points=1.5,
                        add_outliers=10,noplot=F,
                        add_smooth='lm',add_cor=T,show_cor='right',
                        x_as_exp10=F,y_as_exp10=F,
                        ymax,ymin){

  dat_annot = left_join(dat,ANNOT,by=id)
  yavg = mean_(dat_annot[[Y]])
  yavg_line = geom_hline(yintercept = yavg - (centerY*yavg), col=col_lm,linetype=2,size=0.5) # mean
  dat_annot[[Y]] = dat_annot[[Y]] - (yavg*centerY)
  OUTY = get_extremes(dat_annot,X,n=add_outliers)
  OUTX = get_extremes(dat_annot,Y,n=add_outliers)
  if(missing(ymin)){
    ymin = min_(dat_annot[[Y]])
  }
  if(missing(ymax)){
    ymax = max_(dat_annot[[Y]])
  }

  F1A=ggplot(dat_annot,aes_string(y=Y,x=X)) +
      th_txt_size + theme(aspect.ratio = 1) +
      ggpubr::grids(axis = 'xy') +
      scale_y_continuous(n.breaks = 10,limits = c(ymin,ymax))

  if(noplot){ return(F1A + yavg_line )}

  F1A = F1A +
    ggiraph::geom_point_interactive(aes(tooltip=FUNCTION, data_id=ORF),size=size_points,shape=shape_points,alpha=0.5,color=col_points,stroke=stroke_points) +
    stat_density2d(size=0.5,color=col_stat2d) +
    yavg_line +
    ylab(YLAB) + xlab(XLAB)

  F1A = F1A + geom_smooth(method=add_smooth,color=col_lm)
  if(add_cor){
    cor_xy = spearman.toplot(dat_annot[[X]],dat_annot[[Y]])
    if(show_cor=='right'){
     F1A = F1A + geom_text(data=cor_xy,aes(label=toshow),x=Inf,y=-Inf,hjust='inward',vjust='inward',size=3.5)
    }else if(show_cor=='left'){
      F1A = F1A + geom_text(data=cor_xy,aes(label=toshow),x=-Inf,y=-Inf,hjust='inward',vjust='inward',size=3.5)
    }else if(show_cor=='topleft'){
      F1A = F1A + geom_text(data=cor_xy,aes(label=toshow),x=-Inf,y=Inf,hjust='inward',vjust='inward',size=3.5)
    }else if(show_cor=='topright'){
      F1A = F1A + geom_text(data=cor_xy,aes(label=toshow),x=Inf,y=Inf,hjust='inward',vjust='inward',size=3.5)
    }else{
      F1A = F1A + geom_text(data=cor_xy,aes(label=toshow),x=Inf,y=Inf,hjust='inward',vjust='inward',size=3.5)
    }
  }

  if( x_as_exp10 ){
    F1A = F1A + scale_x_continuous(labels = scales::math_format(10^.x),breaks = 1:6)
  }
  if( y_as_exp10 ){
      F1A = F1A + scale_y_continuous(labels = scales::math_format(10^.x),breaks = 1:6)
  }

  if(add_outliers>0){
    F1A  = F1A +
      geom_text_repel(data=OUTX, aes(label = GENENAME),max.overlaps = 20,col='blue',size=TEXT_SIZE) +
      geom_text_repel(data=OUTY, aes(label = GENENAME),max.overlaps = 20,col='red',size=TEXT_SIZE) +
      geom_point(data=OUTX, col='blue',size=0.5) +
      geom_point(data=OUTY, col='red',size=0.5)
  }
  #F1A + facet_wrap(~FunCat)
  return(F1A)
}

make_plot_1B = function(clade1='schizo',clade2='sacch.wgd',control_var='ppm',
                        use_residuals=F,tolog=T,force_intercept=T, text_size=AXIS_TTITLE_SIZE){

  data2plot = load.clade(clade1,clade2)
  xvar='Kc1'
  yvar='Kc2'
  labx=paste('Kc',clade1)
  laby=paste('Kc',clade2)
  if(tolog){
    xvar=paste0(xvar,'.log10')
    yvar=paste0(yvar,'.log10')
    labx=paste(labx,'(log10)')
    laby=paste(laby,'(log10)')
  }

  ctrl1=paste0(control_var,"1.log10")
  ctrl2=paste0(control_var,"2.log10")
  if(use_residuals){
    cat(sprintf("Controlling for variable '%s'\n",control_var))
    cat(sprintf("Controlling branch length in %s by corresponding protein expression '%s'\n",clade1,ctrl1))
    first = data2plot %>%
      broom::augment_columns(x=lm(reformulate(response = xvar, termlabels = ctrl1, intercept = force_intercept),data=.)) %>%
      dplyr::rename(.resid.1 = .resid, .fitted.1=.fitted, .se.fit.1=.se.fit, .hat.1=.hat, .sigma.1=.sigma, .cooksd.1=.cooksd,.std.resid.1=.std.resid)

    cat(sprintf("Controlling branch length in %s by corresponding protein expression '%s'\n",clade2,ctrl2))
    data2plot = first %>%
      broom::augment_columns(x=lm(reformulate(response = yvar, termlabels = ctrl2, intercept = force_intercept),data=.)) %>%
      dplyr::rename(.resid.2 = .resid, .fitted.2=.fitted, .se.fit.2=.se.fit, .hat.2=.hat, .sigma.2=.sigma, .cooksd.2=.cooksd,.std.resid.2=.std.resid) %>%
      mutate( rY=.resid.2, rX=.resid.1)

    xvar='rX'
    yvar='rY'
    labx=paste('resid.',labx)
    laby=paste('resid.',laby)
  }


  all_defined = !is.na(data2plot[[xvar]]) & !is.na(data2plot[[yvar]]) & !is.na(data2plot[[ctrl1]])
  data2plot = data2plot[all_defined,]
  rho=spearman.toplot(data2plot[[xvar]],data2plot[[yvar]])

  F1B = ggplot(data=data2plot,aes_string(x=xvar,y=yvar)) +
    geom_point(size=0.7, alpha=0.8,shape=19) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour=NA,size=0.25,alpha=0.15,show.legend = F)+
    geom_smooth(method = 'lm',col='black',se=F,size=1) +
    geom_text(data=rho,aes(x=Inf,y=-Inf,label=toshow),hjust='inward',vjust='inward',size=4) +
    xlab(labx) + ylab(laby) +
    scale_fill_distiller(palette = 'Spectral',direction = -1) +
    ggpubr::grids() +
    theme(axis.title = element_text(size=text_size),
          axis.text = element_text(size=text_size))
  return(F1B)
}

control_var = function(input_data, target, control_var, suffix="", force_intercept=T){
  cat(sprintf("Controlling '%s' by '%s'\n",target,control_var))
  cols_fit = c('.resid','.fitted','.se.fit','.hat','.sigma','.cooksd','.std.resid')
  controlled = input_data %>%
    broom::augment_columns(x=lm(reformulate(response=target, termlabels=control_var, intercept=force_intercept),data=.)) %>%
    dplyr::rename_with(.cols = all_of(cols_fit),.fn = paste0, sep="", suffix)
  return(controlled)
}

make_plot_1C = function(data2plot=EVOLUTION,Y='EVO.FULL',X='SNP.FULL',control_var='PPM',
                        use_residuals=F,force_intercept=T, text_size=AXIS_TTITLE_SIZE,
                        col_lm = 'black', col_stat2d='black',
                        fill_points='#CCCCCC', col_points='#CCCCCC',
                        shape_points=19, stroke_points=1, size_points=1.5){

  if(use_residuals){
    first = control_var(data2plot, X, control_var, "_x", force_intercept)
    data2plot = control_var(first, Y, control_var, "_y", force_intercept)
    X='.resid_x'
    Y='.resid_y'
  }

  all_defined = !is.na(data2plot[[X]]) & !is.na(data2plot[[Y]]) & !is.na(data2plot[[control_var]])
  data2plot = data2plot[all_defined,]
  rho=spearman.toplot(data2plot[[X]],data2plot[[Y]])

  F1C = ggplot(data=data2plot,aes_string(x=X,y=Y)) +
    geom_point(size=size_points,shape=shape_points,alpha=0.5,color=col_points,stroke=stroke_points)+
    #geom_point(size=0.7, alpha=0.8,shape=19) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour=NA,size=0.25,alpha=0.15,show.legend = F)+
    geom_smooth(method = 'lm',col='black',se=F,size=1) +
    geom_text(data=rho,aes(x=Inf,y=-Inf,label=toshow),hjust='inward',vjust='inward',size=4) +
    scale_fill_distiller(palette = 'Spectral',direction = -1) +
    ggpubr::grids() +
    theme(axis.title = element_text(size=text_size),
          axis.text = element_text(size=text_size))
  return(F1C)
}


# 4. LINEAR FIT ----------------------------------------------------------------

get_XY_data = function(input,x=X,y=Y,noNA=T){
  # GETTING XY DATA FOR CURVE FITTING
  res=dplyr::lst(XX=input[[x]], YY=input[[y]],
                 varnames=c('x'=x,'y'=y),
                 n = c('x'=sum(!is.na(XX)), 'y' = sum(!is.na(YY)), 'xy' = sum(complete.cases(YY,XX)) ),
                 mu = c('x'=mean_(XX),'y'=mean_(YY)),
                 md = c('x'=median_(XX),'y'=median_(YY)),
                 var = c('x'=var_(XX),'y'=var_(YY))
  )

  res$df = input %>% ungroup()
  if( noNA ){
    df.noNA = res$df %>% dplyr::filter(complete.cases(!!sym(y),!!sym(x)))
    res$df = df.noNA
    res$YY=res$df[[y]]
    res$XX=res$df[[x]]

  }
  return(res)
}

get_model_params = function(m,x,y){
  # EXTRACT CURVE FITTING PARAMETERS
  nXY=sum( complete.cases(x,y) )
  mu.y=mean_(y)
  m.params = dplyr::lst(
    fit  = m,
    xx = x,
    yy = y,
    yfit = fitted(fit), # Y-Fitted
    yres = residuals(fit), # Y-Residual
    pfit = coef(fit), # Fitted parameters (intercept, PPM)
    TSS = sum( (yy-mu.y)^2 ),
    ESS = sum( (yfit-mu.y)^2 ), # Explained variance
    RSS = TSS-ESS, # Deviance (Unexplained variance)
    RSE = sqrt( (1 / (nXY-2)) * RSS ), # Residual standard error
  )
  return(m.params)
}

get_fit_data = function(d,m){
  # EXTRACT CURVE FITTING DATA (fitted/residuals quantities)
  mu.y = d$mu['y']
  var.y=d$var['y']
  nXY = d$n['xy']
  xname=d$varnames['x']
  yname=d$varnames['y']
  fit = d$df %>%
    broom::augment_columns(x=m) %>%
    mutate(ESS=sum_( (.fitted-mu.y)^2 ), TSS=sum_( (!!sym(yname)-mu.y)^2 ), RSS=TSS-ESS,
           s2=TSS/(nXY-1), s2.y = var.y, RS=sum(.resid),
           RSE=sqrt( (1 / (nXY-2)) * RSS), AIC = AIC(m), BIC=BIC(m), LL = readr::parse_number(as.character(logLik(m))))
  return(fit)
}


make_linear_fit = function(input,    # Input data
                           x=X, y=Y, # X/Y Variables to fit
                           only.params=T){
  # APPLY LINEAR REGRESSION (Y~X) ON INPUT DATA
  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']
  # model #
  model.name='Linear'
  f=as.formula(paste0(y,"~",x))
  m = lm(formula =f , data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model = model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

make_logistic_fit = function(input,    # Input data
                             y=Y, x=X, # X/Y Variables to fit
                             only.params=T){
  # APPLY LOGISTIC REGRESSION (Y~ 1/[1+e(-X)]) ON INPUT DATA
  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']

  # model #
  model.name='Sigmoid\n(logisitic)'
  f=as.formula(paste0(y,"~","SSlogis(",x,",Asym,xmid,scal)"))
  init.params = getInitial( f,data=xydata$df)
  m = nls(formula = f, start = init.params, data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model =model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

decompose_variance = function(MODEL,to.df=F){

  # DECOMPOSE VARIANCE FROM MODEL
  D = MODEL$model %>% as_tibble()

  N=sum(complete.cases(D))
  TSS = var( D[,1,drop=T] ) * (N-1)
  df.var = broom::tidy(aov(MODEL))

  df.ss = df.var %>% dplyr::filter(term != "Residuals")
  nvar = nrow(df.ss)
  df.res = df.var %>% dplyr::filter(term == "Residuals")
  RSS = sum(df.res$sumsq) # sum( residuals(LM)^2 )

  rss =  deviance(MODEL)
  if( !all.equal(RSS,rss) ){ warning(sprintf("RSS (%.2f) not equal to deviance of the model (%.2f)...",RSS,rss)) }

  ess = TSS-RSS
  ESS = sum(df.ss$sumsq)
  # Maximum ESS (Fitted - Ymean)
  ess.max = sum( (fitted(MODEL) - mean_(D[,1]))^2 )
  if( !all.equal(ESS,ess) ){ warning(sprintf("ESS (%.2f) not equal to sum of squares from individual variables (%.2f)...",ESS,ess)) }

  ESS.pc =100*ESS/TSS
  RSS.pc =100*RSS/TSS

  one_line_formula =  paste(deparse1(formula(MODEL))) %>%
                      str_trunc(side = 'center', width = 80)
  if( class(MODEL) == 'lm'){
    nterms = n_distinct(labels(MODEL))
  }else{
    nterms = n_distinct(attr(terms(MODEL),'term.labels'))
  }

  cat(sprintf("%s\n",one_line_formula))
  cat(sprintf("(%s predictor variables)\n",nterms))
  cat(sprintf("TSS %.1f (n=%s)\n",TSS,N))
  cat(sprintf("--> ESS %.1f (%.0f%%)\n",ESS, ESS.pc))
  cat(sprintf("--> RSS %.1f (%.0f%%)\n", RSS,  RSS.pc))
  if(to.df){
    res=tibble(N=N,nterms=nterms,TSS=TSS,ESS=ESS,RSS=RSS,
               RSS_rel = RSS.pc, ESS_rel=ESS.pc
               )
    return(res)
  }
}


fit_linear_regression = function(INPUT=EVOLUTION, X='PPM', Y="log10.EVO.FULL",
                                 PREDVAR=PREDICTORS, xcor_max=0.6, ycor_max=0.6,
                                 ADD.VARIABLES = 'log10.SNP.FULL',
                                 key_id=c("ORF","UNIPROT"),
                                 key_filter=c("IS_FUNGI","IS_STRAINS"),
                                 force_elim_x=T,
                                 min_obs=5){
  txt_section_break = repchar("-",50)
  #INPUT = EVOLUTION
  #Y = "log10.EVO.FULL" # mean Evolutionary rate (full sequence)
  #X = "MPC" # median Molecules Per Cell
  #X = "PPM" # Protein Abundance (log10 ppm) # ALTERNATIVELY

  LINREG = INPUT %>% dplyr::select(all_of(c(X,Y,ADD.VARIABLES)),any_of(c(key_id,key_filter)))

  XYDATA = get_XY_data(INPUT,x=X,y=Y)
  YY = XYDATA$YY
  XX = XYDATA$XX
  n.xy = XYDATA$n['xy']
  mu.y = XYDATA$mu['y']
  var.y = XYDATA$var['y']
  input = XYDATA$df
  rg.y = range_(YY)

 #starting(colnames(PREDVAR),'cat_')

  Y_X=make_linear_fit(LINREG ,x=X, y=Y, only.params = F)
  var_names = PREDVAR %>% dplyr::select(where(is.numeric) | where(is.logical), -any_of(colnames(Y_X))) %>% colnames
  M0=left_join(Y_X,PREDVAR)
  M0$yavg   = mu.y
  M0$yvar   = var.y
  M0$ymin = rg.y[1]
  M0$ymxa = rg.y[2]
  M0$nxy = n.xy


  # excluding variables with correlation to X/Y
  cor_x = cor(x=M0 %>% dplyr::select(all_of(var_names)),y=M0[[X]],use = 'pairwise') %>% as_vector %>% na.omit()

  x_control_var = rownames(cor_x)[abs(cor_x) >= xcor_max]
  cat("controlling variables correlated to ",X,"...\n")
  x_excluded_var=c()
  for(xx in x_control_var){
    xx_formula = reformulate(termlabels=X, response=xx)
    m_xx = lm(formula = xx_formula,data=M0)
    xx_resid =residuals(m_xx)
    cor_xx = cor(xx_resid, M0[[X]], use="pairwise")
    cor_yy = cor(xx_resid, M0[[Y]], use="pairwise")

    if(abs(cor_xx) < xcor_max & abs(cor_yy) > 0.1 & !force_elim_x){
      all_cor=sprintf("r(%s) = %.2f vs. r(resid_%s_%s) = %.2f and r(%s) = %.2f\n",X,cor_x[xx,1],X,xx,cor_xx,Y,cor_yy)
      cat(all_cor)
      M0[,paste0("resid_",X,"_",xx)] = xx_resid
    }else{
      x_excluded_var=c(x_excluded_var,xx)
    }
  }

  n_xout = length(x_excluded_var)
  cat(sprintf("Excluding %s predictors with cor. to X > %s (%s)\n",n_xout,xcor_max,X))
  cor_y = cor(x=M0[,var_names],y=M0[[Y]],use = 'pairwise') %>% as_vector %>% na.omit()

  y_excluded_var = rownames(cor_y)[abs(cor_y) >= ycor_max]
  n_yout = length(y_excluded_var)
  cat(sprintf("Excluding %s predictors with cor. to Y > %s (%s)\n",n_yout,ycor_max,Y))

  excluded_var = unique(c(x_excluded_var,y_excluded_var))
  n_out = length(excluded_var)
  cat(sprintf("In total, %s predictors are excluded:\n",n_out))
  cat(txt_section_break)
  cat(sprintf("\n%3s. %s",1:n_out,excluded_var))
  cat("\n",txt_section_break,"\n")

  return(M0 %>% dplyr::select(-all_of(excluded_var)))
}

# CHECK EACH VARIABLE
fit_lm_one_var = function(varname,target='.resid',inputdata, verbose=F, .pb=NULL, .pb.toprint=NULL){
  if(verbose){ catn(varname) }
  if(!is.null(.pb) && !.pb$finished){ .pb$tick(tokens=list(what=.pb.toprint)) }
  shortname = str_split_fixed(pattern='\\.',varname,n = 2)[,2]
  not_na = !is.na(inputdata[,varname])
  not_0  = not_na & inputdata[,varname] != 0
  n_var = sum(not_na & not_0)

  formula_var = paste0(target,' ~ ',varname)
  linreg = lm(formula_var, data=inputdata[not_na,])
  #decompose_variance(linreg)
  N=sum(complete.cases(linreg$model))
  TSS = var( linreg$model[,target] ) * (N-1)
  RSS = deviance(linreg)
  ESS = TSS - RSS
  rss_rel = 100*RSS/TSS
  ess_rel = 100*ESS/TSS

  #if(is.na(n_var) | n_var<10){
  #  return(NULL)
    # return(
    #   tibble(var=shortname, nprot_var=n_var,
    #          tss=TSS,rss=RSS,ess=ESS,
    #          pc_ess=ess_rel, pc_rss=rss_rel,
    #          tss_var=NA,rss_var=NA,ess_var=NA,
    #          pc_ess_var=NA, pc_rss_var=NA)
    # )
  #}

  linreg_var = lm(formula_var, data=inputdata[not_na & not_0,])
  #decompose_variance(linreg_var)
  N_var=sum(complete.cases(linreg_var$model))
  TSS_var = var( linreg$model[not_na & not_0,target] ) * (n_var-1)
  RSS_var = deviance(linreg_var)
  ESS_var = TSS_var - RSS_var

  rss_var = 100*RSS_var/TSS_var
  ess_var = 100*ESS_var/TSS_var
  # mu.y = mean_(inputdata[,target][has_var])
  # .fitted.prop = inputdata$`.fitted`[has_var] + mu.y
  # .resid.prop = inputdata[,target][has_var] + mu.y

  return(tibble(var=shortname, nprot_var=n_var,
                tss=TSS,rss=RSS,ess=ESS,
                pc_ess=ess_rel, pc_rss=rss_rel,
                tss_var=TSS_var,rss_var=RSS_var,ess_var=ESS_var,
                pc_ess_var=ess_var, pc_rss_var=rss_var))
                #pc_tss_var = tss_rel,pc_rss_var=rss_rel,))
}



# 5. NULL-MODEL AND VARIABLE SELECTION -----------------------------------------

fit_m0 = function(INPUT_LM,XCOL,YCOL,PREDICTORS,ZCOL,IDCOLS,
                  MAX_XCOR = 0.6, MAX_YCOR=0.8, MIN_N = 1){

  # Remove rare and low-variance predictors/data
  LMDATA = INPUT_LM %>% remove_rare_vars()
  PRED = PREDICTORS %>% remove_rare_vars()

  m0_vars = c(".fitted",".se.fit",".resid",".hat",".sigma",".cooksd",".std.resid","ESS","TSS","RSS","s2","s2.y","RS","RSE","AIC","BIC","LL","model")
  xy_vars = c("yavg","yvar","ymin","ymxa","nxy")
  m0 = fit_linear_regression(INPUT=LMDATA, X=XCOL,Y=YCOL, PREDVAR=PRED,
                             ADD.VARIABLES=ZCOL,key_id = IDCOLS,
                             xcor_max = MAX_XCOR,ycor_max = MAX_YCOR, min_obs=MIN_N, force_elim_x = T
  )


  pred_vars = m0 %>% dplyr::select(where(is.numeric) | where(is.logical), -any_of(c(m0_vars,xy_vars))) %>%
    colnames %>% setdiff(NA)
  if( str_detect( hutils::longest_prefix(pred_vars), 'cat_') ){
    pred_vars = str_subset(pred_vars,'^cat_')
  }

  predictors = m0 %>%
    dplyr::select(all_of(c(IDCOLS,m0_vars,XCOL,YCOL, ZCOL)), all_of(pred_vars))%>%
    mutate(across(.cols = where(is.numeric), .fns=~na_if(., is.infinite(.)))) %>%
    distinct()


  npred = pred_vars %>% n_distinct()
  ngene = predictors[,IDCOLS] %>% distinct() %>% nrow
  message("Number of genes | Number of predictors")
  print(c(ngene,npred))
  message("Type of predictor variables:")
  print(sapply(predictors[pred_vars],class) %>% janitor::tabyl())

  # LINEAR MODEL M0
  predictors = predictors %>% dplyr::filter( !is.na(YCOL) & !is.na(XCOL) )
  if(!is.null(ZCOL)){ predictors %>% dplyr::filter( !is.na(ZCOL) )  }

  formula_M0=reformulate(response = YCOL, termlabels = XCOL, intercept = T )
  LM0 = lm(data=predictors, formula_M0)
  df0=decompose_variance(LM0,T)

  coef0=coef(LM0)
  #pred_vars = colnames(predictors) %>% str_extract("cat_.+")
  x0 = sprintf("offset(%s*%s)",coef0[2],names(coef0)[2])

  if(!is.null(ZCOL)){
    formula_SNP=reformulate(response = YCOL, termlabels = ZCOL, intercept = T )

    LM_SNP = lm(data=predictors,formula_SNP)
    decompose_variance(LM_SNP)
    coef_snp=coef(LM_SNP)

    predictors_snp = predictors

    formula_M0_SNP =  reformulate(response = YCOL, termlabels =c(XCOL,ZCOL), intercept = T )
    LM0_SNP = lm(data=predictors_snp,formula_M0_SNP)
    decompose_variance(LM0_SNP)
    coef0_snp=coef(LM0_SNP)
    x0_snp = paste0( sprintf("offset(%s*%s)",coef0_snp[2:3],names(coef0_snp)[2:3]), collapse=' + ')
  }
  return(lst(P=predictors,P_vars=pred_vars,m0=m0,LM0=LM0))
}

select_variable = function(fitted_data=fit0, response='.resid',
                           raw=F, min_ess = 1, min_ess_frac = 1){

  pred_vars= fitted_data$P_vars %>% setdiff(NA)
  n_pred_vars = pred_vars  %>%  n_distinct()

  message(sprintf('testing each of the %s predictors for the correlation with y-residuals (%s)...',n_pred_vars,response))
  task='test variable prediction...'
  tic(msg = task)
  pb = progress::progress_bar$new(total = n_pred_vars, width = 100, clear=T,
                                  format = " (:spin) :what [:bar] :percent (:current/:total # :elapsed eta: ~:eta)")

  TSS = var( fitted_data$LM0$model[[YCOL]] ) * (nrow(fitted_data$LM0$model)-1)
  RSS = deviance(fitted_data$LM0)
  ESS = TSS-RSS
  mu_y  = mean_(fit0$P[[YCOL]])
  #response='ER'

  lm_single=tibble(variable=pred_vars,tss0=TSS,rss0=RSS,ess0=ESS) %>%
    group_by(variable) %>%
    mutate( fit_lm_one_var(inputdata = fitted_data$P , variable,response,.pb=pb, .pb.toprint=task) ) %>%
    ungroup() %>%
    mutate(pc_ess0 = 100*ess_var / tss0,
           pc_rss0 = 100*rss_var / tss0)

  # res=list()
  # for(i in 4:nrow(lm_single)){
  #    res[[i]] = fit_lm_one_var(inputdata = fitted_data$P , varname = lm_single$variable[i], target = response)
  #    if(i %% 100 == 0){
  #      print(i)
  #    }
  # }

  var_name = str_split_fixed(pattern="[_\\.]",lm_single$variable,n=4)
  lm_single =  lm_single %>% ungroup %>%
    mutate(type=var_name[,2],
           source=var_name[,3],
           varname=str_replace_all(var_name[,4],"[-_\\.]"," "))

  if(raw){
    warning("returning all tested variable!")
    return(lm_single)
  }
  best = lm_single %>% dplyr::filter(pc_ess >  min_ess_frac & ess > min_ess)
  return(best)
}


# 6. ELASTIC REGRESSION --------------------------------------------------------
fit_elastic = function(P=predictors, ftrain=0.66, kfold=100,
                       target=".resid", add_var = NULL,
                       alphas= seq(0,1,by=0.01), best="mse_l_1se"){
  set.seed(42)  # Set seed for reproducibility
  n=nrow(P)
  train_rows = sample(n,size=ftrain*n)
  cat("Fitting elastic net regression with 10-fold cross validation...\n")
  RESPONSE = P[[target]]
  if(!is.null(add_var)){
    cat("Using additional variables:\n",paste0(add_var,"\n"))
  }
  VARIABLES = P %>% dplyr::select(starts_with("cat_"),add_var)

  # training set
  x.train = VARIABLES[train_rows,]
  y.train = RESPONSE[train_rows]
  # test set
  x.test = VARIABLES[-train_rows,]
  y.test = RESPONSE[-train_rows]

  cat("------------------------\n")
  cat(sprintf("TARGET : %s\n",target))
  cat(sprintf("NOBS : %s (%s train %s test)\n",n, length(train_rows), n-length(train_rows)))
  cat(sprintf("# VARIABLES : %s\n",ncol(VARIABLES)))
  cat("------------------------\n")


  library(glmnet)
  cat(sprintf("Testing %s values for parameter alpha (lasso/ridge)...\n",length(alphas)))
  tic(sprintf("Tuning alpha parameters on elastic model (%s-CV)",kfold))
  alpha_fits = lapply(alphas,function(a){
    fit = cv.glmnet(as.matrix(x.train), as.matrix(y.train), type.measure="mse", alpha=a, family="gaussian",trace.it = F)
    return(fit)
  })

  results = purrr::map_dfr(1:101,
    function(i){
      y.fitted.lmin = predict(alpha_fits[[i]],s=alpha_fits[[i]]$lambda.min, newx=as.matrix(x.test))
      y.fitted.l1se = predict(alpha_fits[[i]],s=alpha_fits[[i]]$lambda.1se, newx=as.matrix(x.test))
      mse_lmin <- mean_((y.test - y.fitted.lmin)^2)
      mse_l1se <- mean_((y.test - y.fitted.l1se)^2)
      return(c(alpha=alphas[i],mse_l_1se=mse_l1se,mse_l_min=mse_lmin))
  })
  toc()

  ## View the results
  # Best fit (=with lowest mean square error) parameters
  lowest_mse = results[which.min(results[[best]]),best]
  best_alpha = alpha_fits[[which.min(results[[best]])]]
  alpha_lowest_mse = alphas[which.min(results[[best]])]
  lambda_lowest_mse = best_alpha$lambda.min

  cat("Best fit parameters:\n")
  cat("------------------------\n")
  cat(sprintf(" - min. mse = %.4f (lowest mean square error)\n",lowest_mse))
  cat(sprintf(" - alpha    = %.4f (ridge/lasso)\n",alpha_lowest_mse))
  cat(sprintf(" - lambda   = %.4f (shrinkage)\n",lambda_lowest_mse))

  y.trained = predict(best_alpha,s=best_alpha$lambda.min,newx=as.matrix(x.train))
  r2_train = cor(y.train,y.trained)^2

  y.predicted = predict(best_alpha,s=best_alpha$lambda.min,newx=as.matrix(x.test))
  r2_test = cor(y.test,y.predicted)^2

  cat("\nR2 (explained variance):\n")
  cat("------------------------\n")
  cat(sprintf(" - training = %.1f%%\n",100*r2_train))
  cat(sprintf(" - test     = %.1f%%\n",100*r2_test))

  best_elastic <- glmnet(as.matrix(VARIABLES), RESPONSE, alpha = alpha_lowest_mse, lambda = lambda_lowest_mse, standardize = TRUE,trace.it = F)
  explicative_var = best_elastic$beta[abs(best_elastic$beta[,1])>0,]
  npred_lowest_mse=length(explicative_var)
  pc_ess_lowest_mse = 100*best_elastic$dev.ratio
  rss_lowest_mse = deviance(best_elastic)
  pc_rss_lowest_mse = 100*(1-best_elastic$dev.ratio)
  cat("------------------------\n")
  cat(sprintf("Best elastic model explains %.2f%% of the total variance (%.1f)...\n",pc_ess_lowest_mse,best_elastic$nulldev))
  cat(sprintf("It has %s explicative variables...\n",best_elastic$df))
  cat(sprintf("remaining unexplained variance: RSS = %.2f%% (%.1f)\n",pc_rss_lowest_mse,rss_lowest_mse))

  return(list(all=alpha_fits,min_mse=best_alpha, best=best_elastic,alpha_min_mse=alpha_lowest_mse))
}

# 7. EVALUATE PERFORMANCE OF PREDICTION  _____----------------------------------

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  NVAR = nrow(df)
  # Model performance metrics
  perf = data.frame( RMSE = RMSE, r2 = R_square, ESS=SSE, TSS=SST, nvar=NVAR)
  return(perf)
}

aov_model = function(model,min_pv=1e-10){
  aov_data = tidy(aov(model)) %>%
    mutate(ss_tot=sum(sumsq), ss_rel = sumsq/ss_tot) %>%
    mutate(variable = str_replace(term,"cat_[^\\.]+\\.",""))

  aov_sig = aov_data %>%
    dplyr::filter(p.value < min_pv) %>%
    arrange(ss_rel) %>% mutate(ess=cumsum(ss_rel))
  aov_sig$nvar=nrow(aov_sig)
  aov_sig$max_ess = 100*max(aov_sig$ess)
  return(aov_sig)
}

aov_plot=function(aov_data, name){
  library(tidytext)
  data2plot = aov_data %>%
    mutate(input=as.factor(name),
           ord_var=reorder_within(x = variable, by=ess, within = input)
    )
  max_ess = 100*max(data2plot$ess)
  nvar=n_distinct(data2plot$ord_var)
  p = ggplot(data2plot,aes(y=ord_var,x=ess)) +
    geom_vline(xintercept=seq(0,0.2,by=0.01),linetype='dotted',size=0.5,col='gray')+
    geom_point(col='red',size=0.5) +
    geom_line(aes(x=ess,group=1),col='red') +
    geom_bar(aes(x=ss_rel),stat='identity') +
    theme(axis.text = element_text(size=8),legend.position = 'top')+
    scale_x_continuous(minor_breaks = c(0.05,0.1,0.15),breaks = seq(0,0.8,by=0.1)) +
    scale_y_reordered() +
    ylab('variables (ordered by %SS)') + xlab('Sum of Squares (%SS)') +
    ggtitle(sprintf('%s variables - %%ESS %.1f',nvar,max_ess)) +
    facet_wrap(~input,scales='free_y')

  plot(p)
  return(p)
}

elastic_dev = function(fit_elastic){

  alpha = fit_elastic$alpha_min_mse
  lambda = log(min(fit_elastic$min_mse$lambda))
  df.glmnet = tibble(nvar=fit_elastic$min_mse$glmnet.fit$df,
                     explained.var=fit_elastic$min_mse$glmnet.fit$dev.ratio,
                     l=fit_elastic$min_mse$glmnet.fit$lambda)
  df.best= tibble(nvar=fit_elastic$best$df,explained.var=fit_elastic$best$dev.ratio)

  p_elastic = ggplot(df.glmnet,aes(y=100*explained.var,x=nvar)) +
    geom_line() + geom_point(size=2) +
    #geom_point(data=df.best,col='red') +
    scale_x_continuous(breaks=seq(0,120,by=5)) +
    scale_y_continuous(breaks=seq(0,80,by=5),limits = c(0,80))+
    ylab('Explained variance ESS (%)') + xlab('# of variables') +
    ggtitle(label=sprintf('Elastic Net alpha=%.2f lambda=%.1f',alpha,lambda),
                        subtitle=sprintf("%s variables - %%ESS %.1f",df.best$nvar,  100*df.best$explained.var))
  p_elastic
  return(p_elastic)
}

