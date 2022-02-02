source(here::here("src","__setup_yeastomics__.r"))
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

  annotation = full_join(sgd_desc,uni_feat,by=c("SGD","UNIPROT"='UNIPROTKB')) %>%
               full_join(biofunc,by='ORF')  %>%
               relocate(SGD,GENENAME,ORF,UNIPROT,PNAME,L,FAMILIES,FUNCTION,ROLE,
                        BIOPROCESS_all,LOC,COMPLEX,ORTHO,OTHER,KEYWORDS,
                        EXISTENCE,SCORE) %>%
               filter(!is.na(ORF) | is.na(UNIPROT))

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
    filter(!(is.dup(UNIPROT) & is.na(EVO.FULL))) %>%
    #dplyr::mutate(across(starts_with("EVO."),function(x){ x/mean_(x) },.names="norm.{.col}")) %>% # scale and center all evo rate
    dplyr::mutate(across(starts_with("log10.EVO."),log10,.names = "log10.{.col}")) # # apply log10 to rate4site
    return(fungi)
}

load.strains.evo = function(){
  ### STRAINS CONTAIN ALL DATA ABOUT THE FUNGI LINEAGE EVOLUTIONARY RATE AND THE S.CEREVISIAE POPULATION
  strains = readRDS(here("data","PROTEIN-EVO-FUNGI-SNP.rds")) %>%
    dplyr::select(c(starts_with(c('EVO.','SNP.')),'PPM','ORF','UNIPROT','IS_FUNGI','IS_STRAINS')) %>%
    filter(!(is.dup(UNIPROT) & is.na(EVO.FULL))) %>%
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

normalize_features=function(feat){

  # FEATURES NORMALIZATION
  feat_norm = feat

  # a. Normalize codons counts (0-100) => 64 COLUMNS ---------------------------
  col_codons = grep("cat_transcriptomics.sgd.[ATCG]{3}",colnames(feat))
  feat_norm[,col_codons] =  100 * feat[,col_codons] / (feat$cat_transcriptomics.sgd.prot_size+1)

  # b. Normalize amino acid frequencies (0-100) => 31 COLUMNS  -----------------
  # 04.01.22 changed to use sgd sequences instead of uniprot
  col_f_aa = grep("cat_biophysics.sgd.f_",colnames(feat))
  feat_norm[,col_f_aa] =  100 * feat[,col_f_aa]

  # c. Normalize all other fraction (0-100) => 4 COLUMNS  ----------------------
  col_frac = c("cat_genomics.sgd.pGC","cat_genomics.byrne2005.RLEN","cat_transcriptomics.coRdon.CU_fop",
               "cat_biophysics.d2p2.f")
  #"cat_biophysics.dubreuil2019.IUP20_f","cat_biophysics.dubreuil2019.IUP30_f","cat_biophysics.dubreuil2019.IUP40_f")
  feat_norm[,col_frac] =  100 * feat[,col_frac]

  # d. Normalize count variables => 19 COLUMNS ---------------------------------
  col_count = c("cat_transcriptomics.paxdb.ppm_4932", "cat_transcriptomics.paxdb.ppm_214684","cat_transcriptomics.paxdb.ppm_4896","cat_transcriptomics.paxdb.ppm_5061",
                "cat_transcriptomics.paxdb.ortho_ppm_avg","cat_transcriptomics.paxdb.ortho_ppm_sd",
                "cat_transcriptomics.paxdb.ortho_ppm_max","cat_transcriptomics.paxdb.ortho_ppm_min",
                "cat_interactions.string.cent_pagerank",
                "cat_interactions.string.cent_eigen","cat_interactions.string.cent_authority",
                "cat_interactions.string.cent_hub","cat_interactions.string.cent_subgraph","cat_interactions.intact.cent_deg",
                "cat_interactions.intact.cent_pagerank","cat_interactions.intact.cent_eigen",
                "cat_interactions.intact.cent_authority",
                "cat_interactions.intact.cent_hub","cat_interactions.intact.cent_subgraph")
  feat_norm[,col_count] =  log10(feat[,col_count]+1)

  # e. Normalize count variables => 13 COLUMNS ---------------------------------
  col_doublings = c("cat_genomics.sgd.len","cat_transcriptomics.sgd.prot_size",
                    "cat_genomics.byrne2005.S","cat_genomics.byrne2005.N","cat_genomics.byrne2005.G",
                    "cat_transcriptomics.geisberg2014.HL_mrna", "cat_biophysics.villen2017.HL_prot",
                    "cat_biophysics.d2p2.nseg","cat_biophysics.d2p2.L","cat_biophysics.d2p2.Lsegmax",
                    "cat_biophysics.dubreuil2019.IUP20_L","cat_biophysics.dubreuil2019.IUP30_L","cat_biophysics.dubreuil2019.IUP40_L")
  feat_norm[,col_doublings] =  log2(feat[,col_doublings]+1)

  return(feat_norm)
}


# 2. FIX MISSING VALUES --------------------------------------------------------
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
  rare_vars = binary_vars[ colSums(binary_vars) < min_obs ]
  n_rare = length(rare_vars)
  cat(sprintf("Excluding %s/%s predictors with less than %s observations\n",n_rare,ncol(df),min_obs))
  df_fixed = df %>% dplyr::select(-all_of(names(rare_vars)))
  return(df_fixed)
}


#### B. IN ROWS ####
#### 1. codons counts ####
get_codons_col = function(df,col_prefix='cat_transcriptomics.sgd.'){
  regex_codons=paste0(col_prefix,get.codons4tai(),"$")
  res = df %>% dplyr::select(matches(regex_codons,ignore.case = F))
  return(res)
}

retrieve_missing_codons = function(orf_missing){
  # Retrieve CSD and count codons for orf with missing values
  cod = paste0(get.codons4tai(),"$")
  cds = load.sgd.CDS()
  cat(sprintf("Retrieving missing codons counts for %s ORF...\n",n_distinct(orf_missing)))
  codon_counts = coRdon::codonCounts(coRdon::codonTable(cds[orf_missing]))
  codon_missing = tibble(ORF=orf_missing, as_tibble(codon_counts))
  return(codon_missing)
}

fix_missing_codons = function(df,col_prefix='cat_transcriptomics.sgd.'){
  # Replace orf with missing values with retrieved codons counts from CDS
  orf_missing = df %>% column_to_rownames('ORF') %>% get_codons_col(col_prefix) %>% find_na_rows() %>% rownames()
  cat("Replace columns of codons counts with missing values...\n")
  df_na_codon = retrieve_missing_codons(orf_missing) %>%
    dplyr::rename_with(.cols=matches(get.codons4tai(),"$"),.fn=Pxx, px=col_prefix, s='')
  df_fixed = coalesce_join(x = df, y=df_na_codon, by = "ORF")
  return(df_fixed)
}

#### 2. protein length / Average Molecular Weight ####
fix_missing_peptide_stats = function(df,
                                     col_len='cat_transcriptomics.sgd.prot_size',
                                     col_mw='cat_biophysics.pepstats.mw',
                                     col_mw_avg='cat_transcriptomics.pepstats.mean_MW',
                                     col_charge='cat_biophysics.pepstats.netcharge',
                                     col_pi='cat_biophysics.pepstats.pI'
                                     ){
  col_pep = c(col_len,col_mw,col_mw_avg,col_charge,col_pi)
  # Replace orf with missing values for average molecular weight
  orf_missing = df %>% column_to_rownames('ORF') %>% dplyr::select(all_of(col_pep)) %>% find_na_rows() %>% rownames()
  if(length(orf_missing)>1){
    cat("Replace columns with missing values for protein length / average molecular weight...\n")
    prot_missing = load.sgd.proteome()[orf_missing] %>% as.character

    library(Peptides)
    df_na_pep = tibble(orf_missing,
                      missing_len = Peptides::lengthpep(prot_missing),
                      missing_mw  = Peptides::mw(prot_missing),
                      missing_mw_avg = missing_mw / missing_len,
                      missing_charge = Peptides::charge(prot_missing),
                      missing_pi = Peptides::pI(prot_missing),
                      cat_transcriptomics.pepstats.AA_costly = missing_mw_avg > 118,
                      cat_transcriptomics.pepstats.AA_cheap = missing_mw_avg <= 105) %>%
      dplyr::rename(ORF := orf_missing,
                    !!col_len:=missing_len,
                    !!col_mw:=missing_mw,
                    !!col_mw_avg:=missing_mw_avg,
                    !!col_charge:=missing_charge,
                    !!col_pi:=missing_pi
                    )

    df_fixed = coalesce_join(x = df, y=df_na_pep, by = "ORF")
    return(df_fixed)
  }else{
    warning("No missing values found for any of the peptide stats columns!")

    return(df)
  }
}

#### 3. network centrality ####
get_centrality_col = function(df,col_prefix="cat_interactions.string."){
  regex_centrality = paste0("^",col_prefix,"cent_")
  res = df %>% dplyr::select(matches(regex_centrality,ignore.case = F))
  return(res)
}

retrieve_missing_centrality = function(orf_missing,type='string'){
  interactions=load.network(type)
  cent=interactions %>%
        filter(ORF1 %in% orf_missing | ORF2 %in% orf_missing) %>%
        dplyr::select(ORF1,ORF2) %>%
        network.centrality(fromTo = ., namenet = toupper(type)) %>%
        filter(ids %in% orf_missing)
  return(cent)
}

fix_missing_centrality = function(df,col_prefix='cat_interactions.string.',NA_default_val=0){
  # Replace orf with missing values for centrality with 0's
  net_type = str_extract(col_prefix,'(string|intact)')
  if(net_type=='string'){
    orf_missing = df %>% column_to_rownames('ORF') %>% get_centrality_col(col_prefix) %>% find_na_rows() %>% rownames()
  }else if(net_type=='intact'){
    orf_missing = df %>% filter(!is.dup(UNIPROTKB)) %>% column_to_rownames('UNIPROTKB') %>% get_centrality_col(col_prefix) %>% find_na_rows() %>% rownames()
  }

  df_missing = tibble(ORF = orf_missing)

  df_na_centrality = retrieve_missing_centrality(orf_missing,net_type) %>%
                     dplyr::rename(ORF=ids) %>%
                     right_join(df_missing,by='ORF') %>%
                     dplyr::rename_with(.cols=starts_with('cent_'), Pxx, px=col_prefix, s='')
  # Persistent NAs are modified to a predefined value (e.g. 0)
  warning(sprintf("Remaining missing centrality measures are replaced by %s\n",NA_default_val))
  df_na_centrality[is.na(df_na_centrality)] = NA_default_val

  cat("Replace columns of network centrality with missing values...\n")

  if(net_type=='intact'){
    df_na_centrality = df_na_centrality %>% dplyr::rename(UNIPROTKB=ORF)
    df_fixed = coalesce_join(x = df, y=df_na_centrality, by = 'UNIPROTKB')
  }else{
    df_fixed = coalesce_join(x = df, y=df_na_centrality, by = 'ORF')
  }
  return(df_fixed)
}

# 3. LINEAR FIT ----------------------------------------------------------------

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

decompose_variance = function(LM){
  # DECOMPOSE VARIANCE FROM LINEAR REGRESSION

  N=sum(complete.cases(LM$model))
  TSS = var( LM$model[,1] ) * (N-1)
  df.var = summary(aov(LM))[[1]]
  nvar = nrow(df.var)
  RSS = sum(df.var$`Sum Sq`[nvar])
  ESS = sum(df.var$`Sum Sq`[-nvar])

  ess.max = sum( (fitted(LM) - mean_(LM$model[,1]))^2 )
  ess = TSS-RSS
  #RSS =  deviance(LM)
  #rss = sum( residuals(LM)^2)
  #ESS =  TSS - RSS
  cat(sprintf("TSS %.1f (n=%s)\n",TSS,N))
  cat(sprintf("--> ESS total %.1f (%.0f%%) max. %.1f ==> gain=%.1f (+%.0f%%)\n",ess, 100*ess/TSS, ess.max, ESS, 100*ESS/TSS))
  cat(sprintf("--> RSS %.1f (%.0f%%)\n", RSS,  100*RSS/TSS))
}


fit_linear_regression = function(INPUT=EVOLUTION, X='PPM', Y="log10.EVO.FULL",
                                 PREDVAR=PREDICTORS, xcor_max=0.6, ycor_max=0.6,
                                 min_obs=5){
  txt_section_break = repchar("-",50)
  #INPUT = EVOLUTION
  #Y = "log10.EVO.FULL" # mean Evolutionary rate (full sequence)
  #X = "MPC" # median Molecules Per Cell
  #X = "PPM" # Protein Abundance (log10 ppm) # ALTERNATIVELY
  key_id=c("ORF","UNIPROT")
  key_filter=c("IS_FUNGI","IS_STRAINS")

  LINREG = INPUT %>% dplyr::select(all_of(c(key_id,key_filter,X,Y)))

  XYDATA = get_XY_data(INPUT,x=X,y=Y)
  YY = XYDATA$YY
  XX = XYDATA$XX
  n.xy = XYDATA$n['xy']
  mu.y = XYDATA$mu['y']
  var.y = XYDATA$var['y']
  input = XYDATA$df
  rg.y = range_(YY)
  var_names = starting(colnames(PREDVAR),'cat_')

  Y_X=make_linear_fit(LINREG ,x=X, y=Y, only.params = F)
  M0=left_join(Y_X,PREDVAR)
  M0$yavg   = mu.y
  M0$yvar   = var.y
  M0$ymin = rg.y[1]
  M0$ymxa = rg.y[2]
  M0$nxy = n.xy

  # excluding variables with correlation to X/Y
  cor_x = cor(x=M0[,var_names],y=M0[[X]],use = 'pairwise') %>% as_vector %>% na.omit()
  x_excluded_var = rownames(cor_x)[abs(cor_x) >= xcor_max]
  n_xout = length(x_excluded_var)
  cat(sprintf("Excluding %s predictors with cor. to X > %s (%s)\n",n_xout,xcor_max,X))
  cor_y = cor(x=M0[,var_names],y=M0[[Y]],use = 'pairwise') %>% as_vector %>% na.omit()

  y_excluded_var = rownames(cor_y)[abs(cor_y) >= ycor_max]
  n_yout = length(y_excluded_var)
  cat(sprintf("Excluding %s predictors with cor. to Y > %s (%s)\n",n_yout,ycor_max,Y))

  excluded_var = unique(c(x_excluded_var,y_excluded_var))
  n_out = length(excluded_var)
  cat(sprintf("In total, %s predictors are excluded:\n",n_out))
  cat(txt_section_break,"\n ")
  cat(sprintf("%3s. %s\n",1:n_out,excluded_var))

  return(M0 %>% dplyr::select(-all_of(excluded_var)))
}

# 4. PLOTS FIT ----------------------------------------------------------------
make_plot_1A = function(dat=EVOLUTION, X='PPM', Y="log10.EVO.FULL",
                        ANNOT=ANNOTATION, id=c('ORF','UNIPROT'),
                        add_outliers=10){
  dat_annot = left_join(dat,ANNOT,by=id)
  OUTY = get_extremes(dat_annot,X,n=add_outliers)
  OUTX = get_extremes(dat_annot,Y,n=add_outliers)
  yavg = mean(dat_annot[[Y]])
  F1A=ggplot(dat_annot,aes_string(y=Y,x=X)) +
    ggiraph::geom_point_interactive(aes(tooltip=FUNCTION, data_id=ORF),size=2,shape=19,alpha=0.5,color='gray70',stroke=0) +
    stat_density2d(size=0.5,color='gray20') +
    geom_hline(yintercept = yavg, col='red',linetype=2,size=0.5) + # mean
    ylab('mean Evolutionary Rate (log10)') + xlab('Mean Protein Abundance (log10 ppm)') +
    ggpubr::grids()
    if(add_outliers>0){
      F1A  = F1A +
        geom_text_repel(data=OUTX, aes(label = GENENAME),max.overlaps = 20,col='blue') +
        geom_text_repel(data=OUTY, aes(label = GENENAME),max.overlaps = 20,col='red') +
        geom_point(data=OUTX, col='blue',size=0.5) +
        geom_point(data=OUTY, col='red',size=0.5)
    }
  return(F1A)
}

make_plot_1B = function(clade1='schizo',clade2='sacch.wgd',control_var='ppm',use_residuals=F,tolog=T,force_intercept=T){

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

  if(use_residuals){
    cat(sprintf("Controlling for variable '%s'\n",control_var))
    ctrl1=paste0(control_var,"1.log10")
    ctrl2=paste0(control_var,"2.log10")
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

  rho=spearman.toplot(data2plot[[xvar]],data2plot[[yvar]])

  F1B = ggplot(data=data2plot,aes_string(x=xvar,y=yvar)) +
    geom_point(size=0.7, alpha=0.8,shape=19) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour=NA,size=0.25,alpha=0.25,show.legend = F)+
    geom_smooth(method = 'lm',col='black',se=F,size=1) +
    geom_text(data=rho,aes(x=-Inf,y=-Inf,label=toshow),hjust='inward',vjust='inward',size=4) +
    xlab(labx) + ylab(laby) +
    scale_fill_distiller(palette = 'Spectral',direction = -1) +
    ggpubr::grids()
  return(F1B)
}

make_plot_1C = function(data2plot=EVOLUTION,Y='EVO.FULL',X='SNP.FULL',control_var='PPM',use_residuals=F,force_intercept=T){

  if(use_residuals){
    cat(sprintf("Controlling %s by protein expression '%s'\n",X,control_var))
    first = data2plot %>%
      broom::augment_columns(x=lm(reformulate(response=X, termlabels=control_var, intercept=force_intercept),data=.)) %>%
      dplyr::rename(.resid_x = .resid, .fitted_x=.fitted, .se.fit_x=.se.fit, .hat_x=.hat, .sigma_x=.sigma, .cooksd_x=.cooksd,.std.resid_x=.std.resid)

    cat(sprintf("Controlling %s by protein expression '%s'\n",Y,control_var))
    data2plot = first %>%
      broom::augment_columns(x=lm(reformulate(response=Y, termlabels=control_var, intercept=force_intercept),data=.)) %>%
      dplyr::rename(.resid_y = .resid, .fitted_y=.fitted, .se.fit_y=.se.fit, .hat_y=.hat, .sigma_y=.sigma, .cooksd_y=.cooksd,.std.resid_y=.std.resid)
    X='.resid_x'
    Y='.resid_y'
  }

  rho=spearman.toplot(data2plot[[X]],data2plot[[Y]])

  F1C = ggplot(data=data2plot,aes_string(x=X,y=Y)) +
    geom_point(size=0.7, alpha=0.8,shape=19) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour=NA,size=0.25,alpha=0.25,show.legend = F)+
    geom_smooth(method = 'lm',col='black',se=F,size=1) +
    geom_text(data=rho,aes(x=-Inf,y=Inf,label=toshow),hjust='inward',vjust='inward',size=4) +
    scale_fill_distiller(palette = 'Spectral',direction = -1) +
    ggpubr::grids() + theme()
  return(F1C)
}


