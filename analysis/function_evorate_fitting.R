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

load.clade = function(){
  # BRANCH LENGTH IN FUNGI CLADES
  # Normalized branch length (Kc):
  #  Kc = SUM(branch lengths in clade subtree Tc) / SUM(branch lengths in species tree Ts)
  enog.cols = c('Tax','NOG','nog.1to1','RAXML','FunCat','STRING','sp1','sp2','uni1','uni2','ppm1','ppm2')
  clade.cols = c('ppm1.log10','ppm2.log10','clade1','clade2','XX','YY')
  clade = get_clade_data(g1='schizo',g2='sacch.wgd',rate = 'ratio') %>%
    dplyr::select(all_of(enog.cols), all_of(clade.cols), starts_with('schizo'), starts_with('sacch.wgd')) %>%
    dplyr::rename(Kc1.log10=XX,Kc2.log10=YY)
  return(clade)
}

load.fungi.evo = function(){
  ### FUNGI CONTAIN ALL DATA ABOUT THE FUNGI LINEAGE EVOLUTIONARY RATE
  fungi = load.dubreuil2021.data(1) %>%
    dplyr::select('ORF','UNIPROT','PPM',c(starts_with('EVO.'))) %>%
    dplyr::mutate(across(starts_with("EVO."),function(x){ x/mean_(x) },.names="norm.{.col}")) %>% # scale and center all evo rate
    dplyr::mutate(across(starts_with("norm.EVO."),log10,.names = "log10.{.col}")) # # apply log10 to rate4site
  return(fungi)
}

load.strains.evo = function(){
  ### STRAINS CONTAIN ALL DATA ABOUT THE FUNGI LINEAGE EVOLUTIONARY RATE AND THE S.CEREVISIAE POPULATION
  strains = readRDS(here("data","PROTEIN-EVO-FUNGI-SNP.rds")) %>%
    dplyr::select(c(starts_with(c('EVO.','SNP.')),'PPM','ORF','UNIPROT','IS_FUNGI','IS_STRAINS')) %>%
    filter(!(is.dup(UNIPROT) & is.na(EVO.FULL))) %>%
    group_by(ORF) %>% dplyr::mutate(across(starts_with(c("EVO.","SNP.")),log10,.names = "log10.{.col}")) %>% # apply log10 to SNP rate4site
    group_by(ORF) %>% dplyr::mutate(across(starts_with(c("EVO.","SNP.")),function(x){ x/mean_(x) }))  # scale and center all EVO/SNP rate
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
  feat = readRDS(get.last.file(here("output"),"proteome-features")) %>% ungroup()

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
  col_f_aa = grep("cat_biophysics.uniprot.f_",colnames(feat))
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

# 2. LINEAR FIT ----------------------------------------------------------------

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
