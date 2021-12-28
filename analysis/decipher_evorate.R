# SETUP DATA -------------------------------------------------------------------
source(here::here("analysis","function_evorate_fitting.R"))
txt_section_break = repchar("-",50)

tic("Load data")
# EXPRESSION DATA
ABUNDANCE = load.abundance()
dim(ABUNDANCE)
# EVOLUTION DATA
CLADE = load.clade()
FUNGI = load.fungi.evo()
STRAINS = load.strains.evo()
EVOLUTION = full_join(STRAINS,CLADE,by=c('ORF'='orf'))

# PROTEOME QUALITATIVE AND QUANTITATIVE VARIABLES
PROP = load.properties()
FEAT = load.features() %>% normalize_features()
PREDICTORS = full_join(PROP,FEAT)

# ANNOTATION DATA
SGD_DESC = read_rds(here('data','uniprot-sgd-annotation.rds'))
UNI_FEAT = read_rds(here('data','uniprot-features.rds'))
BIOFUNC = load.vanleeuwen2016.data()
dim(SGD_DESC)
dim(UNI_FEAT)
dim(BIOFUNC)
ANNOTATION=full_join(SGD_DESC,BIOFUNC,by="ORF") %>% full_join(UNI_FEAT,by='SGD')
toc()


# ANALYZE EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ---------------------
fit_linear_regression = function(INPUT=EVOLUTION, X='PPM', Y="log10.EVO.FULL",
                                 PREDVAR=PREDICTORS, xcor.max=0.6, ycor.max=0.6){
  txt_section_break = repchar("-",50)
  # INITIAL PARAMETERS -----------------------------------------------------------
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
  var_names = starting(colnames(PREDICTORS),'cat_')

  Y_X=make_linear_fit(LINREG ,x=X, y=Y, only.params = F)
  M0=left_join(Y_X,PREDICTORS)
  M0$yavg   = mu.y
  M0$yvar   = var.y
  M0$ymin = rg.y[1]
  M0$ymxa = rg.y[2]
  M0$nxy = n.xy

  # excluding variables with correlation to X/Y
  cor_x = cor(x=M0[,var_names],y=M0[[X]],use = 'pairwise') %>% as_vector %>% na.omit()
  x_excluded_var = rownames(cor_x)[abs(cor_x) >= xcor.max]
  n_xout = length(x_excluded_var)
  cat(sprintf("Excluding %s predictors with cor. to X > %s (%s)\n",n_xout,xcor.max,X))
  cor_y = cor(x=M0[,var_names],y=M0[[Y]],use = 'pairwise') %>% as_vector %>% na.omit()

  y_excluded_var = rownames(cor_y)[abs(cor_y) >= ycor.max]
  n_yout = length(y_excluded_var)
  cat(sprintf("Excluding %s predictors with cor. to Y > %s (%s)\n",n_yout,ycor.max,Y))

  excluded_var = unique(x_excluded_var,y_excluded_var)
  n_out = length(excluded_var)
  cat(sprintf("In total, %s predictors are excluded:\n",n_out,section))
  cat(txt_section_break,"\n ")
  cat(sprintf("%3s. %s\n",1:n_out,excluded_var))

  return(M0 %>% dplyr::select(-all_of(excluded_var)))
}

m0 = fit_linear_regression(EVOLUTION, X='PPM', Y="log10.EVO.FULL", PREDICTORS, 0.6, 0.8) %>% left_join(SGD_DESC)


### _FIGURE 1A: EVOLUTION vs EXPRESSION -------------------------------------------
make_figure_1 = function(dat=EVOLUTION, X='PPM', Y="log10.EVO.FULL",ANNOT=SGD_DESC, id=c('ORF','UNIPROT')){

  dat_annot = left_join(dat,ANNOT,by=id)
  OUTLIERS.Y = get_outliers(dat_annot,X)
  OUTLIERS.X = get_outliers(dat_annot,Y)
  yavg = mean(dat_annot[[Y]])
  F1A=ggplot(dat_annot,aes_string(y=Y,x=X)) +
    ggiraph::geom_point_interactive(aes(tooltip=FUNCTION, data_id=ORF),size=2.5,shape=19,alpha=0.5,color='gray70',stroke=0) +
    stat_density2d(size=0.5,color='gray20') +
    geom_hline(yintercept = yavg, col='red',linetype=2,size=0.5) + # mean
    ylab('mean Evolutionary Rate (log10)') + xlab('Mean Protein Abundance (log10 ppm)') +
    ggpubr::grids() +
    geom_text_repel(data=OUTLIERS.Y, aes(label = GENENAME),max.overlaps = 20,col='blue') +
    geom_text_repel(data=OUTLIERS.X, aes(label = GENENAME),max.overlaps = 20,col='red')

  return(F1A)
}
make_figure_1(dat=EVOLUTION,X=PPM,Y=log10.EVO.FULL)
ggiraph::girafe(ggobj = F1A)

### _FIGURE 1B: BRANCH LENGTH vs EXPRESSION ---------------------------------------
F1B

### _FIGURE 1C: SNP EVOLUTION vs EXPRESSION ---------------------------------------
plot(y=TMP$Kc2.log10, x=TMP$ppm2.log10)

TMP = left_join(CLADE,INPUT, by=c('uni2'='UNIPROT.x'))
m=lm(data=TMP,Kc2.log10~ppm2.log10)
decompose_variance(m)
dim(TMP)

spearman(TMP[[Y]], TMP$MPC)

INPUT %>%
slice.iqr
pXY.fit = pXY+
  geom_line(M,mapping=aes(y=.fitted,col=model),size=1,show.legend = F) +
  geom_text(M.x1 ,mapping=aes(label=model,col=model,y=.fitted-0.2),x=0.8,size=5,check_overlap = T,show.legend = F) +
  stat_smooth(method = 'loess', fullrange=T, span=1.2, se = F, col='gray40',size=1) + annotate('text',y=1.9,x=0.8,label='Loess',col='gray40',size=5)
#coord_cartesian(xlim = c(-2.5, 4.5), ylim = c(0, 3.5),expand = F) +
pXY.fit


# Filter variables (abs cor MPC<0.4)
# Make M0 and M3
# Refine variables selections





