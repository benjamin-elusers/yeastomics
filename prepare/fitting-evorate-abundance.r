#### SETUP ####
#### _Functions ####
yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/"
source(file.path(yeastomics,"src/utils.r"))
source(file.path(yeastomics,"src/function_annotation.r"))
source(file.path(yeastomics,"src/function_sequence.r"))
source(file.path(yeastomics,"src/function_phylogenetic.r"))
source(file.path(yeastomics,"src/function_analysis.r"))
source(file.path(yeastomics,"src/function_datalocal.r"))
source(file.path(yeastomics,"src/function_datapub.r"))
library(tidyverse)
library(tictoc)
library(hablar)
#### _Graphics ####
library(cowplot)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(ggrepel)
library(ggpubr)
library(see)
mytheme =   theme_clean(base_size = 14) +
  theme(line=element_line(size=1),
        panel.grid.major = element_line(size=1,lineend='round',linetype='22',color='gray50'),
        panel.grid.minor = element_line(size=0.5,lineend='round',linetype='12',color='gray50'),
        axis.ticks.length = unit(1.5,'mm'),
        axis.ticks = element_line(size=0.5),
        legend.position = 'none')
source("https://raw.githubusercontent.com/clauswilke/dviz.supp/master/R/dviz.supp.R")
source("https://raw.githubusercontent.com/clauswilke/dviz.supp/master/R/themes.R")
empty_theme <- theme_dviz_open(12, rel_small = 1, rel_large = 1) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = grid::unit(0, "pt")
  )
theme_set( theme_modern(base_size = 16, axis.text.angle = 45,legend.position = 'none') + theme(aspect.ratio=1) )

#### _Analysis ####
show_density = function(input, # Input data
                        var,   # Variable name
                        nsd=2  # Number of standard deviation to show
){

  V = input[[var]]

  # Statistics
  mu = mean_(V) # Y mean
  md = median_(V) # Y median
  s  = sd_(V) # Y standard deviation
  rg = range_(V)
  symmetry = function(x){ return( sort(c(x,x)) * (-1)^(1:length(x)) ) }
  sd_x = mu + s* symmetry(x=1:nsd)
  # Density
  D=density(V,bw=0.1,na.rm=T)
  dy.ind = sapply(sd_x,function(X){ nearest(X,D$x) })
  sd_y  = D$y[dy.ind]

  df_sd =tibble( xx = sd_x, yy = sd_y)

  A = ggplot(input) +
    ggpubr::grids() +
    stat_density(bw = 0.1, aes_string(x=V)) +
    geom_vline(xintercept = mu, col='white', linetype=1,size=1) + # mean
    geom_segment(data=df_sd,aes(x=xx,xend=xx,y=rep(0,nsd*2),yend=yy), col='white', linetype=2,size=1) + # median
    geom_segment(data=df_sd[1:2,],aes(x=xx,xend=xx,y=yy,yend=c(1,1)), col='red', linetype=2,size=0.25) + # standard deviation (outside)
    geom_errorbarh(data=df_sd,aes(xmin=xx[1],xmax=xx[2],y=0.8,height=0.05),col='red',linetype=1,size=1) + # standard deviation
    ylim(0,1) + scale_y_continuous(position = "left") + xlab("Mean Evolutionary Rate") +
    theme( axis.line.x = element_blank(),axis.ticks = element_blank()) + #aspect.ratio = 3/2,
    coord_flip(xlim=rg, ylim=c(1,0))
  return(A)
}

show_sample = function(input, # Input data
                       pop.mean,pop.range,
                       name='property', # column for sample name
                       id='ORF', # column for id of single observation
                       value=c('MF_nucleotide_binding','MF_molecular_function'), # sample names
                       var){

  which.sample = input[[name]] %in% value
  selected = input %>% dplyr::filter(which.sample) %>% dplyr::select(!!sym(name), !!sym(id), !!sym(var))
  mu.sample = mean_(selected[[var]])

  df.mu = group_by(selected,!!sym(name)) %>%
    summarise(n=n(),mu.sample=mean_(!!sym(var)),sd.sample=sd_(!!sym(var)),
              sdmin=mu.sample-sd.sample, sdmax=mu.sample+sd.sample) %>%
    mutate(MU = pop.mean)


  #selected.prop = c('MF_nucleotide_binding','MF_molecular_function','essential_core','essential_dispensable','pangenome_rare','pangenome_cloud')
  #ER.prop = propfit %>% dplyr::filter(property %in% selected.prop) %>% dplyr::select('property','EVO.FULL','.resid.evo')
  library(ggbeeswarm)
  library(see)

  B = ggplot(selected,aes(y=!!sym(var),fill=!!sym(name),col=!!sym(name))) +
    geom_violinhalf(aes_string(x=name),color=NA,alpha=0.9,show.legend = F) +
    geom_beeswarm(aes_string(x=name),size = 2,  shape = 21, stroke = 0, alpha=0.5,groupOnX = T,dodge.width=0) +
    geom_pointrange_borderless(data=df.mu,mapping = aes(x=!!sym(name),y=mu.sample,ymin=sdmin,ymax=sdmax),size=1,fatten=6, position=position_dodge2(width=1)) +
    #geom_crossbar(data=df.mu,mapping = aes(x=!!sym(name),y=mu.sample,ymin=sdmin,ymax=sdmax),alpha=0.1,size=0.5,fatten=0) +
    geom_hline(df.mu,mapping = aes(yintercept = MU), col='black',linetype=1,size=1) + # mean
    geom_text(data=df.mu,mapping = aes(label=n,x=!!sym(name),y=2.5,vjust='inward'),show.legend = F) +
    xlab(name) + scale_y_continuous(name='',limits = pop.range, expand = c(0,0)) + ylab('') +
    #theme(legend.position = 'none',legend.direction = 'vertical',axis.text.x = element_text(angle = 45)) +
    ggpubr::grids()+
    #scale_fill_material_d(palette = "ice") + scale_color_material_d(palette = "ice")
    scale_fill_material_d(palette="full") + scale_color_material_d(palette="full")
  B
  return(B)
}

get_XY_data = function(input,x=X,y=Y,noNA=T){
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

make_poly_fit = function(input,    # Input data
                         y=Y, x=X, # X/Y Variables to fit
                         deg=3,
                         only.params=T){
  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']

  # model #
  model.name=sprintf('Polynomial\n(d = %s)',deg)
  f=as.formula(paste0(y,"~poly(",x,",degree=",deg,",raw=T)"))
  m = glm(formula =f , data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model = model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

make_expo_fit = function(input,    # Input data
                         y=Y, x=X, # X/Y Variables to fit
                         only.params=T){
  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']

  expo = function(alpha=1,beta=0,x){ return(alpha * exp(beta*x)) }
  # model #
  model.name='Exponential'
  f=as.formula(paste0(y,"~ expo(alpha,beta,",x,")"))
  init.params=list(alpha=1,beta=-1)
  m=nls(f,start = init.params,data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model = model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

decompose_variance = function(LM){
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
#### LOADING EVOLUTIONARY DATA (FUNGI & STRAINS) ####
tic("Load data")
CLADE = get_clade_data(g1='schizo',g2='sacch.wgd',rate = 'ratio')
FUNGI = readRDS(here("data","fungi-evodata.rds")) %>% filter( !is.na(PPM) & !is.na(MPC) )
STRAINS = readRDS(here("data","yeast_strains-evodata.rds")) %>%
          filter( !is.na(PPM) & !is.na(MPC) ) %>% distinct %>%
          filter(!(is.dup(UNIPROT) & is.na(EVO.FULL)))
PROP=readRDS(get.last.file(here("output"),"proteome-properties")) %>%
  dplyr::select(-contains("chemsig_X")) %>% ungroup()

prop = PROP %>%
  pivot_longer(cols = starts_with('cat_'),
               names_to = c('categories','source','property'),
               names_pattern="cat_(.+)\\.(.+)\\.(.+)",
               values_to = "has_prop") %>%
  mutate(col_prop = paste0(categories,'.',source,'.',property)) %>%
  dplyr::filter(has_prop)

FEAT = readRDS(get.last.file(here("output"),"proteome-features")) %>% ungroup()
features = FEAT %>%
  pivot_longer(cols = starts_with('cat_'),
               names_to = c('categories','source','feature'),
               names_pattern="cat_(.+)\\.(.+)\\.(.+)",
               values_to = "value") %>%
  mutate(col_feat = paste0(categories,'.',source,'.',feature)) %>%
  dplyr::filter(!is.na(value)) %>% add_count(col_feat,name='size')
toc()

FUNGI.NORM = FUNGI %>% mutate(rel_EVO.FULL=log2(EVO.FULL/median_(EVO.FULL)),
                              rel_EVO.DISORDER=log2(EVO.DISORDER/median_(EVO.DISORDER)),
                              rel_EVO.NOT_DISORDER=log2(EVO.NOT_DISORDER/median_(EVO.NOT_DISORDER)),
                              rel_EVO.DOMAINS=log2(EVO.DOMAINS/median_(EVO.DOMAINS)),
                              rel_EVO.ANTI_DOMAIN=log2(EVO.ANTI_DOMAIN/median_(EVO.ANTI_DOMAIN)),
                              rel_EVO.PDB=log2(EVO.PDB/median_(EVO.PDB)),
                              rel_EVO.SURFACE=log2(EVO.SURFACE/median_(EVO.SURFACE)),
                              rel_EVO.BURIED=log2(EVO.BURIED/median_(EVO.BURIED))
                              )


INPUT = FUNGI.NORM
evocols = INPUT %>% dplyr::select(starts_with('rel_EVO.')) %>% colnames
#EVO = left_join(INPUT,PROP) %>% ungroup()

#### EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ####
Y = "rel_EVO.FULL" # mean Evolutionary rate (full sequence)
#X = "PPM" # Protein Abundance (log10 ppm)
X = "MPC" # median Molecules Per Cell
XYDATA = get_XY_data(INPUT,x=X,y=Y)
YY = XYDATA$YY
XX = XYDATA$XX
n.xy = XYDATA$n['xy']
mu.y = XYDATA$mu['y']
var.y = XYDATA$var['y']
input = XYDATA$df
rg.y = range_(YY)
dim(input)
#INPUT[[X]][ INPUT[[X]] < 0 ] = 0
#### _Distribution ####

samples = c('MF_nucleotide_binding','MF_molecular_function')#,'essential_core','essential_dispensable')
input.byprop = left_join(INPUT,prop) %>% ungroup()
A=show_density(input,var = Y) + xlim(rg.y) +scale_y_continuous(name = 'density', expand = c(0,0))
B=show_sample(input=input.byprop ,name='property',id='ORF',value=samples,var=Y, pop.mean=mu.y, pop.range=rg.y) + ylim(rg.y) + theme(aspect.ratio=1)

ggsave(plot=A, filename=here("output","F4.2A_evo_density.pdf"),useDingbats=F)
ggsave(plot=B, filename=here("output","F4.2B_samples_evo_density.pdf"),useDingbats=F)
f4.2A = plot_grid(A,B, align = 'hv',nrow=1)
ggsave(plot=f4.2A, filename = here('output','f4.2A-evofit.pdf'), scale = 1, useDingbats = F)

#### _Fiting model ####
M1=make_linear_fit(input,only.params = F)
M2=make_logistic_fit(input,only.params = F)
M3=make_poly_fit(input,only.params = F,deg=3)
M4=make_expo_fit(input,only.params = F)
M = bind_rows(M1,M2,M3,M4)

selected.par =c('model','RSS','ESS','TSS')
fit.params = bind_rows(
    make_linear_fit(input,only.params = T)[selected.par],
    #make_logistic_fit(input,only.params =T)[selected.par],
    make_poly_fit(input,only.params = T,deg=3)[selected.par],
    make_expo_fit(input,only.params = T)[selected.par]
  )
fit.params


#### _ EVO.vs.MPC SCATTERPLOT WITH FIT   ####
x1=nearest(0.5,XX,n=1,v=T)
M.x1 = M[ M[X] == x1 ,]
spearman(input[[X]],input[[Y]])

pXY=ggplot(input,aes_string(y=Y,x=X)) +
  geom_point(size=2.5,shape=19,alpha=0.5,color='gray70',stroke=0) +
  stat_density2d(size=0.5,color='gray20') +
  #  geom_text(data=CC.ER, aes(label=sprintf(" r %.2f \n p %.1e \n N %s ",estimate,p.value,n)),
  #            x=Inf,y=Inf,hjust='inward',vjust='inward',size=3) +
  ylab('') + xlab('Median Protein Abundance (log10 mpc)') +
  geom_hline(yintercept = mu.y, col='black',linetype=1,size=1) + # mean
  ylim(rg.y) + theme(legend.position = 'none',aspect.ratio = 1) +
  ggpubr::grids()


M.x1$model %in% c('linear')
pXY.fit = pXY+
  geom_line(M,mapping=aes(y=.fitted,col=model),size=1,show.legend = F) +
  geom_text(M1 ,mapping=aes(label=model,col=model,y=.fitted-0.2),x=0.8,size=5,check_overlap = T,show.legend = F) +
  stat_smooth(method = 'loess', fullrange=T, span=1.2, se = F, col='gray40',size=1) + annotate('text',y=1.9,x=0.8,label='Loess',col='gray40',size=5)
  #coord_cartesian(xlim = c(-2.5, 4.5), ylim = c(0, 3.5),expand = F) +
pXY.fit
ggsave(plot=pXY, filename=here("output","F4.2C_evo_mpc_fit.pdf"))

MF_nuc_bind = PROP$ORF[PROP$cat_functions.go.MF_nucleotide_binding]
MF_unknown  = PROP$ORF[PROP$cat_functions.go.MF_molecular_function]
linfit = geom_line(M1,mapping=aes(y=.fitted,col=model),size=1,show.legend = F)
PXY=plot_grid(pXY.fit,
          pXY + geom_point(data=input %>% filter(ORF %in% MF_unknown), color='green',size=1)+linfit,
          pXY + geom_point(data=input %>% filter(ORF %in% MF_nuc_bind) , color='blue',size=1) +linfit,
          nrow=3) + theme_clean()

#### _ RESIDUAL EVO.vs.MPC SCATTERPLOT (EVO.FULL) ####
col.samples = PROP %>% dplyr::select(contains(samples)) %>% colnames
library(ggiraph)
library(ggiraphExtra)
library(gghighlight)
fitdata = left_join(M1,PROP, by='ORF')# %>%
          #bind_rows(M1 %>% mutate(property='all_orthologs', has_prop=T) ) %>%
          #dplyr::filter( property %in% c(samples,'all_orthologs')) %>% mutate(property = as.factor(property))
#res.ortho = res.data %>% filter(property=='all_orthologs')
f0=lm(data=fitdata, formula = EVO.FULL ~MPC)
f1=lm(data=fitdata, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) + cat_functions.go.MF_nucleotide_binding + cat_functions.go.MF_molecular_function)

fa=lm(data=fitdata, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) +  cat_functions.go.MF_molecular_function)
fb=lm(data=fitdata, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) + cat_functions.go.MF_nucleotide_binding)


fa=lm(data=fitdata[ fitdata$cat_functions.go.MF_molecular_function,], formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) +  cat_functions.go.MF_molecular_function)
f0a=lm(data=fitdata[ fitdata$cat_functions.go.MF_molecular_function,], formula = EVO.FULL ~ offset(coef(f0)[2]*MPC))

var(fa$model$EVO.FULL)*nrow(fa$model)
var(f0a$model$EVO.FULL)*nrow(f0a$model)

deviance(fa)
deviance(f0a)


fitdata$fitted.intercept = predict(f1)
fitdata$resid.intercept = residuals(f1)
fitdata.toplot = pivot_longer(fitdata,
                              cols = starts_with('cat_'),
                              names_to = c('categories','source','property'),
                              names_pattern="cat_(.+)\\.(.+)\\.(.+)",values_to = "value")


fitdata.toplot.samples = fitdata.toplot %>% filter(property %in% samples & value==T)

yfits=ggplot(fitdata.toplot.samples ,aes_string(y=Y,x=X,col='property',fill='property')) +
  stat_density2d(size=0.5) +
  #geom_point(size=1.5,shape=19,stroke=0,col='gray70',) +
  geom_point(size=1.5,shape=19,stroke=0,alpha=0.5) +
  geom_line(aes(y=.fitted,col=property),col='gray40',size=1) +
  geom_line(aes(y=fitted.intercept,col=property),size=1) +
  ylab('Mean Evolutionary Rate') + xlab('') +#xlab('Median Protein Abundance (log10 mpc)') +
  ylim(rg.y) +  facet_wrap(~property,nrow = 1) +
  #geom_smooth(data=res.prop,method='lm', formula = str(fit$call), color='red') +
  gghighlight(property %in% c('all_orthologs',samples),use_group_by = F, use_direct_label = F) +
  theme(legend.position = 'none',aspect.ratio = 1) +
  ggpubr::grids() +
  geom_hline(yintercept = mu.y, col='black',linetype=1,size=1) +  # mean
  scale_fill_material_d(palette="full") + scale_color_material_d(palette="full")
  #  geom_text(data=CC.ER, aes(label=sprintf(" r %.2f \n p %.1e \n N %s ",estimate,p.value,n)),
  #            x=Inf,y=Inf,hjust='inward',vjust='inward',size=3) +
yfits

yres1= pXY+
  geom_point(color='gray70',shape=19,alpha=0.3) +
  geom_line(data=fitdata.toplot.samples,aes_string(y='.fitted',col='property',fill='property'),linetype=1,size=1) +
  geom_point(data=fitdata.toplot.samples,aes_string(col='property',fill='property'),shape=19,alpha=0.7,size=1.5) +
  geom_segment(data=fitdata.toplot.samples,aes_string(xend = X, yend = '.fitted',col='property',fill='property'),linetype=1,size=0.75) +
  facet_wrap(~property) +
  geom_hline(yintercept = mu.y,size=1,linetype=2)+
  scale_fill_material_d(palette="full") + scale_color_material_d(palette="full")
yres1

rg.resy = range_(M1$.resid)
pXYres=ggplot(M1,aes_string(y='.resid',x=X)) +
  geom_point(size=2.5,shape=19,alpha=0.5,color='gray70',stroke=0) +
  ylab('residual mean Evolutionary Rate') + xlab('Median Protein Abundance (log10 mpc)') +
  geom_hline(yintercept = 0, col='black',linetype=1,size=1) + # mean
  ylim(rg.resy) + theme(legend.position = 'none',aspect.ratio = 1) +
  ggpubr::grids()

yres2 = pXYres +
   geom_line(data=fitdata.toplot.samples,aes_string(y=0,col='property',fill='property'),linetype=1,size=0.05) +
   geom_point(data=fitdata.toplot.samples,aes_string(col='property',fill='property'),shape=19,alpha=0.5,size=2) +
   geom_segment(data=fitdata.toplot.samples,aes_string(xend = X,y=0,yend='.resid',col='property',fill='property'),linetype=1,size=0.75) +
   facet_wrap(~property) +
   scale_fill_material_d(palette="full") + scale_color_material_d(palette="full")
yres2


#### _RESIDUAL SNP.vs.MPC SCATTERPLOT WITH RESIDUAL (SNP.FULL) ####
cat("==> Get residuals of long/short term evolutionary rate from abundance <==\n")

r4s.out = c("YBL003C","YBR009C","YBR082C","YDL137W","YDL192W","YDR039C","YDR225W","YDR385W","YEL026W","YFL039C",
            "YGL189C","YIL173W","YJR145C","YLR164W","YLR293C","YLR340W","YNL030W","YNL031C","YOR185C","YPR016C")

res.evo = STRAINS %>%
    broom::augment_columns(x=lm(SNP.FULL~MPC,data=.)) %>%
    dplyr::rename(snpfull.resid = .resid, snpfull.fitted=.fitted, snpfull.se=.se.fit, snpfull.hat=.hat, snpfull.sigma=.sigma, snpfull.cooksd=.cooksd,snpfull.std.resid=.std.resid)  %>%
    broom::augment_columns(x=lm(EVO.FULL~MPC,data=.)) %>%
    dplyr::rename(evofull.resid = .resid, evofull.fitted=.fitted, evofull.se.fit=.se.fit, evofull.hat=.hat, evofull.sigma=.sigma, evofull.cooksd=.cooksd,evofull.std.resid=.std.resid) %>%
    mutate(outlier.snp = snpfull.cooksd > (20 * mean(snpfull.cooksd)) )

FXA = make_scatterplot(res.evo,xvar='SNP.FULL',yvar='EVO.FULL',
                       labx='SNP Evolutionary Rate',laby='mean Evolutionary Rate',
                       theme2use = theme(),txtcorner = 'topleft') + geom_point(data=subset(res.evo,outlier.snp), aes(text=paste0(ORF)),col='red') +
      theme(legend.position='none',aspect.ratio = 1)
plotly::ggplotly(FXA)

FXB = make_scatterplot(res.evo,xvar='snpfull.resid',yvar='evofull.resid',
                       labx='residual SNP Evolutionary Rate',laby='residual Evolutionary Rate',
                       theme2use = theme(),txtcorner = 'bottomright') +  theme(legend.position='none',aspect.ratio = 1)

res.evo.clade = res.evo %>% filter(ORF %in% CLADE$orf)

cor(res.evo$snpfull.resid,res.evo$evofull.resid)^2
FXB.S1 = make_scatterplot(res.evo,xvar='SNP.FULL',yvar='evofull.resid',
                       labx='SNP Evolutionary Rate',laby='residual Evolutionary Rate',
                       theme2use = theme(),txtcorner = 'bottomright') +  theme(legend.position='none',aspect.ratio = 1)


FXC = make_scatterplot(res.evo.clade,xvar='SNP.FULL',yvar='EVO.FULL',
                       labx='SNP Evolutionary Rate',laby='mean Evolutionary Rate',
                       theme2use = theme(),txtcorner = 'topleft') +  theme(legend.position='none',aspect.ratio = 1)
FXD = make_scatterplot(res.evo.clade,xvar='snpfull.resid',yvar='evofull.resid',
                       labx='residual SNP Evolutionary Rate',laby='residual Evolutionary Rate',
                       theme2use = theme(),txtcorner = 'bottomright') +  theme(legend.position='none',aspect.ratio = 1)
plot_grid(FXA,FXB,FXC,FXD)
save_plot(plot_grid(FXA,FXB),filename = here('output','Fig3.3-evo.vs.snp-with-residuals.pdf'))

####
fig4.2B = plot_grid(PXY,resi,align = 'hv')
ggsave(plot=fig4.2B, filename = here('output','f4.2B-evofit.pdf'), scale = 1,useDingbats=F)
#### FIGURE ####
library(cowplot)
fig4.2.top = plot_grid(A,B,pXY.fit,nrow=1,ncol=3,align='hv', axis = 'tb')
fig4.2.bottom = plot_grid(plot_grid(yres1,yres2,nrow=2,align = 'hv'),yfits,align = 'hv',ncol=2)
ggsave(fig4.2.bottom, filename = here('output','f4.2C-residuals.pdf'), scale = 2,useDingbats=F)

fig4.2=plot_grid(fig4.2.top,resi,nrow=2)
ggsave(fig4.2,filename = here('output','f4.2-evofit.pdf'),scale=2,useDingbats=F)



### REPRESENTING PROPERTIES AND FEATURES DATA MATRIX
bigprop  = left_join(INPUT,PROP) %>%
            ungroup() %>%  dplyr::select(where( ~ is_logical(.x) && sum(.x) >150)) %>%
            rowwise %>% filter( sum(c_across(everything()))> 5) %>%
            as.matrix
rankcols = order( colSums(bigprop), decreasing = T)
rankrows = order( rowSums(bigprop), decreasing = T)
nc = ncol(bigprop)
nr = nrow(bigprop)
yeastQR = bigprop[rankrows,rankcols][sample(1:300,50,replace=F),sample(1:nc,50,replace=F)]
pheatmap::pheatmap(mat=1-as.matrix(1*yeastQR),
                   color = c('white','black'),
                   show_rownames = F, show_colnames = F,
                   cluster_rows = F, cluster_cols =T,legend = F,clustering_distance_rows = 'binary',
                   clustering_distance_cols = 'binary',treeheight_col=0, treeheight_row=0,
                   cutree_rows = 20, cutree_cols = 1,
                   cellwidth = 12, cellheight = 12, border=NA,
                   width=10,height=10, filename = here('output','fig3.6-yeast-qr-code-properties.pdf'))

bigfeat = left_join(INPUT,features.norm) %>%  dplyr::select(where(~ is.numeric(.x) && max_(.x)< 1e4)) %>% as.matrix()
nr=nrow(bigfeat)
nc=ncol(bigfeat)
pheatmap::pheatmap(mat=bigfeat[sample(1:nr,50),sample(1:nc,50)], scale='column',
                   show_rownames = F, show_colnames = F,
                   cluster_rows = F, cluster_cols =T,legend = F,
                   clustering_distance_rows = 'binary',
                   clustering_distance_cols = 'binary',
                   treeheight_col=0, treeheight_row=0,
                   cutree_rows = 20, cutree_cols = 1,
                   cellwidth = 12, cellheight = 12, border=NA,
                   width=10,height=10, filename = here('output','fig3.7-yeast-features-heatmap.pdf'))



#### FINAL FIT AND RESIDUALS####
colab = c('PPM','MPC','gfp','ms')

#### FEATURES NORMALIZATION ####
features = FEAT
features.norm = features

##### a. Normalize codons counts (0-100) => 64 COLUMNS #####
col_codons = grep("cat_transcriptomics.sgd.[ATCG]{3}",colnames(features))
features.norm[,col_codons] =  100 * features[,col_codons] / (features$cat_transcriptomics.sgd.prot_size+1)
##### b. Normalize amino acid frequencies (0-100) => 31 COLUMNS #####
col_f_aa = grep("cat_biophysics.uniprot.f_",colnames(features))
features.norm[,col_f_aa] =  100 * features[,col_f_aa]
##### c. Normalize all other fraction (0-100) => 4 COLUMNS #####
col_frac = c("cat_genomics.sgd.pGC","cat_genomics.byrne2005.RLEN","cat_transcriptomics.coRdon.CU_fop",
  "cat_biophysics.d2p2.f")
#"cat_biophysics.dubreuil2019.IUP20_f","cat_biophysics.dubreuil2019.IUP30_f","cat_biophysics.dubreuil2019.IUP40_f")
features.norm[,col_frac] =  100 * features[,col_frac]
##### d. Normalize count/doubling variables => 30 COLUMNS #####
col_count = c("cat_transcriptomics.paxdb.ppm_4932", "cat_transcriptomics.paxdb.ppm_214684","cat_transcriptomics.paxdb.ppm_4896","cat_transcriptomics.paxdb.ppm_5061",
              "cat_transcriptomics.paxdb.ortho_ppm_avg","cat_transcriptomics.paxdb.ortho_ppm_sd",
              "cat_transcriptomics.paxdb.ortho_ppm_max","cat_transcriptomics.paxdb.ortho_ppm_min",

              "cat_interactions.string.cent_pagerank",
              "cat_interactions.string.cent_eigen","cat_interactions.string.cent_authority",
              "cat_interactions.string.cent_hub","cat_interactions.string.cent_subgraph","cat_interactions.intact.cent_deg",
              "cat_interactions.intact.cent_pagerank","cat_interactions.intact.cent_eigen",
              "cat_interactions.intact.cent_authority",
              "cat_interactions.intact.cent_hub","cat_interactions.intact.cent_subgraph")
features.norm[,col_count] =  log10(features[,col_count])

col_doublings = c("cat_genomics.sgd.len","cat_transcriptomics.sgd.prot_size",
                  "cat_genomics.byrne2005.S","cat_genomics.byrne2005.N","cat_genomics.byrne2005.G",
                  "cat_transcriptomics.geisberg2014.HL_mrna", "cat_biophysics.villen2017.HL_prot",
                "cat_biophysics.d2p2.nseg","cat_biophysics.d2p2.L","cat_biophysics.d2p2.Lsegmax",
                "cat_biophysics.dubreuil2019.IUP20_L","cat_biophysics.dubreuil2019.IUP30_L","cat_biophysics.dubreuil2019.IUP40_L")
features.norm[,col_doublings] =  log2(features[,col_doublings]+1)


#### GET FULL DATASET WITH PROPERTIES AND FEATURES TOGETHER ####
load(here::here('output','checkpoint-fitting-evorate.rdata'))

uniq_orf = unique(fitdata$ORF)
PROP_FEAT = left_join(PROP,features.norm)
YEASTOMICS = left_join(M1,PROP_FEAT, by='ORF') %>% dplyr::filter(ORF %in% uniq_orf & !duplicated(ORF))

YEASTOMICS = left_join(STRAINS,PROP_FEAT, by='ORF') %>% dplyr::filter(!duplicated(ORF))
dim(YEASTOMICS)
saveRDS(YEASTOMICS, file=here::here('output','yeastOmics-290921.rds'))

unknown = c("cat_functions.go.BP_biological_process","cat_functions.go.MF_molecular_function")
GOBP = YEASTOMICS %>% dplyr::select(contains("go.BP"))
table(GOBP$cat_functions.go.BP_biological_process)
GOBP$cat_functions.go.BP_biological_process[ rowSums(GOBP) > 1 & GOBP$cat_functions.go.BP_biological_process ] = F
YEASTOMICS[,colnames(GOBP)] = GOBP
table(YEASTOMICS$cat_functions.go.BP_biological_process)

GOMF = YEASTOMICS %>% dplyr::select(contains("go.MF"))
table(GOMF$cat_functions.go.MF_molecular_function)
GOMF$cat_functions.go.MF_molecular_function[ rowSums(GOMF) > 1 & GOMF$cat_functions.go.MF_molecular_function] = F
YEASTOMICS[,colnames(GOMF)] = GOMF
table(YEASTOMICS$cat_functions.go.MF_molecular_function)

dim(YEASTOMICS)

ALL_PROP =  YEASTOMICS %>%
       dplyr::select( colnames(PROP), colnames(M1) ) %>%
       pivot_longer(cols = starts_with('cat_'),
               names_to = c('categories','source','property'),
               names_pattern="cat_(.+)\\.(.+)\\.(.+)",
               values_to = "has_prop") %>%
       mutate(col_prop = paste0(categories,'.',source,'.',property)) %>%
       dplyr::filter(has_prop) %>%
       group_by(col_prop) %>%
       mutate( N=n(),
              .fitted.prop = .fitted + mean(.resid),
              .resid.prop = .resid + mean(.resid),
              RS_avg = mean_(.resid),
              EXPECTED = sign(RS_avg),
              TSS.p  = sum( (rel_EVO.FULL - mu.y)^2 ),
              ESS.p  = sum( (.fitted.prop-mu.y)^2),
              RSS.p = TSS.p - ESS.p,

              tss.pc = 100*TSS.p/TSS,
              ess = sum( (.fitted-mu.y)^2),
              rss = TSS.p - ess,
              ess.pc  = (100*ESS.p / TSS.p) * sign(RS_avg),
              rss.pc  = (100*RSS.p / TSS.p) * sign(RS_avg),
              dRSS  = RSS.p-rss,
              dRSS.rss  = 100*dRSS / RSS,
              dRSS.tss  = 100*dRSS / TSS,
              dESS  = ESS.p-ess,
              dESS.ess = 100*dESS / ESS,
              dESS.tss = 100*dESS / TSS,
  ) %>%
  dplyr::select(categories,source,property, col_prop,
                N, RS_avg, EXPECTED, TSS,ESS,RSS,
                TSS.p, tss.pc,
                ESS.p, ess , ess.pc,
                RSS.p, rss, rss.pc,
                dRSS.rss, dRSS.tss,
                dRSS, dESS, dESS.ess, dESS.tss ) %>%
  distinct()

BEST_PROP =  ALL_PROP %>% mutate(tss.1pc = tss.pc > 0.01*TSS) %>%
  dplyr::filter( abs(dRSS.tss)>1 & abs(dESS.tss)>1)

ALL_FEAT = YEASTOMICS %>%
  dplyr::select( colnames(features.norm), colnames(M1) ) %>%
  pivot_longer(cols = starts_with('cat_'),
               names_to = c('categories','source','feature'),
               names_pattern="cat_(.+)\\.(.+)\\.(.+)",
               values_to = "value") %>%
  mutate(col_feat = paste0(categories,'.',source,'.',feature)) %>%
  dplyr::filter(!is.na(value)) %>% add_count(col_feat,name='size') %>%
  mutate( is.codon  = feature %in% get.codons4tai(), is.aa = feature %in% paste0('f_',get.AA1())) %>%
  rowwise %>% mutate( feature=ifelse(is.codon ,sprintf('%s (%s)',tolower(feature), get.AA3()[Biostrings::GENETIC_CODE[feature]]),feature)) %>%
  rowwise %>% mutate( feature=ifelse(is.aa ,sprintf('f(%s)',get.AA3()[gsub("f_","",feature)]),feature) ) %>%
  group_by(col_feat,categories,source,feature) %>%
  summarise( r=scor(.resid,value,met='pearson')$estimate,
             rho=scor(.resid,value)$estimate,
             nfeat = sum.na(value,notNA=T) ) %>%
  ungroup() %>% mutate( rk = rank(-rho,ties.method = 'first') ) %>%
  arrange(rk) %>% relocate(rk)

BEST_FEAT = ALL_FEAT %>%
            filter(abs(r) > 0.2  & !grepl("byrne2005",col_feat) & !grepl("snp_",col_feat) ) %>%
            mutate(SIGN=sign(r))

manual.selection = c("cat_biophysics.uniprot.f_polar_DEHKNQRSTZ",
                     "cat_biophysics.uniprot.f_turnlike_ACDEGHKNQRST",
                     "cat_biophysics.uniprot.f_alcohol_ST",
                     "cat_biophysics.uniprot.f_aliphatic_ILV",
                     "cat_biophysics.dubreuil2019.IUP20_f",
                     "cat_biophysics.d2p2.Lsegmax",
                     "cat_interactions.string.cent_closeness",
                     "cat_interactions.string.cent_pagerank",
                     "cat_functions.go.MF_nucleotide_binding",
                     "cat_biophysics.pfam.HMM_none",
                     "cat_biophysics.superfam.supfam_none",
                     "cat_biophysics.uniprot.f_S",
                     "cat_biophysics.uniprot.f_N",
                     "cat_biophysics.uniprot.f_G",
                     "cat_biophysics.uniprot.f_V",
                     "cat_transcriptomics.sgd.GGT",
                     "cat_biophysics.leuenberger2017.Tm_medium",
                     "cat_transcriptomics.brar2012.TE_metaphase_I",
                     "cat_transcriptomics.barton2010.cost_yeast_sulphur_relative",
                     "cat_functions.costanzo2010.Metabolism_mitochondria",
                     "cat_functions.go.BP_ribosome_biogenesis",
                     "cat_functions.go.BP_biological_process",
                     "cat_phenotypes.vanleeuwen2020.essential_core",
                     "cat_biophysics.dubreuil2019.uniprot_long",
                     "cat_biophysics.superfamilies.supfam_beta_beta_alpha_zinc_fingers")
                     #"cat_genomics.peter2018.snp_full")

library(see)
col_means <- lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
col_zeros = lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), function(x){ return(0) })
YEASTOMICS.nona <- replace_na(YEASTOMICS, col_means)
#YEASTOMICS.nona <- replace_na(YEASTOMICS, col_zeros)


formula.bestprop = str_c("cat_",BEST_PROP$col_prop, collapse = ' + ')
formula.bestfeat = str_c("cat_",BEST_FEAT$col_feat, collapse = ' + ')
formula.best = paste0(formula.bestprop,"+",formula.bestfeat,collapse = ' + ')
#save.image(here("output","checkpoint-fitting-evorate.rdata"))

options(dplyr.width=Inf)

#### FUNGI DATA: ERfull ~ MPC  ####
m0 = lm(data=YEASTOMICS.nona, EVO.FULL ~ MPC )
decompose_variance(m0)
formula.m0 = as.formula(paste0("EVO.FULL ~ offset(",coef(m0)[2],"*MPC) "))
m00 = lm(data=YEASTOMICS.nona, formula.m0 )
decompose_variance(m00)

m1 = step(m00,scope = as.formula(paste0(". ~ . +",formula.bestprop)), direction = 'forward')
m2 = step(m00,scope = as.formula(paste0('EVO.FULL ~ MPC +',formula.bestfeat)), direction = 'forward',na.action = "na.exclude")
m3=step(m00,scope = as.formula(paste0('EVO.FULL ~ MPC +',formula.best)), direction = 'forward')

decompose_variance(m1)
decompose_variance(m2)
decompose_variance(m3)

Mnull = lm(data=YEASTOMICS.nona,EVO.FULL~1)
m4=step(Mnull,scope = as.formula(paste0('EVO.FULL ~ ',formula.best)), direction = 'forward')
decompose_variance(m4)


theme_update(axis.text.x = element_text(size=8,angle = 90),aspect.ratio=NULL)
Psel = gsub("TRUE$","",names(coef(m1))) %>% grep(pattern='cat_',v=T) %>% gsub(patt='cat_',repl="")
F3 = ggplot(BEST_PROP %>% dplyr::mutate( kept = col_prop %in% Psel) %>% dplyr::filter(kept),
            aes(x = reorder(str_trunc(property,w=50), abs(dRSS.tss*EXPECTED)),
                y = dRSS.tss*EXPECTED, fill=factor(kept))) +
  geom_bar(stat='identity',position='dodge', width = 0.6)  +
  #geom_text(aes(label=N),hjust='inward',size=5,col='black') +
  geom_hline(yintercept = c(-5:5),col='gray50', size=0.25) +
  geom_hline(yintercept = 0,col='black', size=0.25) +
  facet_wrap(~-EXPECTED,scales = 'free_y') + coord_flip() +
  ylab('% of TSS\n(variance in evolutionary rate)') + xlab('Properties') + theme(axis.text.x = element_text(size=10,hjust = 1))
plot(F3)

Fsel = gsub("TRUE$","",names(coef(m2))) %>% grep(pattern='cat_',v=T) %>% gsub(patt='cat_',repl="")
theme_update(axis.text.x = element_text(size=8,angle = 90),aspect.ratio=NULL)
F4 = ggplot(BEST_FEAT %>% dplyr::mutate( kept = col_feat %in% Fsel) %>% dplyr::filter(kept),
            aes(x = reorder(str_trunc(feature,w=50), abs(r*SIGN)),
                y = r, fill=factor(kept))) +
  geom_bar(stat='identity',position='dodge', width = 0.6)  +
  geom_hline(yintercept = 0,col='black', size=0.25) +
  #geom_text(aes(label=N),hjust='inward',size=5,col='black') +
  facet_wrap(~-SIGN,scales = 'free_y') + coord_flip() +
  ylab('Pearson correlation coefficient') + xlab('Features') + theme(axis.text.x = element_text(size=10,hjust = 1))
plot(F4)

BARS = plot_grid(F3,F4,nrow = 2) + theme(aspect.ratio = 1/2)
save_plot(BARS,filename = here('output','f3.6-bars-model.pdf'), base_height = 12, base_width = 20)

#### SNP DATA: ERfull ~ MPC + SNP ####

YEASTOMICS.snp = YEASTOMICS %>% dplyr::filter( !is.na(cat_genomics.peter2018.snp_full) )
col_means <- lapply(YEASTOMICS.snp %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)

col_zeros <- lapply(YEASTOMICS.snp %>% dplyr::select(where(is.numeric)), function(x){ return(0) })
YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_means) %>%
            dplyr::rename(evofull.resid = .resid, evofull.fitted=.fitted, evofull.se.fit=.se.fit, evofull.hat=.hat, evofull.sigma=.sigma, evofull.cooksd=.cooksd,evofull.std.resid=.std.resid) %>%
            broom::augment_columns(x=lm(cat_genomics.peter2018.snp_full~MPC,data=.)) %>%
            dplyr::rename(snpfull.resid = .resid, snpfull.fitted=.fitted, snpfull.se=.se.fit, snpfull.hat=.hat, snpfull.sigma=.sigma, snpfull.cooksd=.cooksd,snpfull.std.resid=.std.resid)

#YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_zeros)

m0.1 = lm(data=YEASTOMICS.snp_nona, paste0("EVO.FULL ~ MPC") )
m0.2 = lm(data=YEASTOMICS.snp_nona, paste0("EVO.FULL ~ cat_genomics.peter2018.snp_full") )
decompose_variance(m0.1)
decompose_variance(m0.2)

M0 = lm(data=YEASTOMICS.snp_nona, "EVO.FULL ~ MPC " )
decompose_variance(M0)
formula.M0 = as.formula(paste0("EVO.FULL ~ offset(",coef(M0)[2],"*MPC) "))
M00 = lm(data=YEASTOMICS.snp_nona, formula.M0)
decompose_variance(M00)

M1 = step(M00,scope = as.formula(paste0(formula.M0," + ",formula.bestprop)), direction = 'forward')
M2 = step(M00,scope = as.formula(paste0(formula.M0,' + ',formula.bestfeat)), direction = 'forward',na.action = "na.exclude")
M3=step(M00,scope = as.formula(paste0(formula.M0,' + ',formula.best)), direction = 'forward')

decompose_variance(M1)
decompose_variance(M2)
decompose_variance(M3)

#### SNP DATA: ERfull ~ MPC (SNP DATA) ####

YEASTOMICS.snp = YEASTOMICS %>% dplyr::filter( !is.na(cat_genomics.peter2018.snp_full) )
col_means <- lapply(YEASTOMICS.snp %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
col_zeros <- lapply(YEASTOMICS.snp %>% dplyr::select(where(is.numeric)), function(x){ return(0) })
YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_means)
#YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_zeros)

M0.0 = lm(data=YEASTOMICS.snp_nona, EVO.FULL ~ MPC)
decompose_variance(M0.0)
formula.M0.0 = as.formula(paste0("EVO.FULL ~ offset(",coef(M0.0)[2],"*MPC)"))
M00.0 = lm(data=YEASTOMICS.snp_nona, formula.M0.0 )
decompose_variance(M00.0)

M1.0 = step(M00.0,scope = as.formula(paste0(formula.M0.0 ," + ",formula.bestprop)), direction = 'forward')
M2.0 = step(M00.0,scope = as.formula(paste0(formula.M0.0 ,' + ',formula.bestfeat)), direction = 'forward')
M3.0=step(lm(data=YEASTOMICS.snp_nona,formula = paste0(formula.M0.0,' + ',formula.best)), direction = 'backward')
#M3.1=step(M00.0,scope = as.formula(paste0(formula.M0.0 ,' + ',formula.best)), direction = 'forward')

decompose_variance(M1.0)
decompose_variance(M2.0)
decompose_variance(M3.0)
decompose_variance(M3.1)


Mnull = lm(data=YEASTOMICS.snp_nona,EVO.FULL~1)
M3.2=step(M00.0,scope = as.formula(paste0('EVO.FULL ~ ',formula.best)), direction = 'forward')

M3.3=step(lm(data=YEASTOMICS.snp_nona,paste0('EVO.FULL ~ ',formula.best)), direction = 'forward')
M3.4=step(Mnull,scope = as.formula(paste0('EVO.FULL ~ MPC + ',formula.best)), direction = 'forward')


M3.5=step(lm(data=YEASTOMICS.snp_nona,paste0('EVO.FULL ~ cat_genomics.peter2018.snp_full +',formula.best)), direction = 'backward')

decompose_variance(M3.2)
decompose_variance(M3.3)
decompose_variance(M3.4)
decompose_variance(M3.5)
summary(aov(M3.5))[[1]]


M3.1=step(M00.0,scope = as.formula(paste0(formula.M0.0 ,' + ',formula.best)), direction = 'forward',)
decompose_variance(M3.1)

DF = summary(aov(M3.1))[[1]]

sigvarsel = head(rownames(DF[DF$`Pr(>F)`<0.05,]),-1) %>% str_trim
sum( sigvarsel %in% paste0('cat_',ALL_PROP$col_prop) )
sum( sigvarsel %in% paste0('cat_',ALL_FEAT$col_feat) )
intersect(sigvarsel,paste0('cat_',ALL_PROP$col_prop))
intersect(sigvarsel,paste0('cat_',ALL_FEAT$col_feat))

varsel = head(rownames(DF),-1) %>% str_trim

sum( varsel %in% paste0('cat_',ALL_PROP$col_prop) )
sum( varsel %in% paste0('cat_',ALL_FEAT$col_feat) )
intersect(varsel,paste0('cat_',ALL_PROP$col_prop))
intersect(varsel,paste0('cat_',ALL_FEAT$col_feat))

rownames(DF)[DF$`Sum Sq` > 1]
cor( predict(M3.0), M3.0$model$EVO.FULL)

df = summary(aov(m3))[[1]]
rownames(df)[df$`Sum Sq` > 1]
cor( predict(m3), m3$model$EVO.FULL)

plot(y=M3.0$model$EVO.FULL,predict(M3.0)) # M3.0$model$EVO.FULL
plot(predict(m3), m3$model$EVO.FULL)

plot(M3.0)
msnp =  lm(data=YEASTOMICS.snp_nona, EVO.FULL ~ cat_genomics.peter2018.snp_full)
decompose_variance(msnp)



fit = step(M00.0,scope = as.formula(paste0(formula.M0.0 ," + ",'cat_functions.go.MF_nucleotide_binding')), direction = 'forward')
fit = step(M00.0,scope = as.formula(paste0(formula.M0.0 ," + ",'cat_phenotypes.vanleeuwen2020.essential_core')), direction = 'forward')


f0=lm(data=YEASTOMICS.nona, formula = EVO.FULL ~MPC)
decompose_variance(f0)

fbest=lm(data=YEASTOMICS.nona, formula = EVO.FULL ~MPC)


f1=lm(data=YEASTOMICS.nona, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) + cat_functions.go.MF_nucleotide_binding + cat_functions.go.MF_molecular_function)
f2=lm(data=YEASTOMICS.nona, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) + cat_functions.go.MF_molecular_function + cat_functions.go.MF_nucleotide_binding)
fa=lm(data=YEASTOMICS.nona, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) + cat_functions.go.MF_molecular_function)
fb=lm(data=YEASTOMICS.nona, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) + cat_functions.go.MF_nucleotide_binding)

rss0=deviance(f0)
rss0-deviance(f1)
rss0-deviance(f2)
rss0-deviance(fa)
rss0-deviance(fb)
summary(aov(fa))[[1]]
summary(aov(fb))[[1]]
summary(aov(f1))[[1]]
summary(aov(f2))[[1]]

ntbind = YEASTOMICS.nona$cat_functions.go.MF_nucleotide_binding
mfunk = YEASTOMICS.nona$cat_functions.go.MF_molecular_function
sum(ntbind)
sum(mfunk)
sum(ntbind | mfunk)
sum(ntbind & mfunk)
table(ntbind,mfunk)
mu = mean_(YEASTOMICS.nona$EVO.FULL)
tss = sum( (YEASTOMICS.nona$EVO.FULL-mu)^2 )

tssa=sum( (YEASTOMICS.nona$EVO.FULL[mfunk]-mu)^2 )
tssb=sum( (YEASTOMICS.nona$EVO.FULL[ntbind]-mu)^2 )
tssc = sum( (YEASTOMICS.nona$EVO.FULL[ntbind|mfunk]-mu)^2 )

ess0 =  sum( (fitted(f0)-mu)^2 )

essa = sum( (fitted(fa)-mu)^2 ) - ess0
essb = sum( (fitted(fb)-mu)^2 )- ess0
essc = sum( (fitted(f1)-mu)^2 )- ess0
essd = sum( (fitted(f2)-mu)^2 )- ess0


essa/tss
essb/tss

(essa+essb)/tss
essc/tss
rss0 =  deviance(f0)
rssa = rss0- deviance(fa)
rssb = rss0 - deviance(fb)
rssc = rss0 - deviance(f1)
rssd = rss0 - deviance(f2)
fit = step(M00.0,scope = as.formula(paste0(formula.M0.0 ," + ",'cat_phenotypes.vanleeuwen2020.essential_core')), direction = 'forward')
var(M00.0$model$EVO.FULL) * nrow(M00.0$model)

sigvar = M3[-42,] %>%
         arrange(desc(`Sum Sq`)) %>%
          filter(`Pr(>F)`<0.05) %>%
          rename(SS=`Sum Sq`) %>%
          separate(variable, sep='\\.', into=c('categories','source','variable')  ) %>%
          mutate(variable=str_trim(variable), type = ifelse(variable %in% BEST_PROP$property, 'property','feature'))

F3.7=ggplot(sigvar, aes(y=reorder(variable,SS), x=SS,fill=type)) +
  geom_bar(stat='identity',orientation = 'y', width=0.5) +
  geom_text(aes(label=reorder(paste0(round(SS,1)," ",variable),SS), x=SS+1),col='black',hjust=0) +
  ylab('Variable') + scale_fill_metro() + scale_x_discrete(expand=c(0,0)) +
  theme_cowplot() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
plot(F3.7)
save_plot(F3.7, filename = here('output','FIG3.7-significant-variables-stepwise-reg.pdf'), base_height = 6, base_width = 8)

library(treemap)
treemap(sigvar, index = 'variable', vSize = 'SS',inflate.labels = T, aspRatio = 1, )



cor(x=YEASTOMICS.snp_nona$EVO.FULL, y=YEASTOMICS.snp_nona$cat_genomics.peter2018.snp_full)
cor(x=YEASTOMICS.snp_nona$evofull.resid, y=YEASTOMICS.snp_nona$snpfull.resid)
cor(x=YEASTOMICS.snp_nona$evofull.resid, y=YEASTOMICS.snp_nona$cat_genomics.peter2018.snp_full)

msnp = lm(data=YEASTOMICS.snp_nona, EVO.FULL~ cat_genomics.peter2018.snp_full)
aov(msnp)
RSS.snp = sum(residuals(msnp)^2)
TSS.snp=var(YEASTOMICS.snp_nona$EVO.FULL) * nrow(YEASTOMICS.snp_nona)
RSS.snp/TSS.snp

decompose_variance(msnp)
sqrt(223.7/671.1275)
0.577^2

snp0 = lm(data=YEASTOMICS.snp_nona, EVO.FULL ~ MPC)
m0.snp = lm(data=YEASTOMICS.snp_nona, paste0("EVO.FULL ~ offset(",coef(snp0)[2],"*MPC) + cat_genomics.peter2018.snp_full"))
mm.snp = lm(data=YEASTOMICS.snp_nona, "EVO.FULL ~ MPC + cat_genomics.peter2018.snp_full")
decompose_variance(m0.snp)
decompose_variance(mm.snp)
cor(y=YEASTOMICS.snp_nona$EVO.FULL, x=coef(snp0)[2] * YEASTOMICS.snp_nona$MPC + coef(mm.snp)[3]*YEASTOMICS.snp_nona$cat_genomics.peter2018.snp_full)



####
MANUAL = YEASTOMICS.nona %>% dplyr::select(EVO.FULL,manual.selection)

mm = lm(MANUAL, formula=EVO.FULL ~ .)
decompose_variance(mm)
SS = var(MANUAL$EVO.FULL) * nrow(MANUAL)
deviance(mm)/SS
(SS-deviance(mm))/SS
df = summary(aov(mm))[[1]]


summary(aov(mm))[[1]]

