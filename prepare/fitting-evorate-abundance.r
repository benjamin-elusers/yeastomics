#### SETUP ####
source("https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/__setup_yeastomics__.r")
source("https://raw.githubusercontent.com/clauswilke/dviz.supp/master/R/dviz.supp.R")
source("https://raw.githubusercontent.com/clauswilke/dviz.supp/master/R/themes.R")
empty_theme <- theme_dviz_open(12, rel_small = 1, rel_large = 1) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = grid::unit(0, "pt")
  )
# theme_set( theme_modern(base_size = 16, axis.text.angle = 45,legend.position = 'none') + theme(aspect.ratio=1) )

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

ABUNDANCE = load.ho2018.data() %>%
  dplyr::select(orf, MPC = mean.mpc, gfp=GFP.avg,ms=MS.avg) %>% # select the average expression
  group_by(orf) %>% dplyr::summarise(across(where(is.numeric),log10)) %>% # apply log10 to expression
  rowwise() %>% mutate( noval = mean(is.na(c_across(where(is.numeric)))) ) %>% # check how many expression values
  dplyr::filter(noval<1 || !is.na(MPC) ) %>% dplyr::select(-noval) # remove no values

CLADE = get_clade_data(g1='schizo',g2='sacch.wgd',rate = 'ratio')
#FUNGI = readRDS(here("data","fungi-evodata.rds")) %>% filter( !is.na(PPM) & !is.na(MPC) )
# EVOLUTIONARY RATE ON FUNGI
FUNGI = load.dubreuil2021.data(1) %>%
        dplyr::select(c(starts_with('EVO.'),'PPM','ORF','UNIPROT')) %>%
        ungroup() %>%
        left_join(ABUNDANCE,by=c('ORF'='orf')) %>%
        filter( !is.na(PPM) & !is.na(MPC) )

FUNGI.NORM = FUNGI %>% mutate(rel_EVO.FULL=log2(EVO.FULL),
                              rel_EVO.DISORDER=log2(EVO.DISORDER),
                              rel_EVO.NOT_DISORDER=log2(EVO.NOT_DISORDER),
                              rel_EVO.DOMAINS=log2(EVO.DOMAINS),
                              rel_EVO.ANTI_DOMAIN=log2(EVO.ANTI_DOMAIN),
                              rel_EVO.PDB=log2(EVO.PDB),
                              rel_EVO.SURFACE=log2(EVO.SURFACE),
                              rel_EVO.BURIED=log2(EVO.BURIED))

# STRAINS = readRDS(here("data","yeast_strains-evodata.rds")) %>%
#           filter( !is.na(PPM) & !is.na(MPC) ) %>% distinct %>%
#           filter(!(is.dup(UNIPROT) & is.na(EVO.FULL)))


STRAINS = readRDS(here("data","PROTEIN-EVO-FUNGI-SNP.rds")) %>%
  dplyr::select(c(starts_with(c('EVO.','SNP.')),'PPM','ORF','UNIPROT','IS_FUNGI','IS_STRAINS')) %>%
  group_by(ORF) %>% dplyr::mutate(across(starts_with("SNP."),log)) %>% # apply log10 to SNP rate4site
  left_join(ABUNDANCE,by=c('ORF'='orf')) %>%
  filter( !is.na(PPM) & !is.na(MPC) ) %>%
  filter(!(is.dup(UNIPROT) & is.na(EVO.FULL)))


### STRAINS CONTAIN ALL DATA ABOUT THE FUNGI LINEAGE EVOLUTIONARY RATE

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

INPUT = FUNGI.NORM
evocols = INPUT %>% dplyr::select(starts_with('rel_EVO.')) %>% colnames
#EVO = left_join(INPUT,PROP) %>% ungroup()
INPUT$EVO.FULL = log10(INPUT$EVO.FULL)

#### EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ####
Y = "EVO.FULL" # mean Evolutionary rate (full sequence)
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
#M2=make_logistic_fit(input,only.params = F)
M3=make_poly_fit(input,only.params = F,deg=3)
M4=make_expo_fit(input,only.params = F)
M = bind_rows(M1,M2,M3,M4)

M = M1

selected.par =c('model','RSS','ESS','TSS')
fit.params = bind_rows(
    make_linear_fit(input,only.params = T)[selected.par],
#    make_logistic_fit(input,only.params =T)[selected.par],
#    make_poly_fit(input,only.params = T,deg=3)[selected.par],
#    make_expo_fit(input,only.params = T)[selected.par]
  )
fit.params


#### _ EVO.vs.MPC SCATTERPLOT WITH FIT   ####
x1=nearest(0.5,XX,n=1,v=T)
M.x1 = M[ M[X] == x1 ,]
pXY=ggplot(input,aes_string(y=Y,x=X)) +
  geom_point(size=2.5,shape=19,alpha=0.5,color='gray70',stroke=0) +
  stat_density2d(size=0.5,color='gray20') +
  #  geom_text(data=CC.ER, aes(label=sprintf(" r %.2f \n p %.1e \n N %s ",estimate,p.value,n)),
  #            x=Inf,y=Inf,hjust='inward',vjust='inward',size=3) +
  ylab('mean Evolutionary Rate (log2)') + xlab('Median Protein Abundance (log10 mpc)') +
  geom_hline(yintercept = mu.y, col='black',linetype=1,size=1) + # mean
  ylim(rg.y) + theme(legend.position = 'none',aspect.ratio = 1) +
  ggpubr::grids()


pXY.fit = pXY+
  geom_line(M,mapping=aes(y=.fitted,col=model),size=1,show.legend = F) +
  geom_text(M.x1 ,mapping=aes(label=model,col=model,y=.fitted-0.2),x=0.8,size=5,check_overlap = T,show.legend = F) +
  stat_smooth(method = 'loess', fullrange=T, span=1.2, se = F, col='gray40',size=1) + annotate('text',y=1.9,x=0.8,label='Loess',col='gray40',size=5)
  #coord_cartesian(xlim = c(-2.5, 4.5), ylim = c(0, 3.5),expand = F) +
pXY.fit
ggsave(plot=pXY.fit, filename="~/Desktop/logER-vs-abundance.png",scale=1.2)

MF_nuc_bind = PROP$ORF[PROP$cat_functions.go.MF_nucleotide_binding]
MF_unknown  = PROP$ORF[PROP$cat_functions.go.MF_molecular_function]
linfit = geom_line(M1,mapping=aes(y=.fitted,col=model),size=1,show.legend = F)
PXY=plot_grid(pXY.fit,
          pXY + geom_point(data=input %>% filter(ORF %in% MF_unknown), color='green',size=1)+linfit,
          pXY + geom_point(data=input %>% filter(ORF %in% MF_nuc_bind) , color='blue',size=1) +linfit,
          nrow=3) + ggiraphExtra::theme_clean2()

#### _ RESIDUAL EVO.vs.MPC SCATTERPLOT (EVO.FULL) ####
col.samples = PROP %>% dplyr::select(contains(samples)) %>% colnames
library(ggiraph)
library(ggiraphExtra)
library(gghighlight)
fitdata = left_join(M1,PROP, by='ORF')# %>%
          #bind_rows(M1 %>% mutate(property='all_orthologs', has_prop=T) ) %>%
          #dplyr::filter( property %in% c(samples,'all_orthologs')) %>% mutate(property = as.factor(property))
#res.ortho = res.data %>% filter(property=='all_orthologs')
f0=lm(data=fitdata, formula = rel_EVO.FULL ~MPC)
f1=lm(data=fitdata, formula = rel_EVO.FULL ~ offset(coef(f0)[2]*MPC) + cat_functions.go.MF_nucleotide_binding + cat_functions.go.MF_molecular_function)

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
res.evo = STRAINS %>%
    broom::augment_columns(x=lm(SNP.FULL~MPC,data=.)) %>%
    dplyr::rename(snpfull.resid = .resid, snpfull.fitted=.fitted, snpfull.se=.se.fit, snpfull.hat=.hat, snpfull.sigma=.sigma, snpfull.cooksd=.cooksd,snpfull.std.resid=.std.resid)  %>%
    broom::augment_columns(x=lm(EVO.FULL~MPC,data=.)) %>%
    dplyr::rename(evofull.resid = .resid, evofull.fitted=.fitted, evofull.se.fit=.se.fit, evofull.hat=.hat, evofull.sigma=.sigma, evofull.cooksd=.cooksd,evofull.std.resid=.std.resid)


FXA = make_scatterplot(res.evo,xvar='SNP.FULL',yvar='EVO.FULL',
                       labx='SNP Evolutionary Rate',laby='mean Evolutionary Rate',
                       theme2use = theme(),txtcorner = 'topleft') +
      theme(legend.position='none',aspect.ratio = 1)


FXB = make_scatterplot(res.evo,xvar='snpfull.resid',yvar='evofull.resid',
                       labx='residual SNP Evolutionary Rate',laby='residual Evolutionary Rate',
                       theme2use = theme(),txtcorner = 'bottomright') +  theme(legend.position='none',aspect.ratio = 1)

res.evo.clade = res.evo %>% filter(ORF %in% CLADE$orf)


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

prop.lin = left_join(M1,prop)
colab = c('PPM','MPC','gfp','ms')

FIT =  prop.lin %>%
  group_by(col_prop) %>%
  mutate(
           .fitted.prop = .fitted + mean(.resid),
           .resid.prop = .resid + mean(.resid),
           N=n(),
           N_RS_pos = sum(.resid>0),
           N_RS_neg = sum(.resid<0),
           RS     = sum(.resid),
           RS.log = sign(RS) * log(abs(RS)),
           RS_pos = sum(.resid[.resid>0]),
           RS_neg = sum(.resid[.resid<0]),
           RS_avg = mean_(.resid),
           RSS.pc  = RSS / TSS,
           ESS.pc  = ESS / TSS,
           TSS.p  = sum( (EVO.FULL - mu.y)^2 ),
           ESS.p  = sum( (.fitted.prop-mu.y)^2),
           RSS.p = TSS.p - ESS.p,
           tss.pc = 100*TSS.p/TSS,
           ess = sum( (.fitted-mu.y)^2),
           rss = TSS.p - ess,
           dRSS  = RSS.p-rss,
           dRSS.pc  = 100*dRSS / TSS,
           rss.pc  = (100*RSS.p / TSS.p) * sign(RS_avg),
           dESS  = ESS.p-ess,
           dESS.pc = 100*dESS / TSS,
           ess.pc  = (100*ESS.p / TSS.p) * sign(RS_avg),
           EXPECTED = sign(RS_avg)
  ) %>%
  arrange(desc(RSS.pc)) %>% #,desc(.res_snp.RS.pc)) %>%
  dplyr::select(categories,source,property, col_prop, EXPECTED,#colab, evocols,
                N,N_RS_pos,N_RS_neg,
                RSE,LL,TSS,ESS,ESS.pc,RSS,RSS.pc,
                TSS.p, ESS.p,ess , RSS.p,rss, tss.pc, ess.pc, rss.pc, dRSS, dESS, dESS.pc, dRSS.pc,
                starts_with('RS')) %>%
  distinct() %>%
  group_by(categories) %>% mutate(rRS = dense_rank(rss.pc) , ncat=n()) %>%
  arrange(desc(abs(dRSS.pc)),desc(TSS.p))


#FIT %>% dplyr::filter(col_prop %in% selected.best)
# rank.prop = group_by(FIT, property,col_prop) %>%
#             summarise(n=N,
#                       tss=mean_(TSS),
#                       ESS = mean(ESS),
#                       ess = mean(ESS.pc),
#                       RSS = mean(RSS),
#                       rss = mean(RSS.pc),
#                       ess.p=mean_(ESS.p),
#                       rss.p = mean_(RSS.p),
#                       res.avg=mean_(RS_avg),
#                       dRSS = rss.p-RSS,
#                       drss = 100*dRSS/tss,
#                       dESS = ess.p-ESS,
#                       dess = 100*dESS/tss) %>%
#             dplyr::filter(n>2) %>%
#             arrange(desc(dess),desc(rss))
# best.prop = rank.prop %>% dplyr::filter(abs(drss)>1 & dess>1) %>% pull(col_prop) %>% paste0('cat_',.)
df.best.prop = FIT %>% dplyr::filter(tss.pc > 1 & abs(dRSS.pc)>1 & abs(dESS.pc)>1)
# best.prop = df.best.prop %>% pull(col_prop) %>% paste0('cat_',.)
#
# mall.prop=lm(data=fitdata[,c('EVO.FULL','MPC',best.prop)], formula = EVO.FULL ~MPC)
#
#
# m0.prop=lm(data=fitdata[,c('EVO.FULL','MPC',best.prop)], formula = EVO.FULL ~MPC)
# m1.prop <- lm(data=fitdata,  I(EVO.FULL-coef(m0.prop)[1]) ~  offset(coef(m0.prop)[2]*MPC) + cat_functions.go.MF_nucleotide_binding + cat_functions.go.MF_molecular_function)
# m.bestprop=lm(".resid ~ offset(coef(m0.prop)[2] * MPC) + .", data =fitdata[,c('.resid','MPC',best.prop)])
# add.best = str_c(best.prop, collapse = ' + ')
# fitprop.final = step(m0.prop,scope = as.formula(paste0('. ~ . +',add.best)), direction = 'forward')
# selected.best = gsub("TRUE$","",names(coef(fitprop.final))) %>% grep(pattern='cat_',v=T) %>% gsub(patt='cat_',repl="")
#
# theme_update(axis.text.x = element_text(size=8,angle = 90),aspect.ratio=NULL)
# F3 = ggplot(FIT %>% dplyr::filter(col_prop %in% selected.best),
#             aes(x = reorder(str_trunc(property,w=50), -dRSS.pc*EXPECTED),
#                 y = dRSS.pc*EXPECTED)) +
#   geom_bar(stat='identity',position='dodge')  +
#   #geom_text(aes(label=N),hjust='inward',size=5,col='black') +
#   geom_hline(yintercept = c(-5:5),col='gray50', size=0.25) +
#   facet_wrap(~EXPECTED,scales = 'free_x') + coord_flip() +
#   ylab('% of TSS\n(variance in evolutionary rate)') + xlab('Properties') + theme(axis.text.x = element_text(size=10,hjust = 1))
# plot(F3)
# save_plot(plot = F3, here("output",'f3.5E-barplot-best-prop-selected.pdf'), base_height=12, base_width=20)
#





# prop.res = res.data %>% dplyr::select(EVO.FULL, MPC, starts_with('cat_'))
# f0=lm(data=prop.res, formula = EVO.FULL ~ MPC)
# # PRETTY LONG TO FIT ALL PROPERTIES
# tic('Fitting all proteome properties...')
# fprop=lm(data=prop.res, formula = EVO.FULL ~ offset(coef(f0)[2]*MPC) + . )
# deviance(fprop)
# toc()
# final = step(fprop, direction = 'backward')

#
# # BEST PROPERTIES - FIT WITH FIXED COEFFICIENT
# ggplot(prop.res %>% dplyr::filter(property %in% c('all_orhtologs',best.prop)), aes(y=EVO.FULL,x=MPC,col=property)) +
#   geom_point(show.legend=F) +
#   geom_line(aes(y = fitted.f0), size = 1,show.legend=F)
# summary(fprop)

#
# BEST = EVO %>% dplyr::select(EVO.FULL,PPM,contains(best.prop,ignore.case = F))
# MBEST = lm(EVO.FULL~PPM+.,data=BEST)
# COMPACT = step(MBEST)

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
PROP_FEAT = left_join(PROP,features.norm)
YEASTOMICS = left_join(M1,PROP_FEAT)
dim(YEASTOMICS)

BEST_PROP =  YEASTOMICS %>%
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
              TSS.p  = sum( (EVO.FULL - mu.y)^2 ),
              ESS.p  = sum( (.fitted.prop-mu.y)^2),
              RSS.p = TSS.p - ESS.p,
              tss.pc = 100*TSS.p/TSS,
              ess = sum( (.fitted-mu.y)^2),
              rss = TSS.p - ess,
              dRSS  = RSS.p-rss,
              dRSS.pc  = 100*dRSS / TSS,
              rss.pc  = (100*RSS.p / TSS.p) * sign(RS_avg),
              dESS  = ESS.p-ess,
              dESS.pc = 100*dESS / TSS,
              ess.pc  = (100*ESS.p / TSS.p) * sign(RS_avg)
  ) %>%
  dplyr::select(categories,source,property, col_prop,
                N, RS_avg, EXPECTED, TSS,ESS,RSS,
                TSS.p, ESS.p,ess , RSS.p,rss,
                tss.pc, ess.pc, rss.pc,
                dRSS, dESS, dESS.pc, dRSS.pc) %>%
  distinct() %>%
  dplyr::filter(tss.pc > 1 & abs(dRSS.pc)>1 & abs(dESS.pc)>1)

BEST_FEAT = YEASTOMICS %>%
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
  arrange(rk) %>% relocate(rk) %>%
  filter(abs(r) > 0.2  & !grepl("byrne2005",col_feat) & !grepl("snp_",col_feat) ) %>%
  mutate(SIGN=sign(r))
library(see)

col_means <- lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
col_zeros = lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), function(x){ return(0) })
YEASTOMICS.nona <- replace_na(YEASTOMICS, col_means)
#YEASTOMICS.nona <- replace_na(YEASTOMICS, col_zeros)

formula.bestprop = str_c("cat_",BEST_PROP$col_prop, collapse = ' + ')
formula.bestfeat = str_c("cat_",BEST_FEAT$col_feat, collapse = ' + ')
formula.best = paste0(formula.bestprop,"+",formula.bestfeat,collapse = ' + ')

best_var = c(BEST_PROP$col_prop,BEST_FEAT$col_feat)

cor(YEASTOMICS$MPC,YEASTOMICS[,paste0("cat_",best_var)],use='complete','spearman')
tmp=abs(cor(YEASTOMICS$MPC,YEASTOMICS[,paste0("cat_",best_var)],use='complete','spearman'))<0.4
no_mpc_var = colnames(tmp)[!is.na(tmp) & unlist(tmp)==T]

formula.nompc = str_c(no_mpc_var, collapse = ' + ')

#### FUNGI DATA: ERfull ~ MPC  ####

mnull = lm(data=YEASTOMICS.nona, rel_EVO.FULL ~ 1 )
decompose_variance(mnull)

m0 = lm(data=YEASTOMICS.nona, rel_EVO.FULL ~ MPC )
decompose_variance(m0)
formula.m0 = as.formula(paste0("rel_EVO.FULL ~ offset(",coef(m0)[2],"*MPC) "))
m00 = lm(data=YEASTOMICS.nona, formula.m0 )
decompose_variance(m00)

m1 = step(m00,scope = as.formula(paste0(". ~ . +",formula.bestprop)), direction = 'forward')
decompose_variance(m1)

m2 = step(m00,scope = as.formula(paste0('rel_EVO.FULL ~ MPC +',formula.bestfeat)), direction = 'forward',na.action = "na.exclude")
decompose_variance(m2)

m3=step(m00,scope = as.formula(paste0('rel_EVO.FULL ~ MPC +',formula.best)), direction = 'forward')
decompose_variance(m3)


m4=step(mnull,scope = as.formula(paste0('rel_EVO.FULL ~ ',formula.nompc)), direction = 'forward')
decompose_variance(m4)
anova(m4)



theme_update(axis.text.x = element_text(size=8,angle = 90),aspect.ratio=NULL)
Psel = gsub("TRUE$","",names(coef(m1))) %>% grep(pattern='cat_',v=T) %>% gsub(patt='cat_',repl="")
F3 = ggplot(BEST_PROP %>% dplyr::mutate( kept = col_prop %in% Psel) %>% dplyr::filter(kept),
            aes(x = reorder(str_trunc(property,w=50), abs(dRSS.pc*EXPECTED)),
                y = dRSS.pc*EXPECTED, fill=factor(kept))) +
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
YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_means)
#YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_zeros)

m0 = lm(data=YEASTOMICS.snp, rel_EVO.FULL ~ MPC + cat_genomics.peter2018.snp_full )
decompose_variance(m0)
formula.m0 = as.formula(paste0("rel_EVO.FULL ~ offset(",coef(m0)[2],"*MPC) + offset(",coef(m0)[3],"*cat_genomics.peter2018.snp_full)"))
m00 = lm(data=YEASTOMICS.snp_nona, formula.m0 )
decompose_variance(m00)

m1 = step(m00,scope = as.formula(paste0(". ~ . +",formula.bestprop)), direction = 'forward')
decompose_variance(m1)

m2 = step(m00,scope = as.formula(paste0(formula.m0,' + ',formula.bestfeat)), direction = 'forward',na.action = "na.exclude")
decompose_variance(m2)

m3=step(m00,scope = as.formula(paste0(formula.m0,' + ',formula.best)), direction = 'forward')
decompose_variance(m3)


#### SNP DATA: ERfull ~ MPC (SNP DATA) ####

YEASTOMICS.snp = YEASTOMICS %>% dplyr::filter( !is.na(cat_genomics.peter2018.snp_full) )
col_means <- lapply(YEASTOMICS.snp %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
col_zeros <- lapply(YEASTOMICS.snp %>% dplyr::select(where(is.numeric)), function(x){ return(0) })
YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_means)
#YEASTOMICS.snp_nona <- replace_na(YEASTOMICS.snp, col_zeros)

m0 = lm(data=YEASTOMICS.snp, rel_EVO.FULL ~ MPC)
decompose_variance(m0)
formula.m0 = as.formula(paste0("rel_EVO.FULL ~ offset(",coef(m0)[2],"*MPC)"))
m00 = lm(data=YEASTOMICS.snp_nona, formula.m0 )
decompose_variance(m00)

m1 = step(m00,scope = as.formula(paste0(formula.m0," + ",formula.bestprop,'+cat_genomics.peter2018.snp_full')), direction = 'forward')
decompose_variance(m1)

m2 = step(m00,scope = as.formula(paste0(formula.m0,' + ',formula.bestfeat,'+cat_genomics.peter2018.snp_full')), direction = 'forward',na.action = "na.exclude")
decompose_variance(m2)


m3.0=step(m00,scope = as.formula(paste0(formula.m0,' + ',formula.best)), direction = 'forward')
decompose_variance(m3.0)
m3=step(m00,scope = as.formula(paste0(formula.m0,' + ',formula.best,'+cat_genomics.peter2018.snp_full')), direction = 'forward')
decompose_variance(m3)

msnp =  lm(data=YEASTOMICS.snp_nona, rel_EVO.FULL ~ cat_genomics.peter2018.snp_full)
decompose_variance(msnp)

m0snp =  lm(data=YEASTOMICS.snp_nona, rel_EVO.FULL ~ MPC + cat_genomics.peter2018.snp_full)
decompose_variance(m0snp)

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

#### FIND RELATION BETWEEN PARAMETERS
# library(ess)
# g <- fit_graph(m1$model, trace = FALSE)
# plot(g, vertex.size = 1)

selected.best = gsub("TRUE$","",names(coef(fitprop.final))) %>% grep(pattern='cat_',v=T) %>% gsub(patt='cat_',repl="")

library(correlation)
library(see)
featmat = features.norm %>%
  dplyr::select( ORF, where(is.numeric) & !contains(c('snp','byrne2005','IUP20','IUP40'))) %>%
  bind_cols(EVO.SNP_FULL = features.norm$cat_genomics.peter2018.snp_full)

norm_features = featmat %>%
  pivot_longer(cols = starts_with('cat_'),
             names_to = c('categories','source','feature'),
             names_pattern="cat_(.+)\\.(.+)\\.(.+)",
             values_to = "value") %>%
  mutate(col_feat = paste0(categories,'.',source,'.',feature)) %>%
  dplyr::filter(!is.na(value)) %>% add_count(col_feat,name='size') %>%
  mutate( is.codon  = feature %in% get.codons4tai(), is.aa = feature %in% paste0('f_',get.AA1())) %>%
  rowwise %>% mutate( feature=ifelse(is.codon ,sprintf('%s (%s)',tolower(feature), get.AA3()[Biostrings::GENETIC_CODE[feature]]),feature)) %>%
  rowwise %>% mutate( feature=ifelse(is.aa ,sprintf('f(%s)',get.AA3()[gsub("f_","",feature)]),feature) )

#results <- summary(correlation(featmat))
#plot(results)

corfeat.res = left_join(M1,norm_features) %>%
  #dplyr::filter( !grepl('^snp',feature) & !(feature %in% c('pid','PID1','score_B100','S','G','N','RLEN')) ) %>%
  group_by(col_feat,categories,source,feature) %>%
  summarise( r=scor(.resid,value,met='pearson')$estimate,
             rho=scor(.resid,value)$estimate,
             nfeat = sum.na(value,notNA=T) ) %>%
  ungroup() %>% mutate( rk = rank(-rho,ties.method = 'first') ) %>%
  arrange(rk) %>% relocate(rk)
corfeat.res

bestfeat  = corfeat.res %>% filter(abs(rho) > 0.2 & nfeat > 100 ) %>%
            mutate(SIGN=sign(rho))

evompc.features = FUNGI %>%
  left_join(features.norm) %>%
  dplyr::select(ORF,EVO.FULL,MPC,paste0('cat_',bestfeat$col_feat)) %>%
  filter(complete.cases(EVO.FULL,MPC)) %>% distinct %>%
  ungroup() %>% mutate_if(is.numeric, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))

colnames(evompc.features) = gsub("cat_","",colnames(evompc.features))
#evompc.features[is.na(evompc.features)] = 0
summary(evompc.features)

evofeat.best = evompc.features[,c('ORF','EVO.FULL','MPC',bestfeat[,c('feature','col_feat')])]
m0.feat=lm(data=evofeat.best, formula = EVO.FULL ~MPC,na.action=na.omit)
evofeat.best$EVO.FULL = residuals(m0.feat)
m1.feat=lm(formula = as.formula("EVO.FULL ~ offset(coef(m0.feat)[2] * MPC) + ."), data =evofeat.best[,-1],na.action=na.exclude)
add.bestfeat = str_c(bestfeat$col_feat, collapse = ' + ')
fitfeat.final = step(m0.feat,scope = as.formula(paste0('. ~ . +',add.bestfeat)), direction = 'forward')
selected.feat = names(coef(fitfeat.final))[-c(1:2)]

coef(fitfeat.final)
deviance(m1.feat)
deviance(fitfeat.final)
SSTotal <- var( fitfeat.final$model[,1] ) * (nrow(fitfeat.final$model)-1)

theme_update(axis.text.x = element_text(size=8,angle = 90),aspect.ratio=NULL)
F4 = ggplot(bestfeat %>% dplyr::filter(col_feat %in% selected.feat) ,
            aes(x = reorder(feature, -rho2_evo),y = rho2_evo)) +
  geom_bar(stat='identity',position='dodge')  +
  facet_wrap(~SIGN,scales = 'free_x') + coord_flip() +
  ylab('Pearson correlation coefficient') + xlab('Features') + theme(axis.text.x = element_text(size=10,hjust = 1))
plot(F4 )
save_plot(plot = F4, here("output",'f4.5-barplot-feat.pdf'), base_height=12, base_width=20)


SSTotal-deviance(fit.final)

#Fbest = evompc.features[,c('ORF','EVO.FULL','MPC',bestfeat$col_feat)]
#Pbest = fitdata %>% dplyr::select(ORF,EVO.FULL,MPC, paste0('cat_',selected.best))
#FULL = left_join(Fbest,Pbest)




SSTotal <- var( M3.final$model[,1] ) * (nrow(M3.final$model)-1)
SSTotal-deviance(M3.final)


SNP = FUNGI %>%
      left_join(features.norm) %>%
      dplyr::select(ORF,EVO.FULL,MPC,snp_full=cat_genomics.peter2018.snp_full)
M4.0=lm(data=SNP, formula = EVO.FULL ~ MPC + snp_full, na.action=na.omit)
deviance(M4.0) / SSTotal

Fbest = evompc.features[,c('ORF','EVO.FULL','MPC',bestfeat$col_feat)]
Pbest = fitdata %>% dplyr::select(ORF,EVO.FULL,MPC, paste0('cat_',selected.best))
FULL = left_join(Fbest,Pbest)

all.best = str_c(colnames(FULL)[-c(1:2)], collapse = ' + ')
M4.final=step(M3.0,scope = as.formula(paste0('. ~ . +',all.best)), direction = 'forward')

SSTotal <- var( M3.final$model[,1] ) * (nrow(M3.final$model)-1)
SSTotal-deviance(M3.final)



test2 = FUNGI %>%
  left_join(features.norm) %>%
  dplyr::select(EVO.FULL,MPC,paste0('cat_',c(bestfeat$col_feat,'genomics.peter2018.snp_full'))) %>%
  filter(complete.cases(EVO.FULL,MPC)) %>% distinct %>%
  ungroup() %>% mutate_if(is.numeric, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))

M0=lm(data=test2, formula = EVO.FULL ~MPC, na.action=na.omit)
M00=lm(data=test2, formula = EVO.FULL ~MPC+cat_genomics.peter2018.snp_full, na.action=na.omit)
deviance(M0)
deviance(M00)
test2$EVO.FULL = residuals(M0)
M1=lm(formula = as.formula("EVO.FULL~ offset(coef(M0)[2] * MPC) + ."), data =test2,na.action=na.exclude)
add.bestfeat = str_c(c(bestfeat$col_feat,"genomics.peter2018.snp_full"), collapse = ' + ')
fit.final = step(M0,scope = as.formula(paste0('. ~ . +',add.bestfeat)), direction = 'forward')
#coef(fit.final)
deviance(M1)
deviance(fit.final)


MF1 = lm(data=FEATEVO[,c(colnames(features),"MPC")], formula= MPC ~ .)
residuals(MF1)

FEATEVO = left_join(M1,FEAT) %>% ungroup()
featevo = left_join(M1,feat) %>% ungroup()
correlate(FEATEVO) %>% dplyr::filter(var2 == '.resid') %>% arrange(coef_corr) %>% print(n=245)


# additive.model = function(M0,addTerm,fulldata){
#   ymean = mean(M0$model[,1])
#   yname = names(M0$model)[1]
#   nXY = sum(complete.cases(M0$model))
#
#   cat("------------\n")
#   old.formula = deparse(formula(M0))
#   cat("Previous Model: ",old.formula,"\n")
#   print(M0)
# #  cat("Estimated coefficients :",coef(M0),"\n")
#   RSS0 = deviance(M0)
#   ESS0 = sum( (fitted(M0)-ymean)^2 )
#   cat("RSS0: ",RSS0,"\n")
#   cat("ESS0: ",ESS0,"\n")
#
#   fit0 = broom::augment_columns(x=M0,data = fulldata) %>%
#     mutate(ESS=sum_( (.fitted-ymean)^2 ), TSS=sum_( (!!sym(yname)-ymean)^2 ), RSS=TSS-ESS,
#            RSE=sqrt( (1 / (nXY-2)) * RSS), AIC = AIC(M0), BIC=BIC(M0), LL = readr::parse_number(as.character(logLik(M0))))
#
#   cat("------------\n")
#
#   if(!missing(addTerm)){
#     cat("Adding variable :",addTerm,"\n")
#     target = '.resid ~ '
#     depvar = paste0( sprintf("offset(%.4f * %s)",coef(M0)[-1],variable.names(M0)[-1]), collapse=" + " )
#     intercept = sprintf("offset(%s)",coef(M0)[1])
#     new.formula = paste(target, depvar, " + ", addTerm)
#     cat("New Model: ",new.formula,"\n")
#
#     M1 = update(M0, as.formula(new.formula),data=fit0)
#     RSS1 = deviance(M1)
#     ESS1 = sum( (fitted(M1)-ymean)^2 )
#
#     cat("Estimated coefficients :",coef(M1),"\n")
#     cat("RSS1: ",RSS1,"(d=", RSS0-RSS1,")\n")
#     cat("ESS1: ",ESS1,"(d=", ESS1-ESS0,")\n")
#   }else{
#     M1=M0
#   }
#   return(M1)
# }

evompc.features = FUNGI %>%
                  dplyr::select(ORF,EVO.FULL,MPC) %>%
                  left_join(features.norm) %>%
                  filter(complete.cases(EVO.FULL,MPC))
colnames(evompc.features) = gsub("cat_","",colnames(evompc.features))

evompc.features[is.na(evompc.features)] = 0

ranked.features = corfeat.res$feature
f0=lm(data=evompc.features, formula = EVO.FULL ~MPC)
f.start =additive.model(M0 = f0,fulldata = evompc.features)
ivar=1
for( feat in ranked.features[1:10] ){
  if(ivar==1){ M = f.start }
  M = additive.model(M0 = M, addTerm = feat, fulldata = evompc.features)
  ivar=ivar+1
}

f0.fitted.data = get_fit_data(get_XY_data(evompc.features),m = f0)
mycols=c('.resid',colnames(evompc.features)[-1])
add.bestfeat = str_c(best.prop, collapse = ' + ')
fit.final = step(f0,scope = as.formula(paste0('. ~ . +',add.best)), direction = 'forward')

f1=lm(".resid ~ offset(coef(f0)[2] * MPC) + .", data =f0.fitted.data[,mycols])
step(f0, scope = . ~ . , direction = "forward")
