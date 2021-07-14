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
library(here)
#### _Graphics ####
library(ggplot2)
library(ggthemes)
library(ggsci)
library(ggrepel)
library(ggpubr)
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
    ggpubr::grids(color = 'white') +
    stat_density(bw = 0.1, aes_string(x=V)) +
    geom_vline(xintercept = mu, col='white', linetype=1,size=1) + # mean
    geom_segment(data=df_sd,aes(x=xx,xend=xx,y=rep(0,nsd*2),yend=yy), col='white', linetype=2,size=1) + # median
    geom_segment(data=df_sd[1:2,],aes(x=xx,xend=xx,y=yy,yend=c(1,1)), col='red', linetype=2,size=0.25) + # standard deviation (outside)
    geom_errorbarh(data=df_sd,aes(xmin=xx[1],xmax=xx[2],y=0.8,height=0.05),col='red',linetype=1,size=1) + # standard deviation
    ylim(0,1) + scale_y_continuous(position = "left") + xlab("Mean ER (full seq.)") +
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
  B = ggplot(selected,aes(y=!!sym(var),fill=!!sym(name),col=!!sym(name))) +
    geom_beeswarm(aes_string(x=name),size = 1.5,  shape = 21, stroke = 0.2,alpha=0.3,groupOnX = T,dodge.width=0.3) +
    geom_hline(df.mu,mapping = aes(yintercept = MU), col='black',linetype=1,size=1) + # mean
    geom_crossbar(data=df.mu,mapping = aes(x=!!sym(name),y=mu.sample,ymin=sdmin,ymax=sdmax),alpha=0.2,size=1,fatten=1) +
    geom_text(data=df.mu,mapping = aes(label=n,x=!!sym(name),y=Inf,vjust='inward')) +
    xlab(name) + ylab('') + ylim(pop.range) +
    theme(legend.position = 'none',legend.direction = 'vertical')
  B

  B = ggplot(selected,aes(y=!!sym(var),fill=!!sym(name),col=!!sym(name))) +


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
  xydata = get.XY.data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']

  expo = function(alpha=1,beta=0,x){ return(alpha * exp(beta*x)) }
  # model #
  model.name='Exponential'
  f=as.formula(paste0(y,"~ expo(alpha,beta,",x,")"))
  init.params=list(alpha=1,beta=-1)
  m=nls(f,start = init.params,data=xydata$df)
  m.params = get.model.params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model = model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

#### LOADING EVOLUTIONARY DATA (FUNGI & STRAINS) ####
tic("Load data")
FUNGI = readRDS(paste0(here("data"),"fungi-evodata.rds"))
STRAINS = readRDS(paste0(here("data"),"yeast_strains-evodata.rds"))
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
feat = FEAT %>%
  pivot_longer(cols = starts_with('cat_'),
               names_to = c('categories','source','feature'),
               names_pattern="cat_(.+)\\.(.+)\\.(.+)",
               values_to = "value") %>%
  mutate(col_feat = paste0(categories,'.',source,'.',feature)) %>%
  dplyr::filter(!is.na(value)) %>% add_count(col_feat,name='size')
toc()

INPUT = FUNGI
evocols = INPUT %>% dplyr::select(starts_with('EVO.')) %>% colnames
EVO = left_join(INPUT,PROP) %>% ungroup()
evo = left_join(INPUT,prop) %>% ungroup()

#### EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ####
Y = "EVO.FULL" # mean Evolutionary rate (full sequence)
#X = "PPM" # Protein Abundance (log10 ppm)
X = "MPC" # median Molecules Per Cell
INPUT = INPUT %>% dplyr::filter( !is.na( !!sym(X)) & !is.na( !!sym(Y) ))
YY = INPUT[[Y]]
XX = INPUT[[X]]
NXY = sum(complete.cases(XX,YY))
MU = mean_(YY)
RG = range_(YY)
dim(INPUT)
#INPUT[[X]][ INPUT[[X]] < 0 ] = 0
#### _Distribution ####
A=show_density(input = INPUT,var = Y)
samples = c('MF_nucleotide_binding','MF_molecular_function','essential_core','essential_dispensable')
B=show_sample(input=evo,name='property',id='ORF',value=samples,var=Y, pop.mean=mean_(INPUT$EVO.FULL), pop.range=range_(INPUT$EVO.FULL))

#### _Fiting model ####
M1=make_linear_fit(INPUT,only.params = F)
M2=make_logistic_fit(INPUT,only.params = F)
M3=make_poly_fit(INPUT,only.params = F,deg=3)
M4=make_expo_fit(INPUT,only.params = F)
M = bind_rows(M1,M2,M3,M4)

selected.par =c('model','RSS','ESS','TSS')
fit.params = bind_rows(
    make_linear_fit(INPUT,only.params = T)[selected.par],
    make_logistic_fit(INPUT,only.params =T)[selected.par],
    make_poly_fit(INPUT,only.params = T,deg=3)[selected.par],
    make_expo_fit(INPUT,only.params = T)[selected.par]
  )
fit.params

#### _XY-SCATTERPLOT WITH FIT   ####
x1=nearest(0.5,XX,n=1,v=T)
M.x1 = M[ M[X] == x1 ,]
pXY=ggplot(INPUT,aes_string(y=Y,x=X)) +
  geom_point(size=2.5,shape=19,alpha=0.5,color='gray70',stroke=0) +
  stat_density2d(size=0.5,color='gray20') +
  #  geom_text(data=CC.ER, aes(label=sprintf(" r %.2f \n p %.1e \n N %s ",estimate,p.value,n)),
  #            x=Inf,y=Inf,hjust='inward',vjust='inward',size=3) +
  ylab('Mean Evolutionary Rate') + xlab('Median Protein Abundance (log10 mpc)') +
  geom_hline(yintercept = MU, col='black',linetype=1,size=1) + # mean
  geom_line(M,mapping=aes(y=.fitted,col=model),size=1,show.legend = F) +
  geom_text(M.x1 ,mapping=aes(label=model,col=model,y=.fitted-0.2),x=0.8,size=5,check_overlap = T,show.legend = F) +
  stat_smooth(method = 'loess', fullrange=T, span=1.2, se = F, col='gray40',size=1) + annotate('text',y=1.9,x=0.8,label='Loess',col='gray40',size=5)#coord_cartesian(xlim = c(-2.5, 4.5), ylim = c(0, 3.5),expand = F) +
pXY


#### _XY-SCATTERPLOT WITH RESIDUAL   ####
col.samples = PROP %>% dplyr::select(contains(samples)) %>% colnames
f=sprintf("%s + %s",paste0(Y,'~offset(-0.5467*',X,')'),paste0(col.samples,collapse=' + '))
fit0 = lm(formula = paste0(Y,'~',X),data=INPUT)
fit = lm(formula = f,data=EVO)
library(ggiraph)
library(ggiraphExtra)

equation1=function(x){coef(fit)[2]*x+coef(fit)[1]}
equation2=function(x){coef(fit)[2]*x+coef(fit)[1]+coef(fit1)[3]}
evo.samples=evo %>% dplyr::filter(property %in% samples)
test =evo.samples %>% group_by(property) %>%
  broom::augment_columns(x=lm(formula=as.formula(paste0(Y,'~offset(',coef(fit0)[2],"*",X,")")),data=.),data = .) %>%
  mutate(ESS=sum_( (.fitted-mu.y)^2 ), TSS=sum_( (!!sym(Y)-mu.y)^2 ), RSS=TSS-ESS,
         s2=TSS/(nXY-1), s2.y = var.y, RS=sum(.resid),
         RSE=sqrt( (1 / (nXY-2)) * RSS), AIC = AIC(m), BIC=BIC(m), LL = readr::parse_number(as.character(logLik(m))))
ggplot(test,aes_string(X,Y)) +
  geom_smooth(data=INPUT,method='lm',mapping=aes_string(X,Y)) +
  geom_line(aes(y=.fitted)) +
  #geom_smooth(method='lm',mapping=aes(col=property,fill=property),) +
  facet_wrap(~property)

res.sig = left_join(M1,prop, by='ORF') %>%
          bind_rows(M1 %>% mutate(property='all_orthologs'. has_prop=T) ) %>%
          dplyr::filter( property %in% c(samples,'all_orthologs')) %>% mutate(property = as.factor(property))
library(gghighlight)
resi=ggplot(res.sig,aes_string(y=Y,x=X,col='property',fill='property')) +
  geom_hline(yintercept = MU, col='black',linetype=1,size=1) + # mean
  stat_density2d(size=0.5) +
  geom_point(size=1.5,shape=21,stroke=0.1) +
  geom_line(aes(y=.fitted,col=property),col='gray30',size=1) +
  ylab('Mean Evolutionary Rate') + xlab('Median Protein Abundance (log10 mpc)') +
  ylim(RG) +  facet_wrap(~property,nrow = 2) +
  gghighlight(property %in% samples,use_group_by	=F,use_direct_label	=F) +
  geom_smooth(data=left_join(M1,PROP)%>% dplyr::select(X,Y,contains(samples)),
              mapping=aes_string(y=Y,x=X),
              method='lm',
              formula = 'y~x+property',inherit.aes = F)
  #  geom_text(data=CC.ER, aes(label=sprintf(" r %.2f \n p %.1e \n N %s ",estimate,p.value,n)),
  #            x=Inf,y=Inf,hjust='inward',vjust='inward',size=3) +
resi
fit1=lm(EVO.FULL~MPC+property,data=res.sig)
ggPredict(fit1,interactive = T)
#fit1=lm(NTAV~age+sex,data=radial)
ggPredict(fit1,interactive = T)


ggplot(radial,aes(y=NTAV,x=age,color=sex))+geom_point()+
  stat_function(fun=equation1,geom="line",color=scales::hue_pal()(2)[1])+
  stat_function(fun=equation2,geom="line",color=scales::hue_pal()(2)[2])


#### FIGURE ####
library(cowplot)
plot_grid(A,B,pXY,resi,nrow=2,ncol=2,align='h', axis = 'tb')



#### FINAL FIT AND RESIDUALS####
prop.sig = left_join(M1,prop)
colab = c('PPM','MPC','gfp','ms')

FIT =  prop.sig %>%
  group_by(col_prop) %>%
  mutate(  N=n(),
           N_RS_pos = sum(.resid>0),
           N_RS_neg = sum(.resid<0),
           RS     = sum(.resid),
           RS.log = sign(RS) * log(abs(RS)),
           RS_pos = sum(.resid[.resid>0]),
           RS_neg = sum(.resid[.resid<0]),
           RS_avg = mean_(.resid),
           ESS.p  = sum( (.fitted-mfull)^2),
           RSS.p  = TSS - ESS.p ,
           RSS.pc  = RSS.p / TSS,
           ESS.pc  = ESS.p / TSS
  ) %>%
  arrange(desc(RSS.pc)) %>% #,desc(.res_snp.RS.pc)) %>%
  dplyr::select(categories,source,property, col_prop, colab, evocols,
                N,N_RS_pos,N_RS_neg,
                TSS,RSE,LL,ESS.p,ESS.pc,
                starts_with('RS'),.resid,.fitted) %>%
  distinct() %>%
  group_by(categories) %>% mutate(rRS = dense_rank(RSS.pc) , ncat=n())

rank.prop = group_by(FIT, property) %>%
            summarise(n=n(),
                      tss=mean_(TSS),
                      ess=mean_(ESS.pc)*100,
                      rss = mean_(RSS)*100,
                      rss.pc = mean_(RSS.pc)*100,
                      res.sum=mean_(RS),
                      res.avg=mean_(RS_avg)) %>%
            dplyr::filter(n>2) %>%
            arrange(desc(ess),desc(rss))
best.prop = rank.prop %>% dplyr::filter(rss>5) %>% pull(property)

BEST = EVO %>% dplyr::select(EVO.FULL,PPM,contains(best.prop,ignore.case = F))
MBEST = lm(EVO.FULL~PPM+.,data=BEST)
COMPACT = step(MBEST)


# # residuals correlated to features
# feat.rk = left_join(STRAINS,feat) %>%
#   dplyr::filter(!is.na(value) & !is.na(SNP.FULL) & !is.na(EVO.FULL) & !is.na(MPC)) %>%
#   group_by(col_feat) %>%
#   mutate(feat_n = n(), rk.value = dense_rank(value), rk.mpc = dense_rank(MPC))
#
# semifit.res = fit.evo_mpc %>% left_join(feat.rk) %>% dplyr::filter(!is.na(col_feat))
#
# corfeat.res = semifit.res %>%
#   group_by(col_feat) %>%
#   mutate( rho2_evo=scor(.resid.evo,value)$estimate^2 ) %>%
#   group_by(categories) %>% mutate( ncat=n(), rk.rho2_evo=dense_rank(desc(rho2_evo)) ) %>%
#   summarise(categories,source,feature,feat_n, rho2_evo, rk.rho2_evo, ncat) %>%
#   distinct() %>% arrange(rk.rho2_evo)
#
#
#
#
#
#
#
#
#
# #### Analysis of residuals ####
# propsub = prop_evo_mpc %>% dplyr::filter( name  %in% c('essential_core') )
# library(gghighlight)
# MU=mean_(INPUT[[Y]])
# p1a= ggplot(prop, aes_string(y = Y, x = X)) +
#   #geom_ribbon(aes(ymin=.fitted.evo+0.01, ymax=Inf), fill='#BB0033',alpha=0.3)+
#   #geom_ribbon(aes(ymin=-Inf,ymax=.fitted.evo-0.01), fill='#00BB33',alpha=0.3)+
#   geom_abline(slope=0,intercept = MU,size=2,col='gray',linetype=2)+
#   geom_point(fill='lightgray',shape=19,alpha=0.3) +
#   geom_segment(aes(xend = !!sym(X), yend = .fitted.evo), col='gray',linetype='12',size=0.5,alpha=0.5) +
#   gghighlight()
#   # Overlay the protein sharing a particular properties
#   #geom_line(mapping=aes(y=.fitted.evo),size=3) +
#   #geom_segment(data=propsub, aes(xend = !!sym(X), yend = .fitted.evo,color=name),linetype=2,size=1,alpha=0.8) +
#   #geom_point(data=propsub, aes(color=property), shape=19,size=2) + scale_color_simpsons() +
#   #annotate("text",y = 3, x=0.8, label='> expected ',col='#BB0033',vjust='inward',hjust='inward',size=12)+
#   #annotate("text",y = 0, x=0.8, label='< expected',col='#00BB33',vjust='inward',hjust='inward',size=12)+
#   xlab('Protein abundance (log10 mpc)') + ylab('Mean ER (full seq.)') +
#   theme(legend.position = 'top',legend.direction = 'vertical')
#
# plot(p1a)
#
# p1b=     ggplot(data=fit.evo_mpc, aes(y=.resid.evo,x=MPC)) +
#   geom_ribbon(aes(ymin=0.01, ymax=Inf), fill='#BB0033',alpha=0.3)+
#   geom_ribbon(aes(ymin=-Inf,ymax=-0.01), fill='#00BB33',alpha=0.3)+
#   geom_abline(slope=0,intercept = 0,size=2)+
#   geom_point(col='gray30',shape=19,alpha=0.3) +
#   geom_segment(aes(xend = MPC, yend = 0),col='gray',linetype='12',size=0.5,alpha=0.5) +
#   geom_segment(data=propsub,aes(xend = MPC, y=0, yend = .resid.evo,color=property),linetype=2,size=0.5,alpha=0.8) +
#   geom_point(data=propsub, aes(color=property), shape=19,size=2) + scale_color_simpsons()+
#   annotate("text",y = 2, x=0.8, label='> expected',col='#BB0033',vjust='inward',hjust='inward',size=12)+
#   annotate("text",y = -1.2, x=0.8, label='< expected',col='#00BB33',vjust='inward',hjust='inward',size=12)+
#   xlab('Protein abundance (log10 mpc)') + ylab('Residuals (ER)') + theme(legend.position="none")
# plot(p1b)
#
# propsubres = prop_evo_mpc %>% dplyr::filter( property  %in% c('essential_core','essential_dispensable', 'pangenome_rare','pangenome_variable') )
# YMIN = min(0,propsubres$RSS.pc) %>% RoundDownToNearest(roundto = 0.1)
# YMAX = max(0,propsubres$RSS.pc) %>% RoundUpToNearest(roundto = 1)
# p1c =  ggplot(propsubres,
#               aes(x = reorder(str_trunc(property,w=50), -RSS.pc),
#                   y = RSS.pc, fill=property), group=source) +
#   geom_hline(yintercept = 0,size=1,col='gray',linetype=2) +
#   geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=0,ymax=Inf), fill='#BB0033',alpha=0.3)+
#   geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=0), fill='#00BB33',alpha=0.3)+
#   annotate("text",y = Inf, x=Inf, label='> expected',col='#BB0033',vjust='inward',hjust='inward',size=10)+
#   annotate("text",y = -Inf, x=Inf, label='< expected',col='#00BB33',vjust='inward',hjust='inward',size=10)+
#   geom_bar(stat='identity',position='dodge') +
#   geom_text(aes(label=N),vjust='inward',size=12) +
#   theme( axis.text = element_text(size=20)) +
#   ylim(YMIN,YMAX) + scale_fill_simpsons() +  coord_flip()+  theme(legend.position="none",axis.text.y = element_text(angle=45))+
#   ylab("Average residual evolutionary rate") +
#   xlab('Property')
# plot(p1c)
#
# plot_grid(p1a,p1b,p1c,ncol = 3)
# save_plot(plot =  plot_grid(p1a,p1b,p1c,ncol = 3), filename = here("output",'F1.residual-evolutionary-rate.png'),base_width=15,base_height = 5)
