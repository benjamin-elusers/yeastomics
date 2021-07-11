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

make_linear_fit = function(input,    # Input data
                           y=Y, x=X, # X/Y Variables to fit
                           only.params=T){
  YY=input[[y]]
  XX=input[[x]]
  mu = mean_(YY)

  nYY = sum(!is.na(YY))
  nXX = sum(!is.na(XX))
  nXY = sum(!is.na(YY) & !is.na(XX))
  input.noNA= input %>% ungroup() %>% dplyr::filter(!is.na(!!sym(y)) & !is.na(!!sym(x)))

  # model #
  f=as.formula(paste0(Y,"~",X))
  m = lm(formula =f , data=input.noNA)
  m.params = lst(
    m  = m,
    xx = m$model[[x]],
    yy = m$model[[y]],
    yfit = fitted(m), # Y-Fitted
    yres = residuals(m), # Y-Residual
    pfit = coef(m), # Fitted parameters (intercept, PPM)
    ESS = sum( (yfit-mu)^2 ), # Explained variance
    TSS = sum( (yy-mu)^2 ),
    RSS = TSS-ESS, # Deviance (Unexplained variance)
    RSE = sqrt( (1 / (nXY-2)) * RSS ), # Residual standard error
  )

  if(only.params){ return(m.params) }
  # data with model #
  fit.evo_mpc = input.noNA %>%
    broom::augment_columns(x=lm(f,data=.)) %>%
    mutate(nfull=sum_(complete.cases(!!sym(y),!!sym(x))), mfull=mean_(!!sym(y)), vfull=var_(!!sym(y)),
           ESS=sum_( (.fitted-mfull)^2 ), TSS=sum_( (!!sym(y)-mfull)^2 ), RSS=TSS-ESS,
           RSE=sqrt( (1 / (nXY-2)) * RSS), LL = logLik(m), s2=TSS/(nfull-1), RS=sum(.resid))
  return(fit.evo_mpc)
}

make_logistic_fit = function(input,    # Input data
                             y=Y, x=X, # X/Y Variables to fit
                             only.params=T){
  YY=input[[y]]
  XX=input[[x]]
  mu = mean_(YY)

  nYY = sum(!is.na(YY))
  nXX = sum(!is.na(XX))
  nXY = sum(!is.na(YY) & !is.na(XX))
  input.noNA= input %>% ungroup() %>% dplyr::filter(!is.na(!!sym(y)) & !is.na(!!sym(x)))

  # model #
  f=as.formula(paste0(Y,"~","SSlogis(",X,",Asym,xmid,scal)"))
  init.params = getInitial( f,data=input.noNA)
  m = nls(formula = f, start = init.params, data=input.noNA)
  m.params = lst(
    m  = m,
    xx = XX,
    yy = m$m$lhs(),
    yfit = fitted(m), # Y-Fitted
    yres = residuals(m), # Y-Residual
    pfit = coef(m), # Fitted parameters (intercept, PPM)
    ESS = sum( (yfit-mu)^2 ), # Explained variance
    TSS = sum( (yy-mu)^2 ),
    RSS = TSS-ESS, # Deviance (Unexplained variance)
    RSE = sqrt( (1 / (nXY-2)) * RSS ), # Residual standard error
  )

  if(only.params){ return(m.params) }

  # data with model #

  fit.evo_mpc = input.noNA %>%
    broom::augment_columns(x=nls(formula = f, start = init.params, data=.)) %>%
    mutate(nfull=sum_(complete.cases(!!sym(y),!!sym(x))), mfull=mean_(!!sym(y)), vfull=var_(!!sym(y)),
           ESS=sum_( (.fitted-mfull)^2 ), TSS=sum_( (!!sym(y)-mfull)^2 ), RSS=TSS-ESS,
           RSE=sqrt( (1 / (nXY-2)) * RSS), LL = logLik(m), s2=TSS/(nfull-1), RS=sum(.resid))
  return(fit.evo_mpc)
}


make_poly_fit = function(input,    # Input data
                         y=Y, x=X, # X/Y Variables to fit
                         deg=3,
                         only.params=T){
  YY=input[[y]]
  XX=input[[x]]
  mu = mean_(YY)

  nYY = sum(!is.na(YY))
  nXX = sum(!is.na(XX))
  nXY = sum(!is.na(YY) & !is.na(XX))
  input.noNA= input %>% ungroup() %>% dplyr::filter(!is.na(!!sym(y)) & !is.na(!!sym(x)))
  # model #
  f=as.formula(paste0(Y,"~poly(",X,",degree=",deg,")"))
  m = glm(formula =f , data=input.noNA)
  m.params = lst(
    m  = m,
    xx = m$data[[x]],
    yy = m$data[[y]],
    yfit = fitted(m), # Y-Fitted
    yres = residuals(m), # Y-Residual
    pfit = coef(m), # Fitted parameters (intercept, PPM)
    ESS = sum( (yfit-mu)^2 ), # Explained variance
    TSS = sum( (yy-mu)^2 ),
    RSS = TSS-ESS, # Deviance (Unexplained variance)
    RSE = sqrt( (1 / (nXY-2)) * RSS ), # Residual standard error
  )

  if(only.params){ return(m.params) }

  # data with model #

  fit.evo_mpc = input.noNA %>%
    broom::augment_columns(x=glm(formula = f, data=.)) %>%
    mutate(nfull=sum_(complete.cases(!!sym(y),!!sym(x))), mfull=mean_(!!sym(y)), vfull=var_(!!sym(y)),
           ESS=sum_( (.fitted-mfull)^2 ), TSS=sum_( (!!sym(y)-mfull)^2 ), RSS=TSS-ESS,
           RSE=sqrt( (1 / (nXY-2)) * RSS), LL = logLik(m), s2=TSS/(nfull-1), RS=sum(.resid))
  return(fit.evo_mpc)

}


make_expo_fit = function(input,    # Input data
                         y=Y, x=X, # X/Y Variables to fit
                         only.params=T){
  YY=input[[y]]
  XX=input[[x]]
  mu = mean_(YY)

  nYY = sum(!is.na(YY))
  nXX = sum(!is.na(XX))
  nXY = sum(!is.na(YY) & !is.na(XX))
  input.noNA= input %>% ungroup() %>% dplyr::filter(!is.na(!!sym(y)) & !is.na(!!sym(x)))

  expo = function(alpha=1,beta=0,x){
    return(alpha * exp(beta*x))
  }
  # model #
  f=as.formula(paste0(y,"~ expo(alpha,beta,",x,")"))
  init.params=list(alpha=1,beta=-1)
  m=nls(f,start = init.params,data=input.noNA)

  m.params = lst(
    m  = m,
    xx = XX,
    yy = m$m$lhs(),
    yfit = fitted(m), # Y-Fitted
    yres = residuals(m), # Y-Residual
    pfit = coef(m), # Fitted parameters (intercept, PPM)
    ESS = sum( (yfit-mu)^2 ), # Explained variance
    TSS = sum( (yy-mu)^2 ),
    RSS = TSS-ESS, # Deviance (Unexplained variance)
    RSE = sqrt( (1 / (nXY-2)) * RSS ), # Residual standard error
  )

  if(only.params){ return(m.params) }

  # data with model #

  fit.evo_mpc = input.noNA %>%
    broom::augment_columns(x=nls(f,start = init.params,data=.)) %>%
    mutate(nfull=sum_(complete.cases(!!sym(y),!!sym(x))), mfull=mean_(!!sym(y)), vfull=var_(!!sym(y)),
           ESS=sum_( (.fitted-mfull)^2 ), TSS=sum_( (!!sym(y)-mfull)^2 ), RSS=TSS-ESS,
           RSE=sqrt( (1 / (nXY-2)) * RSS), LL = logLik(m), s2=TSS/(nfull-1), RS=sum(.resid))
  return(fit.evo_mpc)
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
XX = INPUT[[Y]]
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
M3=make_poly_fit(INPUT,only.params = F,deg=2)
M4=make_expo_fit(INPUT,only.params = F)

fit.params =
  bind_rows(
    unlist(make_linear_fit(INPUT,only.params = T)[c('RSS','ESS','TSS','dTSS')]),
    unlist(make_logistic_fit(INPUT,only.params =T)[c('RSS','ESS','TSS','dTSS')]),
    unlist(make_poly_fit(INPUT,only.params = T,deg=2)[c('RSS','ESS','TSS','dTSS')]),
    unlist(make_expo_fit(INPUT,only.params = T)[c('RSS','ESS','TSS','dTSS')]),
) %>% bind_cols(model=c('linear','logit','poly','expo'))

fit.params

XX = INPUT[[X]]
YY = INPUT[[Y]]
myfitted = m$pfit[1] / (1 + exp((m$pfit[2]-XX)/m$pfit[3]))
TSS = var_(YY) * (length(YY)-1)
RSS = sum((YY - myfitted)^2)
ESS = sum((myfitted - mean(YY))^2)
ESS+RSS
sqrt(8.8 / length(YY))

head(INPUT$PPM)
y1 = 1.56
MU

TSS.p = (YY - MU)
RSS.p = (YY-myfitted)
ESS.p = (myfitted-MU)
sum( TSS.p - RSS.p-ESS.p )

mlin  = make_linear_fit(INPUT,Y)
myfitted.lin = mlin$pfit[1] + XX * mlin$pfit[2]
TSS.p = (YY - MU)
RSS.p = (YY-myfitted.lin)
ESS.p = (myfitted.lin-MU)
sum( TSS.p - RSS.p-ESS.p )

#### _XY-SCATTERPLOT WITH FIT   ####
pXY=ggplot(INPUT,aes_string(y=Y,x=X)) +
  geom_point(size=1.5,shape=21,alpha=0.7,fill='gray80',stroke=0.1) +
  stat_density2d(size=0.5) +
  #  geom_text(data=CC.ER, aes(label=sprintf(" r %.2f \n p %.1e \n N %s ",estimate,p.value,n)),
  #            x=Inf,y=Inf,hjust='inward',vjust='inward',size=3) +
  ylab('Mean Evolutionary Rate') + xlab('Median Protein Abundance (log10 mpc)') +
  geom_hline(yintercept = MU, col='black',linetype=1,size=1) + # mean
  geom_line(M1,mapping=aes(y=.fitted),col='red',size=1) + annotate('text',y=2,x=-2,label='Linear',col='red',size=5) +
  geom_line(M2,mapping=aes(y=.fitted),col='orange',size=1) + annotate('text',y=1.6,x=-2,label='Logistic',col='orange',size=5) +
  geom_line(M3,mapping=aes(y=.fitted),col='green',size=1) + annotate('text',y=1.2,x=-1.7,label='Polynomial',col='green',size=5) +
  geom_line(M4,mapping=aes(y=.fitted),col='purple',size=1) + annotate('text',y=2.5,x=-1,label='Exponential',col='purple',size=5) +
  ylim(RG)
#coord_cartesian(xlim = c(-2.5, 4.5), ylim = c(0, 3.5),expand = F) +
pXY


#### _XY-SCATTERPLOT WITH RESIDUAL   ####
res.lin = left_join(M1,prop) %>% dplyr::filter( property %in% samples)
library(gghighlight)
resi=ggplot(res.lin,aes_string(y=Y,x=X,fill='property')) +
  geom_point(size=1.5,shape=21,stroke=0.1) +
#  stat_density2d(size=0.5) +
  #  geom_text(data=CC.ER, aes(label=sprintf(" r %.2f \n p %.1e \n N %s ",estimate,p.value,n)),
  #            x=Inf,y=Inf,hjust='inward',vjust='inward',size=3) +
  ylab('Mean Evolutionary Rate') + xlab('Median Protein Abundance (log10 mpc)') +
  geom_hline(yintercept = MU, col='black',linetype=1,size=1) + # mean
  geom_line(mapping=aes(y=.fitted),col='red',size=1)+ #+ annotate('text',y=2,x=-2,label='Linear',col='red',size=5) +
  ylim(rev(RG)) +  facet_wrap(~property,nrow = 2) + gghighlight( property %in% samples) + coord_cartesian(clip='off')


#### FIGURE ####
library(cowplot)
plot_grid(A,B,pXY,resi,nrow=2,ncol=2,align='h', axis = 'tb')



#### FINAL FIT AND RESIDUALS####
res.lin = left_join(M2,prop)
colab = c('PPM','MPC','gfp','ms')

FIT =  res.lin %>%
  group_by(col_prop) %>%
  mutate(  N=n(),
           N_RS_pos = sum(.resid>0),
           N_RS_neg = sum(.resid<0),
           RS     = sum(.resid),
           RS.log = sign(RS) * log(abs(RS)),
           RS_pos = sum(.resid[.resid>0]),
           RS_neg = sum(.resid[.resid<0]),
           RS_avg = mean_(.resid),

           RSS  = sum(.resid^2),
           RSS.pc  = sum(.resid^2)/ TSS,
           RSS_pos = sum(.resid[.resid>0]^2),
           RSS_neg = sum(.resid[.resid<0]^2),

           RSS_delta  = RSS_neg-RSS_pos,
           RSS_ratio  = RSS_neg/RSS_pos,
  ) %>%
  arrange(desc(RSS.pc)) %>% #,desc(.res_snp.RS.pc)) %>%
  dplyr::select(categories,source,property, col_prop, colab, evocols, N,N_RS_pos,N_RS_neg,
                starts_with('RS'),.resid,.fitted) %>%
  distinct() %>%
  group_by(categories) %>% mutate(rRS = dense_rank(RSS.pc) , ncat=n())

rank.prop = group_by(FIT, property) %>% summarise(n=n(),rss = mean_(RSS.pc) ) %>% arrange(desc(rss))
best.prop = rank.prop %>% dplyr::filter(rss >0.05) %>% pull(property)

BEST = EVO %>% dplyr::select(EVO.FULL,PPM,contains(best.prop,ignore.case = F))
MBEST = lm(EVO.FULL~PPM+.,data=BEST)
COMPACT = step(MBEST)

deviance(COMPACT)

var_(YY)
sum_((YY-MU)^2)

# p1a= ggplot(prop, aes_string(y = Y, x = X)) +
#   #geom_ribbon(aes(ymin=.fitted.evo+0.01, ymax=Inf), fill='#BB0033',alpha=0.3)+
#   #geom_ribbon(aes(ymin=-Inf,ymax=.fitted.evo-0.01), fill='#00BB33',alpha=0.3)+
#   geom_abline(slope=0,intercept = MU,size=2,col='gray',linetype=2)+
#   geom_point(fill='lightgray',shape=19,alpha=0.3) +
#   geom_segment(aes(xend = !!sym(X), yend = .fitted.evo), col='gray',linetype='12',size=0.5,alpha=0.5) +
#   gghighlight()
# # Overlay the protein sharing a particular properties
# #geom_line(mapping=aes(y=.fitted.evo),size=3) +
# #geom_segment(data=propsub, aes(xend = !!sym(X), yend = .fitted.evo,color=name),linetype=2,size=1,alpha=0.8) +
# #geom_point(data=propsub, aes(color=property), shape=19,size=2) + scale_color_simpsons() +
# #annotate("text",y = 3, x=0.8, label='> expected ',col='#BB0033',vjust='inward',hjust='inward',size=12)+
# #annotate("text",y = 0, x=0.8, label='< expected',col='#00BB33',vjust='inward',hjust='inward',size=12)+
# xlab('Protein abundance (log10 mpc)') + ylab('Mean ER (full seq.)') +
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
#
# #
# # prop_evo_mpc = fit.evo_mpc %>%
# #   left_join(prop %>% dplyr::filter(has_prop)) %>%
# #   group_by(col_prop) %>%
# #   mutate(  N=n(),
# #            N_RS_pos = sum(.resid.evo>0),
# #            N_RS_neg = sum(.resid.evo<0),
# #            RS     = sum(.resid.evo),
# #            RS.log = sign(RS) * log(abs(RS)),
# #            RS_pos = sum(.resid.evo[.resid.evo>0]),
# #            RS_neg = sum(.resid.evo[.resid.evo<0]),
# #            RS_avg = mean_(.resid.evo),
# #
# #            RSS  = sum(.resid.evo^2),
# #            RSS.pc  = sum(.resid.evo^2)/TSS.evo,
# #            RSS_pos = sum(.resid.evo[.resid.evo>0]^2),
# #            RSS_neg = sum(.resid.evo[.resid.evo<0]^2),
# #
# #            RSS_delta  = RSS_neg-RSS_pos,
# #            RSS_ratio  = RSS_neg/RSS_pos,
# #   ) %>%
# #   arrange(desc(RSS.pc)) %>% #,desc(.res_snp.RS.pc)) %>%
# #   dplyr::select(categories,source,property, col_prop, colab, evocols, N,N_RS_pos,N_RS_neg,
# #                 starts_with('RS'),ends_with('.evo'),) %>%
# #   distinct() %>%
# #   group_by(categories) %>% mutate(rRS = dense_rank(RSS.pc) , ncat=n())
#
#
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
