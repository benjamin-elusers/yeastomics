library(tidyverse)
library(here)
load(here::here('output','evorate-workspace-290422.rdata'))
source(here::here("analysis","function_evorate_fitting.R"))
library(patchwork)
library(scales)

sc_validation = EVOLUTION %>%
             dplyr::filter( !(ORF %in% orf_orthologs) & !is.na(log10.EVO.FULL_R4S) & !is.na(MPC)) %>%
             left_join(PREDICTORS,by='ORF')
#### VARIABLE SELECTION ####
XCOL="MPC"
YCOL="log10.EVO.FULL_R4S"
ZCOL="log10.SNP.FULL_R4S"
IDCOLS = c("UNIPROTKB","SGD","ORF","GNAME","PNAME")
Y_RESID = '.resid'

# SCATTERPLOT EVO RATE - STRUCTURE
{
dom_col ="#3DAA35"
diso_col="#5C99D3"
labx= 'Protein Abundance\n(molecules per cell)'
laby = 'Orthologs Evolutionary Rate (log10)'
labya = 'Disorder Evolutionary Rate (log10)'
labyb= 'Domains Evolutionary Rate (log10)'
F0A =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0, centerY=F,
                    X=XCOL,Y=YCOL, XLAB=labx, YLAB=laby, x_as_exp10 = T,show_cor='left')

F0B =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0, centerY = F,
                    stroke_points=1, col_stat2d='royalblue',
                    col_points = diso_col, col_lm=diso_col,
                    X=XCOL,Y="log10.EVO.DISORDER", XLAB=labx, YLAB=labya, x_as_exp10 = T,show_cor='left')

F0C =   make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0, centerY = F,
                    stroke_points=1, col_stat2d='forestgreen',
                    col_points = dom_col, col_lm=dom_col,
                    X=XCOL,Y="log10.EVO.DOMAINS", XLAB=labx, YLAB=labyb, x_as_exp10 = T,show_cor='left')

F0 = (F0B|F0A|F0C)
ggsave(F0, path = here::here('plots'),
       scale=1.5, width = 12,height=12, bg = 'white',
       filename = 'scatterplot-evo-full_diso_dom.pdf', device = 'pdf')
}

# SCATTERPLOT EVO RATE - TIMESCALES
{
yvals = pretty( range(orthologs[[YCOL]]), n = 8)
yscale = scale_y_continuous(limits=range(yvals), breaks = yvals, expand = c(0,0))

zvals = pretty( range(orthologs[[ZCOL]]), n = 6)
zscale1 = scale_y_continuous(limits=range(zvals), breaks = zvals, expand = c(0,0))
zscale2 = scale_x_continuous(limits=range(zvals), breaks = zvals, expand = c(0,0))

laby = 'Orthologs Evolutionary Rate (log10)'
labz = 'SNP Evolutionary Rate (log10)'
theme_xy = theme_ipsum(base_family = 'Helvetica',base_size=18, axis_text_size = 12) + theme(aspect.ratio = 1)
add_stat_spectral = function(gg){
  gg_2d = gg + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour=NA,size=0.25,alpha=0.15,show.legend = F)+
  scale_fill_distiller(palette = 'Spectral',direction = -1)
  return(gg_2d)
}

F1A =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0,
                      X=XCOL,Y=YCOL, XLAB=labx, YLAB=laby, x_as_exp10 = T,show_cor='left') %>%
  add_stat_spectral(.) + theme_xy + yscale


F1B =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0, ymin=-2.5, centerY = F,
                   X=XCOL,Y=ZCOL, XLAB=labx, YLAB=labz, x_as_exp10 = T,show_cor='left') %>%
  add_stat_spectral(.) + theme_xy + zscale1

F1C =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0,
                    X=ZCOL,Y=YCOL, XLAB=labz, YLAB='', x_as_exp10 = F)  %>%
  add_stat_spectral(.) + theme_xy + yscale + zscale2

c_orthologs = control_var(orthologs, YCOL, XCOL, '_y') %>% control_var(target=ZCOL, XCOL, '_x')

.yvals = pretty( range(c_orthologs[['.resid_y']]), n = 8)
.yscale = scale_y_continuous(limits=range(.yvals), breaks = .yvals, expand = c(0,0))

.zvals = pretty( range(c_orthologs[['.resid_x']]), n = 6)
.zscale = scale_x_continuous(limits=range(.zvals), breaks = .zvals, expand = c(0,0))

F1D =  make_plot_1A(dat=c_orthologs, ANNOT = ANNOTATION, add_outliers = 0, centerY = F,
                    X='.resid_x',Y='.resid_y', XLAB=paste('residual',labz), YLAB=paste('residual',laby),
                    x_as_exp10 = F) %>% add_stat_spectral(.) +
       theme_xy + .yscale + .zscale

F1 = ( (F1A|F1C) / (F1B|F1D) )
ggsave(F1, path = here::here('plots'), family='Helvetica',
       scale=1.5, width = 12,height=12, bg = 'white',
       filename = 'scatterplot-xyz.pdf', device = 'pdf')
ggsave(F1, path = here::here('plots'), scale=1,
       width = 12,height=12, bg = 'white',
       filename = 'scatterplot-xyz.png', device = 'png')
}

#### FILTERING PREDICTORS ####
SC_PREDICTORS = PREDICTORS %>% filter(ORF %in% orthologs$ORF)

## TARGET = ER
fit_ER = fit_m0(orthologs,XCOL,YCOL,SC_PREDICTORS,ZCOL,IDCOLS)
sc_evo_all = preload( here("output","sc-all-lm-evo.rds"), select_variable(fit_ER,response = Y_RESID, raw = T) )#min_ess = 0, min_ess_frac = 0)

# Most explicative variables for evolution
sc_evo = sc_evo_all %>%
  dplyr::filter(pc_ess  > 2 & !variable %in% c(XCOL,YCOL,ZCOL) &
                !str_detect(variable,pattern = 'byrne2005') &
                !str_detect(variable,pattern = 'peter2018.snp_') &
                !str_detect(variable,pattern = 'paxdb.ortho') &
                !str_detect(variable,pattern = 'brar2012') &
                !str_detect(variable,pattern = 'barton2010')
  )

dim(sc_evo)

sc_evo = sc_evo_all %>%
  dplyr::filter(pc_ess  > 0.2 & !variable %in% c(XCOL,YCOL,ZCOL)  )

saveRDS(sc_evo$var,file = here::here('output','sc_features_ess_over_2.rds'))

sc_evo_pred = fit_ER$P[,c(XCOL, YCOL, Y_RESID, ZCOL, sc_evo$variable)]
n_sc_evo = n_distinct(sc_evo$variable)

## TARGET = MPC
fit_EXP = fit_m0(orthologs,YCOL,XCOL,PREDICTORS,ZCOL,IDCOLS,MAX_YCOR = 0.6)
sc_mpc_all = select_variable(fit_EXP,response = XCOL, min_ess = 2, min_ess_frac = 2)
saveRDS(sc_mpc_all,here("output","sc-all-lm-mpc.rds"))
# Most explicative variables for abundance
P_mpc = P_mpc_all %>% dplyr::filter(pc_ess  > 2 )
dim(P_mpc)

mpc_pred = fit_EXP$P[,c(XCOL, YCOL, Y_RESID, ZCOL, P_mpc$variable)]
n_mpc  = n_distinct(P_mpc$variable)

#### CORRELATION BEST FEATURES
cor_evo_pred = cor(PREDICTORS[,P_evo$variable])
image( cor_evo_pred > 0.7 )
#hclust() -> single linkage -> pick best variable per group
# Annotate group by theme
pheatmap::pheatmap( cor_evo_pred  )

cor_mpc_pred = cor(PREDICTORS[,P_mpc$variable])
image( cor_mpc_pred > 0.7 )
#hclust() -> single linkage -> pick best variable per group
# Annotate group by theme
pheatmap::pheatmap( cor_mpc_pred  )

#set.seed(2022)
set.seed(01052022)
#### STEPWISE REGRESSION ####
formula_null_ER = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM_ER = lm(data=sc_evo_pred, formula_null)
decompose_variance(LM_ER)

#formula_best_ER=paste0(YCOL," ~ ",paste0( P_evo$variable,collapse=" + "))
formula_sc_best_ER=reformulate(response = YCOL, termlabels = sc_evo$variable)
m_best_ER = step(object=LM_ER, scope = as.formula(formula_sc_best_ER), direction = 'forward',k=log(n_sc_evo)*2,trace=0)

formula_sc_best_ppm_ER=reformulate(response = YCOL, termlabels = c(XCOL,sc_evo$variable))
#formula_best_ppm_ER=paste0(YCOL," ~ ", XCOL," + ",paste0(P_evo$variable,collapse=" + "))
m_best_ppm_ER = step(object=LM_ER, scope = formula_sc_best_ppm_ER, direction = 'forward', k=log(n_sc_evo+1)*2,trace=0)

formula_sc_best_snp_ER=reformulate(response = YCOL, termlabels = c(ZCOL,sc_evo$variable))
#formula_best_snp_ER=paste0(YCOL," ~ ",ZCOL," + ",paste0( P_evo$variable,collapse=" + "))
m_best_snp_ER = step(object=LM_ER, scope = formula_sc_best_snp_ER, direction = 'forward',k=log(n_sc_evo+1)*2,na.action=na.omit,trace=0)

bind_rows(
  decompose_variance(m_best_ER,T),
  decompose_variance(m_best_ppm_ER,T),
  decompose_variance(m_best_snp_ER,T)
)

dim(sc_evo)

sc_best_lm = lm(reformulate(response = YCOL, termlabels = labels(m_best_ER)),data=sc_evo_pred)
print()
decompose_variance(sc_best_lm,to.df = T)

sc_validation_lm = lm(reformulate(response = YCOL, termlabels = labels(m_best_ER)),data=validation)
sc_validation_lm_ppm = lm(reformulate(response = YCOL, termlabels = labels(m_best_ppm_ER)),data=validation)
sc_validation_lm_snp = lm(reformulate(response = YCOL, termlabels = labels(m_best_snp_ER)),data=validation)

bind_rows(
  decompose_variance(sc_validation_lm,T),
  decompose_variance(sc_validation_lm_ppm,T),
  decompose_variance(sc_validation_lm_snp,T)
)

print(labels(sc_validation_lm))
decompose_variance(sc_validation_lm,to.df = T)

#### ELASTIC NET (LASSO/RIDGE) REGRESSION ####
library(glmnet)
elastic_ER = fit_elastic(P=evo_pred, target = YCOL, add_var = NULL)
elastic_ppm_ER = fit_elastic(P=evo_pred, target = YCOL, add_var =XCOL)
elastic_snp_ER = fit_elastic(P=evo_pred, target = YCOL, add_var =ZCOL)

#### PARTIAL LEAST SQUARES REGRESSION ####
library(pls)
pls_ER <- plsr(formula=as.formula(formula_best_ER), data=evo_pred, scale=TRUE, validation='CV')
pls_ppm_ER <- plsr(formula=as.formula(formula_best_ppm_ER), data=evo_pred, scale=TRUE, validation='CV')
pls_snp_ER <- plsr(formula=as.formula(formula_best_snp_ER), data=evo_pred, scale=TRUE, validation='CV')
decompose_variance(pls_ER)
decompose_variance(pls_ppm_ER)
decompose_variance(pls_snp_ER)

ndims_ER = which.min( RMSEP(pls_ER)$val[estimate = "adjCV", , ]) - 1
ndims_ppm_ER = which.min( RMSEP(pls_ppm_ER)$val[estimate = "adjCV", , ]) - 1
ndims_snp_ER = which.min( RMSEP(pls_snp_ER)$val[estimate = "adjCV", , ]) - 1

pred_ER  = predict(pls_ER, fit_ER$P, ncomp = ndims_ER/6)
pred_val_ER  = predict(pls_ER, validation, ncomp =ndims_ER/6)
cor(fit_ER$P$log10.EVO.FULL_R4S, pred_ER, use='complete')
cor(validation$log10.EVO.FULL_R4S, pred_val_ER, use='complete')

#### SUM OF SQUARES PER VARIABLE ####

library(broom)

# STEPWISE
q1=aov_plot(aov_model(m_best_ER,min_pv = 1e-3), name='best')
q2=aov_plot(aov_model(m_best_ppm_ER), name='ppm')
q3=aov_plot(aov_model(m_best_snp_ER), name='snp')
STEP_PERF = (q1 | q2 | q3)
ggsave(STEP_PERF, path = here::here('plots'),
       scale=1.5, width = 12,height=12, bg = 'white',
       filename = 'performance-stepwise-model.pdf')

# PLS
p1=aov_plot(aov_model(pls_ER), name='best')
p2=aov_plot(aov_model(pls_ppm_ER), name='ppm')
p3=aov_plot(aov_model(pls_snp_ER), name='snp')
PLS_PERF = (p1 / p2 / p3)
ggsave(PLS_PERF, path = here::here('plots'),
       scale=1.5, width = 12,height=12, bg = 'white',
       filename = 'performance-pls-model.pdf')


# ELASTIC
r1 = elastic_dev(elastic)

# USE SET OF VARIABLES FROM STEPWISE REGRESSION FOR ALL ANALYSES
# PREDICT FOR ORTHOLOGS WITHOUT ABUNDANCE/SNP DATA
# FIND VARIABLES PRESENT IN HUMAN AND PREDICT EVORATE IN HUMAN

### CHECK ABUNDANCE
formula_null = reformulate(response=XCOL,termlabels = "1",intercept = T)
LM1.1 = lm(data=mpc_pred, formula_null)
decompose_variance(LM1.1)

formula_MPC=paste0(XCOL," ~ ",paste0(P_mpc$variable,collapse=" + "))
m_mpc = step(object=LM1.1, scope = as.formula(formula_MPC), direction = 'both',k=log(n_mpc))
decompose_variance(m_mpc)
pls_mpc <- plsr(formula=as.formula(formula_MPC), data=mpc_pred, scale=TRUE, validation='CV')
decompose_variance(pls_mpc)
ndims_MPC = which.min( RMSEP(pls_mpc)$val[estimate = "adjCV", , ]) - 1

elastic_mpc = fit_elastic(P=mpc_pred, target = XCOL, add_var = NULL)

nosnp = validation %>% dplyr::filter(is.na(log10.SNP.FULL_R4S))
cor(nosnp$log10.EVO.FULL_R4S,nosnp$MPC)


pred_ER  = predict(m_best_ER, fit_EXP$P, ncomp=1)
pred_val_MPC  = predict(object = LM0, newdata = nosnp)
spearman(nosnp$log10.EVO.FULL_R4S,nosnp$MPC)

pred_val_ER  = predict(object = m_best_ppm_ER, newdata = nosnp)
cor(fit_ER$P$log10.EVO.FULL_R4S, pred_ER, use='complete')^2
cor(nosnp$log10.EVO.FULL_R4S, pred_val_ER, use='complete')^2

test = tidy(aov(pls_mpc)) %>% filter(p.value<1e-5)
ndim_mpc = which.min( RMSEP(pls_mpc)$val[estimate = "adjCV", , ]) - 1

dim(test)
sum(test$sumsq) / decompose_variance(m_mpc,T)$TSS

tidy(aov(m_mpc))

formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM1.0 = lm(data=evo_pred, formula_null)
decompose_variance(LM1.0)
formula_MPC=paste0(YCOL," ~ ",paste0(P_evo$variable,collapse=" + "))
m_mpc.0 = step(object=LM1.0, scope = as.formula(formula_MPC), direction = 'both',k=log(n_evo))
decompose_variance(m_mpc.0)


B=show_sample(input=input.byprop ,name='property',id='ORF',value=samples,var=Y, pop.mean=mu.y, pop.range=rg.y) + ylim(rg.y) + theme(aspect.ratio=1)
samples = c('MF_nucleotide_binding','MF_molecular_function','essential_core','essential_dispensable')
input=fit_ER$P
pop.mean=mean_(fit_ER$P[[Y_RESID]])
pop.range=range_(fit_ER$P[[Y_RESID]])

samples_names = strfind(strings = colnames(input),samples)
property = samples_names[1]
which.sample = fit_ER$P[[property]]
sample = input %>% dplyr::filter(which.sample)
mu.sample = mean_(sample[[Y_RESID]])

df.mu = group_by(input, !!sym(property)) %>%
    summarise(n=n(),mu.sample=mean_(!!sym(Y_RESID)),sd.sample=sd_(!!sym(Y_RESID)),
              sdmin=mu.sample-sd.sample, sdmax=mu.sample+sd.sample) %>%
    mutate(MU = pop.mean)

library(ggbeeswarm)
library(see)
col_mf= '3DA6DE'
col_nb = 'FBBC1D'
B = ggplot(fit_ER$P,aes_string(y=Y_RESID, x=name, fill=name,col=name)) +
  geom_violin(aes_string(x=name),color=NA,alpha=0.9,show.legend = F) +
  geom_pointrange_borderless(data=df.mu,mapping = aes_string(x=name,y='mu.sample',ymin='sdmin',ymax='sdmax'),size=1,fatten=6, position=position_dodge2(width=1)) +
  geom_hline(df.mu,mapping = aes(yintercept = MU), col='black',linetype=1,size=1) + # mean
  geom_text(data=df.mu,mapping = aes(label=n,x=!!sym(name),y=2.5,vjust='inward'),show.legend = F) +
  scale_y_continuous(name='',limits = pop.range, expand = c(0,0)) + ylab('') +
  ggpubr::grids() + xlab('Nucleotide binding') + ylab('Residual Evolutionary Rate (log10)') +
  theme_ipsum() + theme(legend.position = 'none') + scale_fill('TRUE'=)

B


### HUMAN PREDICTION
## best 10 variables from yeast to human
## pick best codon for each AA


