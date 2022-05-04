library(tidyverse)
library(here)
load(here::here('output','evorate-workspace-290422.rdata'))
source(here::here("analysis","function_evorate_fitting.R"))
library(patchwork)
library(scales)


#### VARIABLE SELECTION ####
XCOL="MPC"
YCOL="log10.EVO.FULL_R4S"
ZCOL="log10.SNP.FULL_R4S"
IDCOLS = c("UNIPROTKB","SGD","ORF","GNAME","PNAME")
Y_RESID = '.resid'

laby = 'Evolutionary Rate\nOrthologs (log10)'
labz = 'Evolutionary Rate\nSNP (log10)'
labx= 'Protein Abundance\n(molecules per cell)'
F1A =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0,
                      X=XCOL,Y=YCOL, XLAB=labx, YLAB=laby, x_as_exp10 = T,show_cor='left')
F1B =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0,
                        X=ZCOL,Y=YCOL, XLAB=labz, YLAB=laby, x_as_exp10 = F)
F1C =  make_plot_1A(dat=orthologs,ANNOT = ANNOTATION, add_outliers = 0,
                   X=XCOL,Y=ZCOL, XLAB=labx, YLAB=labz, x_as_exp10 = T,show_cor='left')

F1 = (F1A|F1B|F1C)
ggsave(F1, path = here::here('plots'),
       scale=1.5, width = 12,height=12, bg = 'white',
       filename = 'scatterplot-xyz.pdf')


decompose_variance(LM_SNP)
decompose_variance(LM0)
decompose_variance(LM0_SNP)

#### FILTERING PREDICTORS ####
fit_ER = fit_m0(orthologs,XCOL,YCOL,PREDICTORS,ZCOL,IDCOLS)
P_evo = select_variable(fit_ER,response = Y_RESID, min_ess = 2, min_ess_frac = 2)
evo_pred = fit_ER$P[,c(XCOL, YCOL, y_resid, ZCOL, P_evo$variable)]
n_evo = n_distinct(P_evo$variable)

#### CORRELATION BEST FEATURES
cor_evo_pred = cor(PREDICTORS[,P_evo$variable])
image( cor_evo_pred > 0.7 )
#hclust() -> single linkage -> pick best variable per group
# Annotate group by theme
pheatmap::pheatmap( cor_evo_pred  )

#set.seed(2022)
set.seed(01052022)
#### STEPWISE REGRESSION ####
formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM1 = lm(data=evo_pred, formula_null)
decompose_variance(LM1)

formula_best=paste0(YCOL," ~ ",paste0( P_evo$variable,collapse=" + "))
m_best = step(object=LM1, scope = as.formula(formula_best), direction = 'both',k=log(n_best))
formula_best_ppm=paste0(YCOL," ~ ", XCOL," + ",paste0(P_evo$variable,collapse=" + "))
m_best_ppm = step(object=LM1, scope = as.formula(formula_best_ppm), direction = 'both', k=log(n_best))
formula_best_snp=paste0(YCOL," ~ ",ZCOL," + ",paste0( P_evo$variable,collapse=" + "))
m_best_snp = step(object=LM1, scope = as.formula(formula_best_snp), direction = 'both',k=log(n_best),na.action=na.omit)

decompose_variance(m_best)
decompose_variance(m_best_ppm)
decompose_variance(m_best_snp)

#### ELASTIC NET (LASSO/RIDGE) REGRESSION ####
library(glmnet)
elastic = fit_elastic(P=evo_pred, target = YCOL, add_var = NULL)
elastic_ppm = fit_elastic(P=evo_pred, target = YCOL, add_var =XCOL)
elastic_snp = fit_elastic(P=evo_pred, target = YCOL, add_var =ZCOL)

#### PARTIAL LEAST SQUARES REGRESSION ####
library(pls)
pls <- plsr(formula=as.formula(formula_best), data=evo_pred, scale=TRUE, validation='CV')
pls_ppm <- plsr(formula=as.formula(formula_best_ppm), data=evo_pred, scale=TRUE, validation='CV')
pls_snp <- plsr(formula=as.formula(formula_best_snp), data=evo_pred, scale=TRUE, validation='CV')
decompose_variance(pls)
decompose_variance(pls_ppm)
decompose_variance(pls_snp)


ndims = which.min( RMSEP(pls)$val[estimate = "adjCV", , ]) - 1
ndims_ppm = which.min( RMSEP(pls_ppm)$val[estimate = "adjCV", , ]) - 1
ndims_snp = which.min( RMSEP(pls_snp)$val[estimate = "adjCV", , ]) - 1

#### SUM OF SQUARES PER VARIABLE ####
library(broom)

# PLS
p1=aov_plot(aov_model(pls), name='best')
p2=aov_plot(aov_model(pls_ppm), name='ppm')
p3=aov_plot(aov_model(pls_snp), name='snp')
PLS_PERF = (p1 / p2 / p3)
ggsave(PLS_PERF, path = here::here('plots'),
       scale=1.5, width = 12,height=12, bg = 'white',
       filename = 'performance-pls-model.pdf')

# STEPWISE
q1=aov_plot(aov_model(m_best,min_pv = 1e-3), name='best')
q2=aov_plot(aov_model(m_best_ppm), name='ppm')
q3=aov_plot(aov_model(m_best_snp), name='snp')
STEP_PERF = (q1 / q2 / q3)
ggsave(STEP_PERF, path = here::here('plots'),
       scale=1.5, width = 12,height=12, bg = 'white',
       filename = 'performance-stepwise-model.pdf')

# ELASTIC
r1 = elastic_dev(elastic)

# USE SET OF VARIABLES FROM STEPWISE REGRESSION FOR ALL ANALYSES
# PREDICT FOR ORTHOLOGS WITHOUT ABUNDANCE/SNP DATA
# FIND VARIABLES PRESENT IN HUMAN AND PREDICT EVORATE IN HUMAN

### CHECK ABUNDANCE
YCOL = "MPC"
fit_EXP = fit_m0(orthologs,XCOL,YCOL,PREDICTORS,ZCOL,IDCOLS,MAX_YCOR = 0.6)
P_mpc.0 = select_variable(fit0,response = YCOL, min_ess = 1, min_ess_frac = 1)
P_mpc.1 = select_variable(fit0.1,response = YCOL, min_ess = 1, min_ess_frac = 1)

mpc_pred = fit0.1$P[,c(XCOL, YCOL, ZCOL, P_mpc$variable)]
n_mpc  = n_distinct(P_mpc$variable)

formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM1.1 = lm(data=mpc_pred, formula_null)
decompose_variance(LM1.1)

formula_MPC=paste0(YCOL," ~ ",paste0(P_mpc$variable,collapse=" + "))
m_mpc = step(object=LM1.1, scope = as.formula(formula_MPC), direction = 'both',k=log(n_best))
decompose_variance(m_mpc)

elastic_mpc = fit_elastic(P=mpc_pred, target = YCOL, add_var = NULL)
pls_mpc <- plsr(formula=as.formula(formula_MPC), data=mpc_pred, scale=TRUE, validation='CV')
decompose_variance(pls_mpc)

test = tidy(aov(pls_mpc)) %>% filter(p.value<1e-5)
dim(test)
sum(test$sumsq) / decompose_variance(m_mpc,T)$TSS


tidy(aov(m_mpc))

formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM1.0 = lm(data=evo_pred, formula_null)
decompose_variance(LM1.0)
formula_MPC=paste0(YCOL," ~ ",paste0(P_evo$variable,collapse=" + "))
m_mpc = step(object=LM1, scope = as.formula(formula_MPC), direction = 'both',k=log(n_best))
decompose_variance(m_mpc)
