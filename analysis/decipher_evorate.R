# SETUP DATA -------------------------------------------------------------------
source(here::here("analysis","function_evorate_fitting.R"))
library(tidyverse)
txt_section_break = repchar("-",50)

tic("Load data")
ABUNDANCE = load.abundance() # EXPRESSION DATA
CLADE = load.clade()
FUNGI = load.fungi.evo()
STRAINS = load.strains.evo()
EVOLUTION = full_join(STRAINS,CLADE,by=c('ORF'='orf')) # EVOLUTION DATA

orf_evolution = EVOLUTION %>% pull(ORF) %>% unique
orf_orthologs = EVOLUTION %>% filter(IS_FUNGI) %>% pull(ORF) %>% unique

# ANNOTATION DATA
ANNOTATION=load.annotation()

# PROTEOME QUALITATIVE AND QUANTITATIVE VARIABLES
PROP = load.properties()
FEAT = load.features() %>% normalize_features() %>% distinct()
PROP_FEAT = inner_join(PROP,FEAT) %>% filter(ORF %in% orf_evolution)
PREDICTORS = PROCESS_MISSING_VALUES(MAT=PROP_FEAT, IDS=orf_orthologs)
missing_var = check_missing_var(PREDICTORS)
toc()

# ANALYZE EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ---------------------
### _FIGURE 1A: EVOLUTION vs EXPRESSION -------------------------------------------
F1A=make_plot_1A(dat=EVOLUTION,X='PPM',Y='log10.EVO.FULL',add_outliers = 5,ANNOT = ANNOTATION)
x = ggiraph::girafe(ggobj = F1A)
x <-  ggiraph::girafe_options(x, ggiraph::opts_hover(css = "fill-opacity:1;fill:orange;stroke:red;") )
x
### _FIGURE 1B: BRANCH LENGTH vs EXPRESSION ---------------------------------------
F1B=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = F)
F1B.resi=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = T,force_intercept = T) + xlab("") + ylab("")

### _FIGURE 1C: SNP EVOLUTION vs EXPRESSION ---------------------------------------
F1C=make_plot_1C(EVOLUTION,Y='log10.EVO.FULL',X='log10.SNP.FULL','PPM',use_residuals = F)
F1C.resi=make_plot_1C(EVOLUTION,Y='log10.EVO.FULL',X='log10.SNP.FULL','PPM',use_residuals = T,force_intercept = T)  + xlab("") + ylab("")

library(patchwork)
FIGURE1 = (F1A | (F1B/F1B.resi) | (F1C/F1C.resi)) + plot_layout(widths = 1) +  plot_annotation(tag_levels = 'A') + theme(axis.text = element_text(size=2))
ggsave(FIGURE1,filename = "draft-figure1.png",device = 'png',scale = 2, path = "~/Desktop/")
ggsave(FIGURE1,filename = "draft-figure1.pdf",device = 'pdf',scale = 2, path = "~/Desktop/")


id_vars = c("UNIPROTKB","SGD","ORF","GNAME","PNAME")
m0_vars = c(".fitted",".se.fit",".resid",".hat",".sigma",".cooksd",".std.resid","ESS","TSS","RSS","s2","s2.y","RS","RSE","AIC","BIC","LL","model")


# LINEAR REGRESSION of Protein Expression and Evolutionary rate (M0) -----------
m0 = fit_linear_regression(INPUT=EVOLUTION, X='PPM', Y="log10.EVO.FULL", PREDVAR=PREDICTORS,
                           ADD.VARIABLES ="log10.SNP.FULL",
                           xcor_max = 0.6,ycor_max = 0.6, min_obs=1 ) %>% remove_rare_vars()


predictors = m0 %>% dplyr::select(all_of(c(id_vars,m0_vars)),
                                  starts_with('cat_'),
                                  'PPM', "log10.EVO.FULL", "log10.SNP.FULL",
                                  -contains(c("peter2018","byrne2005","brar2012.TE","paxdb","rna_exp","coRdon")))

# LINEAR MODEL M0
LM0 = lm(data=predictors, log10.EVO.FULL ~ PPM )
decompose_variance(LM0)
coef0=coef(LM0)
pred_vars = colnames(predictors) %>% str_extract("cat_.+") %>% setdiff(NA)
x0 = sprintf("offset(%s*%s)",coef0[2],names(coef0)[2])

LM_SNP = lm(data=predictors, log10.EVO.FULL ~ log10.SNP.FULL)
decompose_variance(LM_SNP)
coef_snp=coef(LM_SNP)

LM0_SNP = lm(data=predictors, log10.EVO.FULL ~ PPM + log10.SNP.FULL)
decompose_variance(LM0_SNP)
coef0_snp=coef(LM0_SNP)
pred_vars = colnames(predictors) %>% str_extract("cat_.+") %>% setdiff(NA)
x0_snp = paste0( sprintf("offset(%s*%s)",coef0_snp[2:3],names(coef0_snp)[2:3]), collapse=' + ')

# Missing value imputation (replacing NAs)
# col_means <- lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
# col_zeros = lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), function(x){ return(0) })
#col_mins = lapply(PREDICTORS %>% dplyr::select(missing_var$skim_variable), function(x){ return(min_(x)) })
#PREDICTORS[is.na(PREDICTORS)] = 0

#f1=lm(data=m0, formula = log10.EVO.FULL ~ offset(coef(LM0)[2]*PPM) + cat_functions.go.MF_nucleotide_binding + cat_functions.go.MF_molecular_function)
#decompose_variance(f1)
#samples=c("cat_functions.go.MF_nucleotide_binding","cat_functions.go.MF_molecular_function")
#tibble(var=samples) %>%  group_by(var) %>%
#  mutate(fit_lm_one_var(var,'.resid',m0))

TSS = var( LM0$model$log10.EVO.FULL ) * (nrow(LM0$model)-1)
RSS = deviance(LM0)
ESS = TSS-RSS
##### VARIABLE SELECTION #####
lm_single=tibble(variable=pred_vars,tss=TSS,rss_PPM=RSS,ess_PPM=ESS) %>%
  group_by(variable) %>%
  mutate( fit_lm_one_var(variable,'.resid',predictors) )

best = lm_single %>% ungroup() %>%
        mutate(pc_ess_var = 100*ess_var / tss_var ) %>%
        dplyr::filter(pc_ess_var > 1)
best_var = best$variable

#### STEPWISE REGRESSION ####
# RESIDUAL LINEAR REGRESSION WITH NO PREDICTORS
predictors[is.na(predictors)] = 0
LM1 = lm(data=predictors, log10.EVO.FULL ~ 1 )
decompose_variance(LM1)

formula_best_pred=paste0("log10.EVO.FULL ~ PPM + ",paste0(best_var,collapse=" + "))
m_best = step(object=LM1, scope = as.formula(formula_best_pred), direction = 'both')

formula_best_pred_no_ppm=paste0("log10.EVO.FULL ~ ",paste0(best_var,collapse=" + "))
m_best_no_ppm = step(object=LM1, scope = as.formula(formula_best_pred_no_ppm), direction = 'both')

formula_best_pred_snp=paste0("log10.EVO.FULL ~ PPM + log10.SNP.FULL + ",paste0(best_var,collapse=" + "))
m_best_snp = step(object=LM1, scope = as.formula(formula_best_pred_snp), direction = 'both')

decompose_variance(m_best)
decompose_variance(m_best_no_ppm)
decompose_variance(m_best_snp)

#### LASSO REGRESSION ####
target=predictors$log10.EVO.FULL
best_noppm = predictors[,best_var]
best_var_ppm = predictors[,c('PPM',best_var)]
best_var_ppm_snp =predictors[,c('PPM','log10.SNP.FULL',best_var)]

library(glmnet)
lambdas <- 10^seq(5, -5, by = -.05)
# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(as.matrix(best_var_ppm), target, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
lambda_best <- lasso_reg$lambda.min
lasso_model <- glmnet(as.matrix(best_var_ppm), target, alpha = 1, lambda = lambda_best, standardize = TRUE)
lasso_vip=caret::varImp(lasso_model,lambda_best)
lasso_vip %>% filter(Overall>0) %>% dim
lasso_pred = predict(lasso_model, s = lambda_best, newx = as.matrix(best_var_ppm))
eval_results(target, lasso_pred, best_var_ppm )


## WITH SNP
lasso_reg <- cv.glmnet(as.matrix(best_var_ppm_snp), target, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
lambda_best <- lasso_reg$lambda.min
lasso_model <- glmnet(as.matrix(best_var_ppm_snp), target, alpha = 1, lambda = lambda_best, standardize = TRUE)
lasso_vip=caret::varImp(lasso_model,lambda_best)
lasso_vip %>% filter(Overall>0) %>% dim
lasso_pred = predict(lasso_model, s = lambda_best, newx = as.matrix(best_var_ppm_snp))
eval_results(target, lasso_pred, best_var_ppm_snp )


## NO ABUNDANCE
lasso_reg <- cv.glmnet(as.matrix(best_noppm), target, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
lambda_best <- lasso_reg$lambda.min
lasso_model <- glmnet(as.matrix(best_noppm), target, alpha = 1, lambda = lambda_best, standardize = TRUE)
lasso_vip=caret::varImp(lasso_model,lambda_best)
lasso_vip %>% filter(Overall>0) %>% dim
lasso_pred = predict(lasso_model, s = lambda_best, newx = as.matrix(best_noppm))
eval_results(target, lasso_pred, best_noppm )

### Partial Least Squares
#install.packages("pls")
library(pls)
library(parallel)
#fit PLSR model
pls.options(parallel = makeCluster(8, type = "PSOCK"))
set.seed(1)
model <- plsr(formula=as.formula(formula_best_pred), data=predictors, scale=TRUE, validation='CV')
stopCluster(pls.options()$parallel)
mpi.exit()

cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1

coefficients = coef(model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])
barplot(tail(coefficients, 10))

decompose_variance(model)
coefficients = sort(coefficients, decreasing = TRUE)



