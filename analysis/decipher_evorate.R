# SETUP DATA -------------------------------------------------------------------
source(here::here("analysis","function_evorate_fitting.R"))
library(tidyverse)
txt_section_break = repchar("-",50)

tic("Load data")
# EXPRESSION DATA
ABUNDANCE = load.abundance()
# EVOLUTION DATA
CLADE = load.clade()
FUNGI = load.fungi.evo()
STRAINS = load.strains.evo()
EVOLUTION = full_join(STRAINS,CLADE,by=c('ORF'='orf'))
orf_evolution = EVOLUTION %>% pull(ORF) %>% unique
orf_orthologs = EVOLUTION %>% filter(IS_FUNGI) %>% pull(ORF) %>% unique

# PROTEOME QUALITATIVE AND QUANTITATIVE VARIABLES
PROP = load.properties()
FEAT = load.features() %>% normalize_features() %>% distinct()
PREDICTORS = inner_join(PROP,FEAT) %>% filter(ORF %in% orf_evolution)
missing_var_0 = PREDICTORS %>%
  dplyr::select(where(~!is.logical(.x))) %>%
  skimr::skim(.) %>% filter(complete_rate<1)

# Predictors with missing values must be corrected
#  I. Fix missing observations for certain genes (rows):
#     A) missing codons counts
PREDICTORS.1 = fix_missing_codons(PREDICTORS,col_prefix="cat_transcriptomics.sgd.")
#     B) missing protein length/MW
PREDICTORS.2 = fix_missing_peptide_stats(PREDICTORS.1,
                                         col_len='cat_transcriptomics.sgd.prot_size',
                                         col_mw='cat_biophysics.pepstats.mw',
                                         col_mw_avg='cat_transcriptomics.pepstats.mean_MW',
                                         col_charge='cat_biophysics.pepstats.netcharge',
                                         col_pi='cat_biophysics.pepstats.pI')
#     C) missing centrality values (STRING and INTACT network are treated individually)
PREDICTORS.3 = fix_missing_centrality(PREDICTORS.2,col_prefix="cat_interactions.string.")
PREDICTORS.4 = fix_missing_centrality(PREDICTORS.3,col_prefix="cat_interactions.intact.")
#  II.  Remove unnecessary variables (columns):
#     A) rarely observed data (less than 2 observations)
PREDICTORS.5 = PREDICTORS.4 %>% filter(ORF %in% orf_orthologs) %>% remove_rare_vars()

missing_var = PREDICTORS.5 %>%
  dplyr::select(where(~!is.logical(.x))) %>%
  skimr::skim(.) %>% filter(complete_rate<1)

#     C) missing protein disorder (D2P2/IUP)
#     D) missing transcriptomics in barton 2010

library(mice)
library(VIM)
# mice_plot <- aggr(PR %>%  dplyr::select(where(~!is.logical(.x))),
#                    col=c('navyblue','yellow'), numbers=TRUE, sortVars=TRUE,
#                    labels=names(PREDICTORS.5), cex.axis=.7, gap=3, ylab=c("Missing data","Pattern"))
# imputed_Data <- mice(PREDICTORS.5, m=5, maxit = 50, method = 'pmm', seed = 500)
#sum(is.na(PREDICTORS.5)) / prod(dim(PREDICTORS.5))

# ANNOTATION DATA
ANNOTATION=load.annotation()
toc()


PREDICTORS.5[is.na(PREDICTORS.5)] = 0

# ANALYZE EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ---------------------
### _FIGURE 1A: EVOLUTION vs EXPRESSION -------------------------------------------
F1A=make_plot_1A(dat=EVOLUTION,X='PPM',Y='log10.EVO.FULL',add_outliers = 5,ANNOT = ANNOTATION)
x = ggiraph::girafe(ggobj = F1A)
x <-  ggiraph::girafe_options(x, ggiraph::opts_hover(css = "fill-opacity:1;fill:orange;stroke:red;") )
x

### _FIGURE 1B: BRANCH LENGTH vs EXPRESSION ---------------------------------------
F1B=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = F)
F1B.y0=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = T,force_intercept = T) + xlab("") + ylab("")
F1B.y=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = T,force_intercept = F)  + xlab("") + ylab("")

### _FIGURE 1C: SNP EVOLUTION vs EXPRESSION ---------------------------------------
F1C=make_plot_1C(EVOLUTION,Y='log10.EVO.FULL',X='log10.SNP.FULL','PPM',use_residuals = F)
F1C.y0=make_plot_1C(EVOLUTION,Y='log10.EVO.FULL',X='log10.SNP.FULL','PPM',use_residuals = T,force_intercept = T)  + xlab("") + ylab("")
F1C.y=make_plot_1C(EVOLUTION,Y='log10.EVO.FULL',X='log10.SNP.FULL','PPM',use_residuals = T,force_intercept = F) + xlab("") + ylab("")

library(patchwork)
FIGURE1 = (F1A | (F1B/F1B.y0/F1B.y) | (F1C/F1C.y0/F1C.y)) + plot_layout(widths = 1) +  plot_annotation(tag_levels = 'A') + theme(axis.text = element_text(size=2))
ggsave(FIGURE1,filename = "draft-figure1.png",device = 'png',scale = 2, path = "~/Desktop/")
ggsave(FIGURE1,filename = "draft-figure1.pdf",device = 'pdf',scale = 2, path = "~/Desktop/")



id_vars = c("UNIPROTKB","SGD","ORF","GNAME","PNAME")
m0_vars = c(".fitted",".se.fit",".resid",".hat",".sigma",".cooksd",".std.resid","ESS","TSS","RSS","s2","s2.y","RS","RSE","AIC","BIC","LL","model")
# 1ST REGRESSION of Protein Expression and Evolutionary rate (M0) -----------
m0 = fit_linear_regression(INPUT=EVOLUTION, X='PPM', Y="log10.EVO.FULL", PREDVAR=PREDICTORS.5,
                           xcor_max = 0.6,ycor_max = 0.6, min_obs=1 )


# PROTEIN FEATURES ENGINEERING -------------------------------------------------
predictors = m0 %>% dplyr::select(all_of(c(id_vars,m0_vars)),starts_with('cat'),'PPM', "log10.EVO.FULL",
                                  -contains(c("peter2018","byrne2005")))
LM0 = lm(data=predictors, log10.EVO.FULL ~ PPM )
decompose_variance(LM0)
coef0=coef(LM0)

pred_vars = colnames(predictors) %>% str_extract("cat_.+") %>% setdiff(NA)
x0 = sprintf("offset(%s*%s)",coef0[2],names(coef0)[2])

# Missing value imputation
na_count = colSums(is.na(predictors))
na_vars = na_count[na_count>0]
sum.na(pred_vars)

# Imputation (replacing NAs)
# col_means <- lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
# col_zeros = lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), function(x){ return(0) })

LM1 = lm(data=predictors, .resid ~ 1 )
decompose_variance(LM1)

formula_all_pred=paste0(".resid ~ ",paste0(pred_vars,collapse="+"))
#formula_all_pred=reformulate(termlabels = c(x0,pred_vars), response = "log10.EVO.FULL",intercept = T)
m1 = step(object=LM1, scope = as.formula(formula_all_pred), direction = 'both')
decompose_variance(m1)




library(tidymodels)      # for the recipes package, along with the rest of tidymodels
set.seed(123)

pred_vars = starting(colnames(m0),'cat_')
id_vars = c("ORF","UNIPROT","SGD","PNAME","GNAME","GENENAME","UNIPROTKB")
desc_vars= c("FUNCTION","ROLE","LOC","ORTHO","COMPLEX","OTHER","IS_FUNGI","IS_STRAINS")
m0_vars = c("PPM","log10.EVO.FULL",
            ".fitted",".se.fit",".hat",".sigma",".cooksd",".std.resid",
            "ESS","TSS","RSS","s2","s2.y","RS","RSE","AIC","BIC","LL","model")

# lm_mod <- linear_reg() %>% set_engine("lm")
# lm_fit <- lm_mod %>%
#   fit(reformulate(termlabels = pred_vars, response = '.resid', intercept = T), data = m0)
#

# evo_rec = recipe(x=m0) %>%
#   update_role(id_vars, new_role = "ID") %>%
#   update_role(desc_vars, new_role = "DESC") %>%
#   update_role(m0_vars, new_role = "none") %>%
#   update_role(pred_vars, new_role = "predictor") %>%
#   update_role(".resid", new_role = "outcome") %>%
#   step_unknown(all_nominal_predictors()) %>%
#   step_dummy(all_nominal_predictors()) %>%
#   step_zv(all_predictors()) %>%
#   step_center(all_predictors(), -all_outcomes()) %>%
#   step_scale(all_predictors(), -all_outcomes())

### Partial Least Squares
#install.packages("pls")
library(pls)
library(parallel)
#fit PLSR model
pls.options(parallel = makeCluster(8, type = "PSOCK"))
set.seed(1)
model <- plsr(formula=as.formula(formula_all_pred), data=predictors, scale=TRUE, validation='CV')
stopCluster(pls.options()$parallel)
mpi.exit()

cv = RMSEP(model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1
formula_all_pred
coefficients = coef(model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])
barplot(tail(coefficients, 5))


coefficients = coef(model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
names(coefficients) = TidyLabels(Labels(dat)[-1])
coefficients = sort(coefficients, decreasing = TRUE)



# CHECK EACH VARIABLE
test_var = function(varname,inputdata,verbose=F){
  if(verbose)
    catn(varname)
  formula_var = paste0('.resid ~ ',varname)
  linreg = lm(formula_var, data=inputdata)
  #decompose_variance(linreg)

  not0 = inputdata[,varname]!=0
  N=sum(complete.cases(linreg$model))
  TSS = var( linreg$model[,1] ) * (N-1)
  RSS = deviance(linreg)
  ESS = TSS - RSS

  mu.y = mean_(inputdata$`.resid`[not0])
  .fitted.prop = inputdata$`.fitted`[not0] + mu.y
  .resid.prop = inputdata$`.resid`[not0] + mu.y
  TSS.p  = sum( (inputdata$`.resid`[not0] - mu.y)^2 )
  ESS.p  = sum( (.fitted.prop-mu.y)^2)
  RSS.p = TSS.p - ESS.p

  return(tibble(var=str_split_fixed(pattern='\\.',varname,n = 2)[,2],
                tss0=TSS,rss0=RSS,ess0=ESS, nprot_var = sum(not0),
                tss_var=TSS.p,rss_var=RSS.p,ess_var=ESS.p))
}


res_lm1=tibble(variable=pred_vars) %>%
        group_by(variable) %>%
        mutate(test_var(variable,predictors)) %>%
        mutate(pc.ess_var = 100*ess_var / tss_var,
               pc0.ess_var = 100*ess_var / tss0)

best_var = res_lm1$variable[res_lm1$pc0.ess_var>0.5]
formula_best_pred=paste0(".resid ~ ",paste0(best_var,collapse="+"))
m_best = step(object=LM1, scope = as.formula(formula_best_pred), direction = 'both')
decompose_variance(m_best)
