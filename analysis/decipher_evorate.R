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

# Missing value imputation (replacing NAs)
# col_means <- lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
# col_zeros = lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), function(x){ return(0) })
PREDICTORS[is.na(PREDICTORS)] = 0
na_count = colSums(is.na(PREDICTORS))
na_vars = na_count[na_count>0]

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

# LINEAR REGRESSION of Protein Expression and Evolutionary rate (M0) -----------
m0 = fit_linear_regression(INPUT=EVOLUTION, X='PPM', Y="log10.EVO.FULL", PREDVAR=PREDICTORS.5,
                           xcor_max = 0.6,ycor_max = 0.6, min_obs=1 )

predictors = m0 %>%
            dplyr::select(all_of(c(id_vars,m0_vars)),
                          starts_with('cat'),
                          'PPM', "log10.EVO.FULL",
                           -contains(c("peter2018","byrne2005")))
# LINEAR MODEL M0
LM0 = lm(data=predictors, log10.EVO.FULL ~ PPM )
decompose_variance(LM0)
coef0=coef(LM0)
pred_vars = colnames(predictors) %>% str_extract("cat_.+") %>% setdiff(NA)
x0 = sprintf("offset(%s*%s)",coef0[2],names(coef0)[2])

yeastomics = read_rds("released-dataset/yeastOmics-290921-v0.rds") %>%
              filter(ORF %in% orf_orthologs) %>% ungroup()
yeastomics_var = yeastomics %>% dplyr::select(starts_with('cat_')) %>% colnames
yeastomics$log10.EVO.FULL = log10(yeastomics$EVO.FULL)
yeastomics_model =make_linear_fit(x='PPM', y="log10.EVO.FULL",input = yeastomics,only.params = F)

##### VARIABLE SELECTION #####
lm_single=tibble(variable=yeastomics_var) %>%
  group_by(variable) %>%
  mutate(fit_lm_one_var(variable,'.resid',yeastomics_model)) %>%
  mutate(pc.ess_var = 100*ess_var / tss_var,
         dess =  ess_var-ess0,
         pc.dess = dess/tss0,
         pc0.ess_var = 100*ess_var / tss0,
         pc0.rss_var = 100*rss_var / tss0,
         pc0.tss = 100*tss_var/tss0)

best_lm_single = lm_single %>% filter(pc0.ess_var > 1 & pc0.rss_var>1 & pc0.tss> 1 & abs(pc.dess) > 0.01)
best_var = best_lm_single$variable


yeastomics_lm_long =  yeastomics_model %>%
  pivot_longer(cols = starts_with('cat_'),
               names_to = c('categories','source','property'),
               names_pattern="cat_(.+)\\.(.+)\\.(.+)",
               values_to = "has_prop")
mu.y = mean_(yeastomics$EVO.FULL)

BEST_PROP = yeastomics_lm_long %>%
            mutate(col_prop = paste0(categories,'.',source,'.',property)) %>%
            dplyr::filter(has_prop!=0) %>%
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


#### LASSO REGRESSION ####
library(glmnet)
lambdas <- 10^seq(5, -5, by = -.05)
# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(as.matrix(m0[,best_var]), m0$`.resid`, alpha = 0, lambda = lambdas, standardize = TRUE, nfolds = 5)
lambda_best <- lasso_reg$lambda.min
lasso_model <- glmnet(as.matrix(m0[,best_var]), m0$`.resid`, alpha = 0, lambda = lambda_best, standardize = TRUE)
summary(lasso_model)
lasso_vip=caret::varImp(lasso_model,lambda_best)

lasso_vip %>% filter(Overall>0) %>% arrange(desc(Overall)) %>% dplyr::slice_head(n = 10)
lasso_pred = predict(lasso_model, s = lambda_best, newx = as.matrix(m0[,best_var]))
eval_results(m0$.resid, lasso_pred, m0[,best_var] )

#### STEPWISE REGRESSION ####
formula_best_pred=paste0(".resid ~ ",paste0(best_var,collapse="+"))
m_best = step(object=LM1, scope = as.formula(formula_best_pred), direction = 'both')
decompose_variance(m_best)
summary(m_best)

length( coef(m_best) )


sigvar = M3[-42,] %>%
  arrange(desc(`Sum Sq`)) %>%
  filter(`Pr(>F)`<0.05) %>%
  rename(SS=`Sum Sq`) %>%
  separate(variable, sep='\\.', into=c('categories','source','variable')  ) %>%
  mutate(variable=str_trim(variable), type = ifelse(variable %in% BEST_PROP$property, 'property','feature'))


# RESIDUAL LINEAR REGRESSION WITH NO PREDICTORS
LM1 = lm(data=predictors, .resid ~ 1 )
decompose_variance(LM1)



formula_all_pred=paste0(".resid ~ ",paste0(pred_vars,collapse="+"))

#m1 = step(object=LM1, scope = as.formula(formula_all_pred), direction = 'both')

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



