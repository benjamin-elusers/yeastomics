# SETUP DATA -------------------------------------------------------------------
source(here::here("analysis","function_evorate_fitting.R"))
library(tidyverse)
txt_section_break = repchar("-",50)

tic("Load data")
ABUNDANCE = load.abundance() # EXPRESSION DATA
CLADE = load.clade()
FUNGI = load.fungi.evo()
STRAINS = load.strains.evo()
EVOLUTION = full_join(STRAINS,CLADE,by=c('ORF'='orf')) %>%
            left_join(ABUNDANCE,by=c('ORF'='orf')) # EVOLUTION DATA

fungi_data= read_rds(here::here("data","evorate-fungi-othologs.rds"))
strains_data=  read_rds(here::here("data","evorate-strains-snp.rds"))
evo = left_join(strains_data,fungi_data,by=c('yk11_id'='fungi_id','yk11_ref_pos'='fungi_ref_pos','yk11_ref_aa'='fungi_ref_aa'))
EVO = evo %>% group_by(yk11_id) %>% summarise(across(yk11_r4s_rate:fungi_leisr_local,mean_))
EVO_AB = left_join(EVO,ABUNDANCE,by=c('yk11_id'='orf')) %>% dplyr::rename(ORF=yk11_id)

orthologs = EVOLUTION %>% filter(IS_FUNGI & !is.na(MPC) & !is.na(PPM) & !is.na(log10.EVO.FULL_R4S))
orf_evolution = EVOLUTION %>% pull(ORF) %>% unique
orf_orthologs = orthologs %>% pull(ORF) %>% unique

# ANNOTATION DATA
ANNOTATION=load.annotation()

# PROTEOME QUALITATIVE AND QUANTITATIVE VARIABLES
PROP = load.properties()
FEAT = load.features() %>% normalize_features() %>% distinct()

#%>% normalize_features() %>% distinct()
PROP_FEAT = inner_join(PROP,FEAT) %>% filter(ORF %in% orf_evolution)
PREDICTORS = PROCESS_MISSING_VALUES(MAT=PROP_FEAT, IDS=orf_orthologs)
missing_var = check_missing_var(PREDICTORS)
toc()

# Add the amino acid next to the codon
codons = sapply(seqinr::SEQINR.UTIL$CODON.AA,toupper) %>% as_tibble %>%
  mutate(L = str_replace(L, pattern = "\\*", replacement = "STOP") )

col_codons = strfind(paste0('sgd.',codons$CODON,"$"),strings = colnames(PREDICTORS))
col_codons_aa = paste0(col_codons ,"_",str_to_title(codons$AA),"_",codons$L)
PREDICTORS = PREDICTORS %>% dplyr::rename_with(.cols=all_of(col_codons), ~col_codons_aa)
NPRED = PREDICTORS %>% dplyr::select(starts_with("cat_")) %>% colnames %>% n_distinct
save_predictors = here::here("output",paste0("PREDICTORS-",NPRED,"vars-missing_processed.rds"))
write_rds(PREDICTORS,save_predictors)
# ANALYZE EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ---------------------


XCOL="MPC"
YCOL="EVO.FULL_R4S"
ZCOL="SNP.FULL_R4S"

### _FIGURE 1A: EVOLUTION vs EXPRESSION -------------------------------------------
cor_ER = spearman.toplot(orthologs[[XCOL]],orthologs[[YCOL]])
F1A=make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 5,ANNOT = ANNOTATION)

F1A0 = make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 0,ANNOT = ANNOTATION,noplot=T) + xlim(range(orthologs$PPM))

F1A_lm=  make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 0,ANNOT = ANNOTATION)  +
    geom_smooth(method='lm',color='limegreen') +
    geom_text(data=cor_ER,aes(label=toshow),x=-Inf,y=-Inf,hjust='inward',vjust='inward')
ggsave(F1A0,filename = "scatterplot-ER-empty.png",device = 'png',scale = 1.2, path = "~/Desktop/")
ggsave(F1A,filename = "scatterplot-ER.png",device = 'png',scale = 1.2, path = "~/Desktop/")
ggsave(F1A_lm,filename = "scatterplot-ER-lm.png",device = 'png',scale = 1.2, path = "~/Desktop/")


rates_col = EVO_AB %>% dplyr::select(contains("r4s"),contains("leisr"),contains("iq_")) %>% colnames()
ER_cor = pivot_longer(EVO_AB,cols = all_of(rates_col), names_to = 'ER_method', values_to = 'rate') %>%
  group_by(ER_method) %>%
  summarise( tmp = spearman.toplot(rate,MPC) ) %>% unnest(tmp) %>% arrange(abs(estimate))
F1A_lm=  make_plot_1A(dat=EVO_AB,X='MPC',Y='fungi_r4s_rate',add_outliers = 0,ANNOT = ANNOTATION,id=c('ORF'))  +
  geom_smooth(method='lm',color='limegreen') +
  geom_smooth(method='loess',color='brown',span = 0.5) +
  geom_text(data=ER_cor %>% dplyr::filter(ER_method=='fungi_r4s_rate'),aes(label=toshow),x=-Inf,y=-Inf,hjust='inward',vjust='inward')


F1.out=make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 5,ANNOT = ANNOTATION)
x = ggiraph::girafe(ggobj = F1.out)
x <-  ggiraph::girafe_options(x, ggiraph::opts_hover(css = "fill-opacity:1;fill:orange;stroke:red;") )
x
### _FIGURE 1B: BRANCH LENGTH vs EXPRESSION ---------------------------------------
F1B=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = F)
F1B.resi=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = T,force_intercept = T)

### _FIGURE 1C: SNP EVOLUTION vs EXPRESSION ---------------------------------------
F1C=make_plot_1C(EVOLUTION,Y=YCOL,X=ZCOL,control_var = XCOL,use_residuals = F)
F1C.resi=make_plot_1C(EVOLUTION,Y=YCOL,X=ZCOL,control_var=XCOL,use_residuals = T,force_intercept = T)

library(patchwork)
FIGURE1 = (F1A_lm | (F1B/F1B.resi) | (F1C/F1C.resi)) + plot_layout(widths = 2) +  plot_annotation(tag_levels = 'A') + theme(axis.text = element_text(size=2))
ggsave(FIGURE1,filename = "draft-figure1.png",device = 'png',scale = 1.5, path = "~/Desktop/")
ggsave(FIGURE1,filename = "draft-figure1.pdf",device = 'pdf',scale = 1.5, path = "~/Desktop/")


id_vars = c("UNIPROTKB","SGD","ORF","GNAME","PNAME")
m0_vars = c(".fitted",".se.fit",".resid",".hat",".sigma",".cooksd",".std.resid","ESS","TSS","RSS","s2","s2.y","RS","RSE","AIC","BIC","LL","model")


# LINEAR REGRESSION of Protein Expression and Evolutionary rate (M0) -----------
m0 = fit_linear_regression(INPUT=orthologs, X=XCOL,Y=YCOL, PREDVAR=PREDICTORS,
                           ADD.VARIABLES =ZCOL,
                           xcor_max = 0.6,ycor_max = 0.6, min_obs=1 ) %>% remove_rare_vars()

predictors = m0 %>%
  dplyr::select(all_of(c(id_vars,m0_vars)), starts_with('cat_'),
                XCOL,YCOL, "SNP.FULL_R4S",
               -contains(c("peter2018","byrne2005","brar2012.TE","paxdb","rna_exp","coRdon","barton2010")) )

# LINEAR MODEL M0
formula_M0=reformulate(response = YCOL, termlabels = XCOL, intercept = T )
LM0 = lm(data=predictors, formula_M0)
df0=decompose_variance(LM0,T)

coef0=coef(LM0)
pred_vars = colnames(predictors) %>% str_extract("cat_.+") %>% setdiff(NA)
x0 = sprintf("offset(%s*%s)",coef0[2],names(coef0)[2])

formula_SNP=reformulate(response = YCOL, termlabels = ZCOL, intercept = T )

LM_SNP = lm(data=predictors,formula_SNP)
decompose_variance(LM_SNP)
coef_snp=coef(LM_SNP)

predictors_snp = predictors %>% filter(!is.na(SNP.FULL_R4S))

formula_M0_SNP =  reformulate(response = YCOL, termlabels =c(XCOL,ZCOL), intercept = T )
LM0_SNP = lm(data=predictors_snp,formula_M0_SNP)
decompose_variance(LM0_SNP)
coef0_snp=coef(LM0_SNP)
pred_vars = colnames(predictors) %>% str_extract("cat_.+") %>% setdiff(NA)
x0_snp = paste0( sprintf("offset(%s*%s)",coef0_snp[2:3],names(coef0_snp)[2:3]), collapse=' + ')


#f1=lm(data=m0, formula = log10.EVO.FULL ~ offset(coef(LM0)[2]*PPM) + cat_functions.go.MF_nucleotide_binding + cat_functions.go.MF_molecular_function)
#decompose_variance(f1)
#samples=c("cat_functions.go.MF_nucleotide_binding","cat_functions.go.MF_molecular_function")
#tibble(var=samples) %>%  group_by(var) %>%
#  mutate(fit_lm_one_var(var,'.resid',m0))

TSS = var( LM0$model[[YCOL]] ) * (nrow(LM0$model)-1)
RSS = deviance(LM0)
ESS = TSS-RSS
mu_y  = mean_(predictors[[YCOL]])
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
formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM1 = lm(data=predictors, formula_null)
decompose_variance(LM1)

formula_best_pred=paste0(YCOL," ~ ", XCOL," + ",paste0(best_var,collapse=" + "))
m_best = step(object=LM0, scope = as.formula(formula_best_pred), direction = 'both')

formula_best_pred_no_ppm=paste0(YCOL," ~ ",paste0(best_var,collapse=" + "))
m_best_no_ppm = step(object=LM1, scope = as.formula(formula_best_pred_no_ppm), direction = 'both')

formula_best_pred_snp=paste0(YCOL," ~ ",XCOL," + ",ZCOL," + ",paste0(best_var,collapse=" + "))
m_best_snp = step(object=LM0_SNP, scope = as.formula(formula_best_pred_snp), direction = 'forward')

decompose_variance(m_best)
decompose_variance(m_best_no_ppm)
decompose_variance(m_best_snp)

#### ELASTIC NET (LASSO/RIDGE) REGRESSION ####
set.seed(42)  # Set seed for reproducibility
n=nrow(predictors)
train_rows = sample(n,size=0.66*n)
target = YCOL
RESPONSE = predictors[[target]]

best_noppm = predictors[,best_var]
best_var_ppm = predictors[,c(XCOL,best_var)]
best_var_ppm_snp =predictors[,c(XCOL,ZCOL,best_var)]

BEST=best_noppm
x.train = BEST[train_rows,]
x.test = BEST[-train_rows,]

y.train = RESPONSE[train_rows]
y.test = RESPONSE[-train_rows]

library(glmnet)
#alpha.fit <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0.5, family="gaussian")
alpha_val = seq(0,1,by=0.01)
alpha_fits = lapply(alpha_val,function(a){
  fit = cv.glmnet(as.matrix(x.train), as.matrix(y.train), type.measure="mse", alpha=a, family="gaussian")
  return(fit)
})

results = sapply(1:101,function(i){
 y.fitted = predict(alpha_fits[[i]],s=alpha_fits[[i]]$lambda.1se, newx=as.matrix(x.test))
 mse <- mean_((y.test - y.fitted)^2)
 return(mse)
})

## View the results
names(results) = paste0("a=",alpha_val)
best_alpha = alpha_fits[[which.min(results)]]
alpha_lowest_mse = alpha_val[which.min(results)]
lambda_lowest_mse = best_alpha$lambda.1se

y.trained = predict(best_alpha,s=best_alpha$lambda.1se,newx=as.matrix(x.train))
cor(y.train,y.trained)^2

y.predicted = predict(best_alpha,s=best_alpha$lambda.1se,newx=as.matrix(x.test))
cor(y.test,y.predicted)^2

best_elastic <- glmnet(as.matrix(best_noppm), RESPONSE, alpha = alpha_lowest_mse, lambda = lambda_lowest_mse, standardize = TRUE)
deviance(best_elastic)
explicative_var = best_elastic$beta[abs(best_elastic$beta[,1])>0,]
npred_lowest_mse=length(explicative_var)
ypred = predict(best_elastic,s=1e-2,newx=as.matrix(best_noppm))
R2_lowest_mse=cor(RESPONSE,ypred)^2

var_final = str_split_fixed(pattern="(_|\\.)",rownames(best_elastic$beta),n=3)
final_model = tibble(type=var_final[,2],var=var_final[,3] ,
                     coeff = best_elastic$beta[,1], imp = caret::varImp(best_elastic,alpha=alpha_lowest_mse,lambda = lambda_lowest_mse)$Overall) %>%
              arrange(imp,coeff)

pp = ggplot(final_model %>% filter(abs(coeff)>0)) +
  geom_bar(aes(y=str_wrap(string=reorder(var,imp),width=15),fill=type,x=imp),stat='identity',orientation = 'y') +
  labs(title=sprintf('elastic net A=%.3f L=%.3f',alpha_lowest_mse,lambda_lowest_mse),
       subtitle=sprintf('%s predictors | ESS %.2f%%',npred_lowest_mse,100*R2_lowest_mse)) +
  scale_x_log10() +
  facet_wrap(~type,scale='free_y',ncol = 2) + theme(text = element_text(size=12), legend.position='bottom')
pp
ggsave(pp,filename = "draft-figure4-elastic-variables.png",device = 'png',scale = 1.5, path = "~/Desktop/")

res = expand.grid(A=seq(0,1,by=0.01),
                  L=10^seq(-5,2, by =0.1))

test = as_tibble(res) %>%
  rowwise() %>%
  mutate( fit=list( glmnet(as.matrix(best_noppm), RESPONSE, alpha = A, lambda = L, standardize = TRUE) ),
          ESS_rel = 100*fit$dev.ratio, nvar= fit$df)

ESS_perc = ggplot(test %>% filter(L<10^2)) +
  geom_raster(aes(x=A,y=L,fill=ESS_rel ),interpolate = T) +
  scale_y_log10() + scale_fill_viridis_c(name = "% expl. var.") + theme(legend.position='top') +
  xlab('Alpha') + ylab('Lambda')

NVAR = ggplot(test %>% filter(L<10^2)) +
  geom_raster(aes(x=A,y=L,fill=nvar ),interpolate = T) +
  scale_y_log10() + scale_fill_viridis_c(option = 'A',name = '# Variables') + theme(legend.position='top') +
  xlab('Alpha') + ylab('Lambda')

ggplot(test %>% filter(ESS_rel>70 )) +
  geom_raster(aes(x=A,y=L,fill=ESS_rel ),interpolate = F) +
  scale_y_log10() + scale_fill_viridis_c(name = "% expl. var.") + theme(legend.position='top') +
  xlab('Alpha') + ylab('Lambda')
F3 = (ESS_perc | NVAR)
ggsave(F3,filename = "draft-figure3-elastic.png",device = 'png',scale = 1.5, path = "~/Desktop/")


lasso_vip=caret::varImp(best_elastic,alpha=0.5,lambda=1e-2,scale=T)
lasso_var = lasso_vip[lasso_vip$Overall>0,]
names(lasso_var)=rownames(lasso_vip)[lasso_vip$Overall>0]
barplot(sort(lasso_var),beside=T,las=2)

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



