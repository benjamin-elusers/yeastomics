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

orthologs = EVOLUTION %>% filter(IS_FUNGI & !is.na(MPC) & !is.na(log10.EVO.FULL_R4S) & !is.na(log10.SNP.FULL_R4S))
orf_evolution = EVOLUTION %>% pull(ORF) %>% unique
orf_orthologs = orthologs %>% pull(ORF) %>% unique

# ANNOTATION DATA
ANNOTATION=load.annotation()

# PROTEOME QUALITATIVE AND QUANTITATIVE VARIABLES
PROP = load.properties()
FEAT = load.features()
PROP_FEAT = inner_join(PROP,FEAT) %>% filter(ORF %in% orf_evolution)
toc()

save_predictors = here::here("output",paste0("PREDICTORS-nomiss.rds"))
save_norm_predictors = here::here("output",paste0("PREDICTORS-nomiss-normalized.rds"))

if(file.exists(save_predictors)){
  PREDICTORS_raw = read_rds(save_predictors)
}else{
  tic("Handling missing values...")
  # HANDLING MISSING VALUES (MANUALLY & AUTOMATICALLY)
  PREDICTORS_raw = PROCESS_MISSING_VALUES(MAT=PROP_FEAT, IDS=orf_orthologs, taxon=4932)
  #missing_var = check_missing_var(PREDICTORS)
  PREDICTORS_raw = PREDICTORS
  write_rds(PREDICTORS_raw,save_predictors)
  toc()
}

if(file.exists(save_norm_predictors)){
  PREDICTORS_norm = read_rds(save_norm_predictors)
}else{
  tic("Normalize variables...")
  # NORMALIZING VARIABLES (MANUALLY & AUTOMATICALLY)
  PREDICTORS_norm = PREDICTORS_raw
  num_pred = PREDICTORS_norm %>% dplyr::select(starts_with("cat_") & where(is.numeric))
  #dplyr::select(-c("cat_transcriptomics.hausse2019.rates_am"),-contains("lahtvee2017.exp_log10_trans"))
  numvar_norm = normalize_features(num_pred,automated=F)
  PREDICTORS_norm[,names(num_pred)] = numvar_norm
  write_rds(PREDICTORS_norm,save_norm_predictors)
  toc()
}

PREDICTORS=PREDICTORS_norm
NPRED = PREDICTORS %>% dplyr::select(starts_with("cat_")) %>% colnames %>% n_distinct
# ANALYZE EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ---------------------
XCOL="MPC"
YCOL="log10.EVO.FULL_R4S"
ZCOL="log10.SNP.FULL_R4S"
IDCOLS = c("UNIPROTKB","SGD","ORF","GNAME","PNAME")

### _FIGURE 1A: EVOLUTION vs EXPRESSION -------------------------------------------

cor_ER = spearman.toplot(orthologs[[XCOL]],orthologs[[YCOL]])

F1A=make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 5,ANNOT = ANNOTATION)

F1A0 = make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 0,ANNOT = ANNOTATION,noplot=T) + xlim(range(orthologs$MPC))

F1A_lm=  make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 0,ANNOT = ANNOTATION,col_points = 'black')  +
    geom_smooth(method='lm',color='limegreen') +
    geom_text(data=cor_ER,aes(label=toshow),x=-Inf,y=-Inf,hjust='inward',vjust='inward')


ggsave(F1A0,filename = "scatterplot-ER-empty.png",device = 'png',scale = 1.2, path = "~/Desktop/")
ggsave(F1A,filename = "scatterplot-ER.png",device = 'png',scale = 1.2, path = "~/Desktop/")
ggsave(F1A_lm,filename = "scatterplot-ER-lm-log10.png",device = 'png',scale = 1.2, path = "~/Desktop/")
ggsave(F1A_lm,filename = "scatterplot-ER-lm-log10.pdf",device = 'pdf',scale = 1.2, path = "~/Desktop/")


rates_col = EVO_AB %>% dplyr::select(contains("r4s"),contains("leisr"),contains("iq_")) %>% colnames()
ER_cor = pivot_longer(EVO_AB,cols = all_of(rates_col), names_to = 'ER_method', values_to = 'rate') %>%
  group_by(ER_method) %>%
  summarise( tmp = spearman.toplot(rate,MPC) ) %>% unnest(tmp) %>% arrange(abs(estimate))
F1A_lm=  make_plot_1A(dat=EVO_AB,X='MPC',Y='fungi_r4s_rate',add_outliers = 0,ANNOT = ANNOTATION,id=c('ORF'))  +
  geom_smooth(method='lm',color='limegreen') +
  geom_smooth(method='loess',color='brown',span = 0.5) +
  geom_text(data=ER_cor %>% dplyr::filter(ER_method=='fungi_r4s_rate'),aes(label=toshow),x=-Inf,y=-Inf,hjust='inward',vjust='inward')


F1.out=make_plot_1A(dat=orthologs,X=XCOL,Y=YCOL,add_outliers = 8,ANNOT = ANNOTATION,centerY = T)
x = ggiraph::girafe(ggobj = F1.out)
x <-  ggiraph::girafe_options(x, ggiraph::opts_hover(css = "fill-opacity:1;fill:orange;stroke:red;") )
x
### _FIGURE 1B: BRANCH LENGTH vs EXPRESSION ---------------------------------------
F1B=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = F) +
  ylab("Schizo. clade\nEvolutionary Rate (log10)") + xlab("Sacch. clade\nEvolutionary Rate (log10)")
F1B.resi=make_plot_1B('schizo','sacch.wgd','ppm',use_residuals = T,force_intercept = T) +
ylab("Schizo. clade\nresidual Evolutionary Rate (log10)") + xlab("Sacch. clade\nresidual Evolutionary Rate (log10)")

### _FIGURE 1C: SNP EVOLUTION vs EXPRESSION ---------------------------------------
F1C=make_plot_1C(EVOLUTION,Y=YCOL,X=ZCOL,control_var = XCOL,use_residuals = F) +
    ylab("Orthologs\nEvolutionary Rate (log10)") + xlab("Strains\nEvolutionary Rate (log10)")
F1C.resi=make_plot_1C(EVOLUTION,Y=YCOL,X=ZCOL,control_var=XCOL,use_residuals = T,force_intercept = T) +
  ylab("Orthologs\nresidual Evolutionary Rate (log10)") + xlab("Strains\nresidual Evolutionary Rate (log10)")

library(patchwork)
FIGURE1 = (F1A_lm | (F1B/F1B.resi) | (F1C/F1C.resi)) +
#  plot_layout(widths = 2) +  plot_annotation(tag_levels = 'A'))
ggsave(FIGURE1,filename = "draft-figure1.png",device = 'png',scale = 1.5, path = "~/Desktop/")
ggsave(FIGURE1,filename = "draft-figure1.pdf",device = 'pdf',scale = 1.5, path = "~/Desktop/")

FIG2=( (F1B/F1B.resi) | (F1C/F1C.resi) )
ggsave(FIG2,filename = "draft-figure2-evo.png",device = 'png',scale = 1, path = "~/Desktop/")
ggsave(FIG2,filename = "draft-figure2-evo.pdf",device = 'pdf',scale = 1, path = "~/Desktop/")


#### VARIABLE SELECTION ####
fit0 = fit_m0(orthologs,XCOL,YCOL,PREDICTORS,ZCOL,IDCOLS)
P_best= select_variable(fit0,response = '.resid', min_ess = 1, min_ess_frac = 1)
best_pred = fit0$P[,c(XCOL, YCOL, y_resid, ZCOL, P_best$variable)]
n_best = n_distinct(P_best$variable)

#### CORRELATION BEST FEATURES
cor_best = cor(PREDICTORS[,P_best$variable])
image( cor_best > 0.5 )
pheatmap::pheatmap( cor_best  )


f1=lm(data=fit0$m0, formula = log10.EVO.FULL_R4S ~ offset(coef(fit0$LM0)[2]*MPC) + cat_functions.go.MF_nucleotide_binding + cat_functions.go.MF_molecular_function)
decompose_variance(f1)
samples=c("cat_functions.go.MF_nucleotide_binding","cat_functions.go.MF_molecular_function")
tibble(var=samples) %>%  group_by(var) %>%
  mutate(fit_lm_one_var(var,'.resid',m0))

#### STEPWISE REGRESSION ####
# RESIDUAL LINEAR REGRESSION WITH NO PREDICTORS
formula_null = reformulate(response=YCOL,termlabels = "1",intercept = T)
LM1 = lm(data=best_pred, formula_null)
LM1_SNP = lm(data=best_pred, formula_null)
decompose_variance(LM1)
decompose_variance(LM1_SNP)

formula_best=paste0(YCOL," ~ ",paste0(best_var,collapse=" + "))
m_best = step(object=LM1_SNP, scope = as.formula(formula_best), direction = 'both',k=log(n_best))

formula_best_ppm=paste0(YCOL," ~ ", XCOL," + ",paste0(best_var,collapse=" + "))
m_best_ppm = step(object=LM1_SNP, scope = as.formula(formula_best_ppm), direction = 'both', k=log(n_best))

formula_best_snp=paste0(YCOL," ~ ",ZCOL," + ",paste0(best_var,collapse=" + "))
m_best_snp = step(object=LM1_SNP, scope = as.formula(formula_best_snp), direction = 'both',k=log(n_best),na.action=na.omit)

formula_best_snp_ppm=paste0(YCOL," ~ ",XCOL," + ",ZCOL," + ",paste0(best_var,collapse=" + "))
m_best_snp_ppm = step(object=LM1_SNP, scope = as.formula(formula_best_snp_ppm), direction = 'both',k=log(n_best),na.action=na.omit)

decompose_variance(fit0$LM0)
decompose_variance(LM_SNP)
decompose_variance(LM0_SNP)

decompose_variance(m_best)
decompose_variance(m_best_ppm)
decompose_variance(m_best_snp_ppm)
decompose_variance(m_best_snp)

m_final = lm(formula(m_best), data=best_pred)
decompose_variance(m_final)
df_final = augment(m_final)

#### ELASTIC NET (LASSO/RIDGE) REGRESSION ####
elastic = fit_elastic(P=best_pred, target = YCOL, add_var = NULL,
                      ftrain = 0.66,alphas = seq(0,1,len=101))

elastic_snp = fit_elastic(P=best_pred, target = YCOL, add_var = c(ZCOL),
                          ftrain = 0.66,alphas = seq(0,1,len=101))

elastic_ppm_snp = fit_elastic(P=best_pred, target = YCOL, add_var = c(XCOL,ZCOL),
                          ftrain = 0.66,alphas = seq(0,1,len=101))

df.glmnet = tibble(nvar=elastic$min_mse$glmnet.fit$df,
                   explained.var=elastic$min_mse$glmnet.fit$dev.ratio,
                   l=elastic$min_mse$glmnet.fit$lambda)
df.best= tibble(nvar=elastic$best$df,explained.var=elastic$best$dev.ratio)

p_perf = ggplot(df.glmnet,aes(y=100*explained.var,x=nvar)) +
  geom_line() + geom_point(size=2) +
  #geom_point(data=df.best,col='red') +
  geom_hline(yintercept = c(60,70,75),linetype='dotted') +
  geom_vline(xintercept = c(10,15,20,40,60,80),linetype='dotted')+
  ylab('Explained variance ESS (%)') + xlab('# of variables') +
  ylim(0,80)
p_perf
#ggsave(p_perf,filename = "draft-figure5-best-elastic-performance.png",device = 'png',scale = 1.5, path = "~/Desktop/")

var_final = str_split_fixed(pattern="[_\\.]",rownames(elastic$best$beta),n=4)
final_model = tibble(type=var_final[,2],source=var_final[,3],varname=str_replace_all(var_final[,4],"[-_\\.]"," "),
                     coeff = elastic$best$beta[,1], imp = caret::varImp(elastic$best,alpha=elastic$min_mse$lambda,lambda = elastic$min_mse$lambda.min)$Overall) %>%
  arrange(imp,coeff)

ggplot(final_model %>% filter(imp>0),aes(y=reorder(varname,imp),fill=type,x=imp)) +
  geom_col(orientation='y',width = 0.5,show.legend=F) +
  geom_text(aes(label=varname,x=0),hjust=0,size=2.5,col='gray20') +
  scale_y_discrete(labels = NULL,'') +
  facet_wrap(~type,scale='free_y')
  #scale_x_log10() + xlab('variable importance (log10)')
  #ypred = predict(elastic$best,s=1e-2,newx=as.matrix(predictors[,best_var]))
  #R2_lowest_mse=cor(predictors[,YCOL],ypred)^2


pp = ggplot(final_model %>% filter(abs(coeff)>0)) +
  geom_bar(aes(y=str_wrap(string=reorder(varname,imp),width=15),fill=type,x=imp),stat='identity',orientation = 'y') +
  labs(title=sprintf('elastic net A=%.3f L=%.3f',elastic$alpha_min_mse,elastic$min_mse$lambda.1se),
       subtitle=sprintf('%s predictors | ESS %.2f%%',elastic$best$df,100*elastic$best$dev.ratio)) +
  scale_x_log10() +
  facet_wrap(~type,scale='free_y',ncol = 2) + theme(text = element_text(size=12), legend.position='bottom')
pp
ggsave(pp,filename = "draft-figure4-elastic-variables.png",device = 'png',scale = 1.5, path = "~/Desktop/")

res = expand.grid(A=seq(0,1,by=0.01),
                  L=10^seq(-5,2, by =0.1))

# test = as_tibble(res) %>%
#   rowwise() %>%
#   mutate( fit=list( glmnet(as.matrix(BEST), RESPONSE, alpha = A, lambda = L, standardize = TRUE) ),
#           ESS_rel = 100*fit$dev.ratio, nvar= fit$df)

# ESS_perc = ggplot(test %>% filter(L<10^2)) +
#   geom_raster(aes(x=A,y=L,fill=ESS_rel ),interpolate = T) +
#   scale_y_log10() + scale_fill_viridis_c(name = "% expl. var.") + theme(legend.position='top') +
#   xlab('Alpha') + ylab('Lambda')
#
# NVAR = ggplot(test %>% filter(L<10^2)) +
#   geom_raster(aes(x=A,y=L,fill=nvar ),interpolate = T) +
#   scale_y_log10() + scale_fill_viridis_c(option = 'A',name = '# Variables') + theme(legend.position='top') +
#   xlab('Alpha') + ylab('Lambda')
#
# ggplot(test %>% filter(ESS_rel>70 )) +
#   geom_raster(aes(x=A,y=L,fill=ESS_rel ),interpolate = F) +
#   scale_y_log10() + scale_fill_viridis_c(name = "% expl. var.") + theme(legend.position='top') +
#   xlab('Alpha') + ylab('Lambda')
# F3 = (ESS_perc | NVAR)
# ggsave(F3,filename = "draft-figure3-elastic.png",device = 'png',scale = 1.5, path = "~/Desktop/")


lasso_vip=caret::varImp(elastic$best,alpha=0.5,lambda=1e-2,scale=T)
lasso_var = lasso_vip[lasso_vip$Overall>0,]
names(lasso_var)=rownames(lasso_vip)[lasso_vip$Overall>0]
barplot(sort(lasso_var),beside=T,las=2)

best_noppm = predictors[,best_var]
best_var_ppm = predictors[,c(XCOL,best_var)]
best_var_ppm_snp =predictors[,c(XCOL,ZCOL,best_var)]


### Partial Least Squares
#install.packages("pls")
library(pls)
library(parallel)
#fit PLSR model
pls.options(parallel = makeCluster(8, type = "PSOCK"))
set.seed(1)
pls_model <- plsr(formula=as.formula(formula_best), data=fit0$m0, scale=TRUE, validation='CV')
stopCluster(pls.options()$parallel)
mpi.exit()

cv = RMSEP(pls_model)
best.dims = which.min(cv$val[estimate = "adjCV", , ]) - 1

coefficients = coef(pls_model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])
barplot(coefficients,las=2)

decompose_variance(pls_model)
coefficients = sort(coefficients, decreasing = TRUE)

#### SUM OF SQUARES PER VARIABLE ####
# PLS
pls_ss = tidy(aov(pls_model))
ggplot(pls_ss) + geom_point(aes(x=sumsq/sum(sumsq),y=term))
# ELASTIC
ela = tidy(elastic$best)
ggplot(ela) + geom_point(aes(x=estimate,y=term))
# STEPwise
lmstep = tidy(aov(m_best))
ggplot(lmstep) + geom_point(aes(x=sumsq/sum(sumsq),y=term))


#save.image(file = here('output','evorate-workspace-290422.rdata'))
