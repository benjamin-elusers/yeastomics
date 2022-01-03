# SETUP DATA -------------------------------------------------------------------
source(here::here("analysis","function_evorate_fitting.R"))
txt_section_break = repchar("-",50)

tic("Load data")
# EXPRESSION DATA
ABUNDANCE = load.abundance()
dim(ABUNDANCE)
# EVOLUTION DATA
CLADE = load.clade()
FUNGI = load.fungi.evo()
STRAINS = load.strains.evo()
EVOLUTION = full_join(STRAINS,CLADE,by=c('ORF'='orf'))
orf_orthologs = EVOLUTION$ORF

# PROTEOME QUALITATIVE AND QUANTITATIVE VARIABLES
PROP = load.properties()
FEAT = load.features() %>% normalize_features()
PREDICTORS = full_join(PROP,FEAT) %>% filter(ORF %in% orf_orthologs)

# ANNOTATION DATA
ANNOTATION=load.annotation()
toc()

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


# Fitting Protein Expression to Evolutionary rate (M0)
na_count = colSums(is.na(PREDICTORS))
na_vars = na_count[na_count>0]

m0 = fit_linear_regression(INPUT=EVOLUTION, X='PPM', Y="log10.EVO.FULL", PREDVAR=PREDICTORS,
                           xcor_max = 0.6,ycor_max = 0.6, min_obs=1 )


# Imputation (replacing NAs)
col_means <- lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), mean, na.rm = TRUE)
col_zeros = lapply(YEASTOMICS %>% dplyr::select(where(is.numeric)), function(x){ return(0) })




library(tidymodels)      # for the recipes package, along with the rest of tidymodels
set.seed(123)

pred_vars = starting(colnames(m0),'cat_')
id_vars = c("ORF","UNIPROT","SGD","PNAME","GNAME","GENENAME","UNIPROTKB")
desc_vars= c("FUNCTION","ROLE","LOC","ORTHO","COMPLEX","OTHER","IS_FUNGI","IS_STRAINS")
m0_vars = c("PPM","log10.EVO.FULL",
            ".fitted",".se.fit",".hat",".sigma",".cooksd",".std.resid",
            "ESS","TSS","RSS","s2","s2.y","RS","RSE","AIC","BIC","LL","model")

lm_mod <- linear_reg() %>% set_engine("lm")
lm_fit <- lm_mod %>%
  fit(reformulate(termlabels = pred_vars, response = '.resid', intercept = T), data = m0)


evo_rec = recipe(x=m0) %>%
  update_role(id_vars, new_role = "ID") %>%
  update_role(desc_vars, new_role = "DESC") %>%
  update_role(m0_vars, new_role = "none") %>%
  update_role(pred_vars, new_role = "predictor") %>%
  update_role(".resid", new_role = "outcome") %>%
  step_dummy(all_nominal_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_center(all_predictors(), -all_outcomes()) %>%
  step_scale(all_predictors(), -all_outcomes())

### Partial Least Squares
#install.packages("pls")
library(pls)
set.seed(1)
#fit PCR model
model <- plsr(formula=reformulate(termlabels = var_names , response = ".resid", intercept = T),data=m0, scale=TRUE)
