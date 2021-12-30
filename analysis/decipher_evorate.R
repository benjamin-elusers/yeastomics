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

# PROTEOME QUALITATIVE AND QUANTITATIVE VARIABLES
PROP = load.properties()
FEAT = load.features() %>% normalize_features()
PREDICTORS = full_join(PROP,FEAT)

# ANNOTATION DATA
SGD_DESC = read_rds(here('data','uniprot-sgd-annotation.rds'))
UNI_FEAT = read_rds(here('data','uniprot-features.rds'))
BIOFUNC = load.vanleeuwen2016.data()
dim(SGD_DESC)
dim(UNI_FEAT)
dim(BIOFUNC)
ANNOTATION=full_join(SGD_DESC,BIOFUNC,by="ORF") %>% full_join(UNI_FEAT,by='SGD')
toc()

# ANALYZE EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ---------------------
m0 = fit_linear_regression(EVOLUTION, X='PPM', Y="log10.EVO.FULL", PREDICTORS, 0.6, 0.8) %>% left_join(SGD_DESC)
### _FIGURE 1A: EVOLUTION vs EXPRESSION -------------------------------------------
F1A=make_plot_1A(dat=EVOLUTION,X='PPM',Y='log10.EVO.FULL',add_outliers = 5)
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

TMP = left_join(CLADE,INPUT, by=c('uni2'='UNIPROT.x'))
m=lm(data=TMP,Kc2.log10~ppm2.log10)
decompose_variance(m)
dim(TMP)

spearman(TMP[[Y]], TMP$MPC)

INPUT %>%
  slice.iqr
pXY.fit = pXY+
  geom_line(M,mapping=aes(y=.fitted,col=model),size=1,show.legend = F) +
  geom_text(M.x1 ,mapping=aes(label=model,col=model,y=.fitted-0.2),x=0.8,size=5,check_overlap = T,show.legend = F) +
  stat_smooth(method = 'loess', fullrange=T, span=1.2, se = F, col='gray40',size=1) + annotate('text',y=1.9,x=0.8,label='Loess',col='gray40',size=5)
#coord_cartesian(xlim = c(-2.5, 4.5), ylim = c(0, 3.5),expand = F) +
pXY.fit


# Filter variables (abs cor MPC<0.4)
# Make M0 and M3
# Refine variables selections





