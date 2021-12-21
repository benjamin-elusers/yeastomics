source(here::here("analysis","function_evorate_fitting.R"))
tic("Load data")
ABUNDANCE = load.abundance()
CLADE = load.clade()
FUNGI = load.fungi.evo()
STRAINS = load.strains.evo()
PROP = load.properties()
FEAT = load.features() #%>% normalize_features()
DESC = read_rds(here('data','uniprot-sgd-annotation.rds'))
BIOFUNC = load.vanleeuwen2016.data()
toc()

dim(ABUNDANCE)
dim(CLADE)
dim(FUNGI)
dim(STRAINS)
dim(PROP)
dim(FEAT)
dim(DESC)
dim(BIOFUNC)

#### EVOLUTIONARY RATE (Y) vs. PROTEIN EXPRESSION (X) ####
INPUT=FUNGI %>% left_join(DESC,by='ORF')

Y = "log10.norm.EVO.FULL" # mean Evolutionary rate (full sequence)
X = "MPC" # median Molecules Per Cell
#X = "PPM" # Protein Abundance (log10 ppm)
XYDATA = get_XY_data(INPUT,x=X,y=Y)
input = XYDATA$df

mu.y = XYDATA$mu['y']
rg.y = range_(XYDATA$YY)
dim(input)
#n.xy = XYDATA$n['xy']
#var.y = XYDATA$var['y']


M1=make_linear_fit(input,only.params = F)
selected.par =c('model','RSS','ESS','TSS')
fit.params = bind_rows(
  make_linear_fit(input,only.params = T)[selected.par],
  #    make_logistic_fit(input,only.params =T)[selected.par],
  #    make_poly_fit(input,only.params = T,deg=3)[selected.par],
  #    make_expo_fit(input,only.params = T)[selected.par]
)
fit.params


Y = "log10.norm.EVO.FULL" # mean Evolutionary rate (full sequence)
OUTLIERS.Y = filter(INPUT, dense_rank(log10.norm.EVO.FULL) < 11 | dense_rank(-log10.norm.EVO.FULL) < 11)
OUTLIERS.X =  filter(INPUT, dense_rank(MPC) < 11 | dense_rank(-MPC) < 11)
#slice.iqr(INPUT,x=MPC, lower = 0.005, upper=0.995,negate = T)

F1=ggplot(INPUT,aes_string(y=Y,x=X)) +
  ggiraph::geom_point_interactive(aes(tooltip=FUNCTION, data_id=ORF),size=2.5,shape=19,alpha=0.5,color='gray70',stroke=0) +
  stat_density2d(size=0.5,color='gray20') +
  geom_hline(yintercept = mean_(INPUT$log10.EVO.FULL), col='red',linetype=2,size=0.5) + # mean
  ylab('mean Evolutionary Rate (log10)') + xlab('Mean Protein Abundance (log10 mpc)') +
  ggpubr::grids() +
  geom_text_repel(data=OUTLIERS.Y, aes(label = GENENAME),max.overlaps = 20,col='blue') +
  geom_text_repel(data=OUTLIERS.X, aes(label = GENENAME),max.overlaps = 20,col='red')

F1
ggiraph::girafe(ggobj = F1)

TMP = left_join(CLADE,INPUT, by=c('uni2'='UNIPROT.x'))
plot(y=TMP$Kc2.log10, x=TMP$ppm2.log10)

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

