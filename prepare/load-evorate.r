source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)

#### YEAST EVORATE ####
# Get R4S results --------------------------------------------------------------
sc_evo=list()

##### Fungi lineage ------------------------------------------------------------
fungi_dir = "/media/WEXAC_data/FUNGI/"
fungi.rds =  here('output','sc-evorate-fungi.rds')

sc_evo$fungi = preload(fungi.rds,
        load.evorate(alndir=file.path(fungi_dir,'fasta'), id_type = 'ORF', ext.seq = 'fasta', ref = "Saccharomyces_cerevisiae"),
        'get evolutionary rates for fungi...')

##### Cerevisiae isolates ------------------------------------------------------
strains_dir = "/media/WEXAC_data/1011G/"
strains.rds =  here('output','sc-evorate-yk11.rds')

sc_evo$yk11 = preload(strains.rds,
        load.evorate(alndir=file.path(strains_dir,'aln_s288c'), id_type = 'ORF', ext.seq = 'fasta', ref = 'S288C'),
        'get evolutionary rates for 1011 isolated yeast...')

#### HUMAN EVORATE ####
mammals_dir = "/media/WEXAC/MAMMALS/"
phylums = c('mammals','mammals_all','euarchontoglires','laurasiatheria','glires','laurasiatheria.1','primates','carnivora')
alndir = normalizePath(file.path(mammals_dir,paste0('hs_',phylums),'ali')) %>% set_names(phylums) %>% as.list

path_ortho = here::here('output','ens_hs_ortho')
r4s.saved = here('output',paste0('hs-r4s-',phylums,'.rds')) %>% set_names(phylums) %>% as.list

# Get R4S results for each phylum ----------------------------------------------
hs_r4s=list()
##### Euarchontoglires (with human) --------------------------------------------
hs_r4s$euarchontoglires = preload(r4s.saved$euarchontoglires,
        load.evorate(alndir=alndir$euarchontoglires, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for euarchontoglires...')

##### Laurasiatheria (without humans) ------------------------------------------
hs_r4s$laurasiatheria = preload(r4s.saved$laurasiatheria,
        load.evorate(alndir=alndir$laurasiatheria, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for laurasiatheria...')

##### Mammals (selected with >79% orthologs with respect with human) -----------
hs_r4s$mammals = preload(r4s.saved$mammals,
        load.evorate(alndir=alndir$mammals, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for mammals...')

##### All Mammals (from Ensembl) -----------------------------------------------
hs_r4s$mammals_all = preload(r4s.saved$mammals_all,
        load.evorate(alndir=alndir$mammals_all, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for all mammals...')

##### Glires -------------------------------------------------------------------
hs_r4s$glires = preload(r4s.saved$glires,
        load.evorate(alndir=alndir$glires, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for all glires...')

##### Laurasiatheria.1 (sub phylum) --------------------------------------------
hs_r4s$laurasiatheria.1 = preload(r4s.saved$laurasiatheria.1,
        load.evorate(alndir=alndir$laurasiatheria.1, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for all laurasiatheria (sub-phylum)...')

##### Primates -----------------------------------------------------------------
hs_r4s$primates = preload(r4s.saved$primates,
        load.evorate(alndir=alndir$primates, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for all primates...')

##### Carnivora ----------------------------------------------------------------
hs_r4s$carnivora = preload(r4s.saved$carnivora,
        load.evorate(alndir=alndir$carnivora, id_type = 'ENSEMBL', ext.seq = 'mu', ref = NULL),
        'get evolutionary rates for all carnivora...')

sapply(hs_r4s,colnames)
sapply(hs_r4s,dim)

# Evolutionary rate ------------------------------------------------------------
BY = c('id','ref_aa','ref_pos')
HS_R4S = left_join(hs_r4s$mammals_all,hs_r4s$mammals, by=BY, suffix=c('','_mammals')) %>%
  left_join(hs_r4s$euarchontoglires, by=BY, suffix=c('','_euarchontoglires')) %>%
  left_join(hs_r4s$laurasiatheria, by=BY, suffix=c('','_laurasiatheria')) %>%
  left_join(hs_r4s$primates, by=BY, suffix=c('','_primates')) %>%
  left_join(hs_r4s$glires, by=BY, suffix=c('','_glires')) %>%
  left_join(hs_r4s$laurasiatheria.1, by=BY, suffix=c('','_laurasiatheria.1')) %>%
  left_join(hs_r4s$carnivora, by=BY, suffix=c('','_carnivora')) %>%
  dplyr::rename_with(starts_with('r4s_rate'),.fn = str_replace, pattern='r4s_rate',replacement='r4s') %>%
  dplyr::rename_with(starts_with('total'),.fn = str_replace, pattern='total',replacement='nsp') %>%
  ungroup()

saveRDS(HS_R4S, here::here('plots','paper_evo','rate4site_mammals_ensembl.rds'))

colnames(HS_R4S)
head(HS_R4S)

NSP=HS_R4S %>% dplyr::select(starts_with('nsp')) %>% distinct() %>% summarize(across(.fn=max_))

HS_EVO = HS_R4S %>% group_by(id) %>% summarize( across(starts_with('r4s'),mean_) )
colnames(HS_EVO)

##### Correlation to Mammals Evolutionary Rate ---------------------------------
library(ggcorrplot)
COR = cor(HS_EVO[,-1] ,use='pairwise',method = 'spearman')
# corrplot::corrplot(COR,method = 'circle',is.corr = T,diag = F,
#                    outline = 'white',order = 'original', addCoef.col = 'white',
#                    number.digits = 2, tl.srt = 45, number.cex = 0.8,
#                    col=corrplot::COL1('Greys',n=100))

pR4S=ggcorrplot::ggcorrplot(COR,method = 'circle',show.diag = F,
                       show.legend = F, outline.color='white',colors = c('white','gray','gray20'),
                       lab = T, lab_col = 'white', lab_size=3, digits = 2)


p0 = ggplot(HS_EVO,aes(x=r4s)) + xlab('R4S Mammals (All from Ensembl)') + scale_x_log10() + scale_y_log10()

c1=spearman.toplot(HS_EVO$r4s,HS_EVO$r4s_mammals)
p1 = p0 + geom_point(aes(y=r4s_euarchontoglires)) + ylab('R4S Euarchontoglires') +
  geom_smooth(aes(y=r4s_euarchontoglires),method='lm') +
  geom_text(c1,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

c2=spearman.toplot(HS_EVO$r4s,HS_EVO$r4s_laurasiatheria)
p2 = p0 + geom_point(aes(y=r4s_laurasiatheria)) + ylab('R4S Laurasiatheria') +
  geom_smooth(aes(y=r4s_laurasiatheria),method='lm') +
  geom_text(c2,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

c3=spearman.toplot(HS_EVO$r4s,HS_EVO$r4s_primates)
p3 = p0 + geom_point(aes(y=r4s_primates)) + ylab('R4S Primates') +
  geom_smooth(aes(y=r4s_primates),method='lm') +
  geom_text(c3,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

c4=spearman.toplot(HS_EVO$r4s,HS_EVO$r4s_glires)
p4 = p0 + geom_point(aes(y=r4s_glires)) + ylab('R4S Glires') +
  geom_smooth(aes(y=r4s_glires),method='lm') +
  geom_text(c4,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

c5=spearman.toplot(HS_EVO$r4s,HS_EVO$r4s_laurasiatheria.1)
p5 = p0 + geom_point(aes(y=r4s_laurasiatheria.1)) + ylab('R4S sub-Laurasiatheria')+
  geom_smooth(aes(y=r4s_laurasiatheria.1),method='lm') +
  geom_text(c5,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

c6=spearman.toplot(HS_EVO$r4s,HS_EVO$r4s_carnivora)
p6 = p0 + geom_point(aes(y=r4s_carnivora)) + ylab('R4S Carnivora')+
  geom_smooth(aes(y=r4s_carnivora),method='lm') +
  geom_text(c6,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

##### Correlation of Phylum Evolutionary Rate ----------------------------------
CC1 = spearman.toplot(HS_EVO$r4s_euarchontoglires,HS_EVO$r4s_laurasiatheria)
pp1 = ggplot(HS_EVO,aes(x=r4s_euarchontoglires,y=r4s_laurasiatheria)) + geom_point() +
  xlab('R4S Euarchontoglires') + ylab('R4S Laurasiatheria') + scale_x_log10() + scale_y_log10() +
  geom_smooth(method='lm') +
  geom_text(CC1,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

CCA = spearman.toplot(HS_EVO$r4s_primates,HS_EVO$r4s_glires)
pA = ggplot(HS_EVO,aes(x=r4s_primates,y=r4s_glires)) + geom_point() +
  xlab('R4S Primates') + ylab('R4S Glires') + scale_x_log10() + scale_y_log10() +
  geom_smooth(method='lm') +
  geom_text(CCA,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

CCB = spearman.toplot(HS_EVO$r4s_primates,HS_EVO$r4s_laurasiatheria.1)
pB = ggplot(HS_EVO,aes(x=r4s_primates,y=r4s_laurasiatheria.1)) + geom_point() +
  xlab('R4S Primates') + ylab('R4S sub-Laurasiatheria') + scale_x_log10() + scale_y_log10()+
  geom_smooth(method='lm') +
  geom_text(CCB,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

CCC = spearman.toplot(HS_EVO$r4s_primates,HS_EVO$r4s_carnivora)
pC = ggplot(HS_EVO,aes(x=r4s_primates,y=r4s_carnivora)) + geom_point() +
  xlab('R4S Primates') +  ylab('R4S Carnivora') + scale_x_log10() + scale_y_log10()+
  geom_smooth(method='lm') +
  geom_text(CCC,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)


CCX = spearman.toplot(HS_EVO$r4s_glires,HS_EVO$r4s_laurasiatheria.1)
pX = ggplot(HS_EVO,aes(x=r4s_glires,y=r4s_laurasiatheria.1)) + geom_point() +
  xlab('R4S Glires') + ylab('R4S sub-Laurasiatheria') + scale_x_log10() + scale_y_log10() +
  geom_smooth(method='lm') +
  geom_text(CCX,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

CCY = spearman.toplot(HS_EVO$r4s_glires,HS_EVO$r4s_carnivora)
pY = ggplot(HS_EVO,aes(x=r4s_glires,y=r4s_carnivora)) + geom_point() +
  xlab('R4S Glires') + ylab('R4S Carnivora') + scale_x_log10() + scale_y_log10() +
  geom_smooth(method='lm') +
  geom_text(CCY,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

CCZ = spearman.toplot(HS_EVO$r4s_laurasiatheria.1,HS_EVO$r4s_carnivora)
pZ = ggplot(HS_EVO,aes(x=r4s_laurasiatheria.1,y=r4s_carnivora)) + geom_point() +
  xlab('R4S sub-Laurasiatheria') + ylab('R4S Carnivora') + scale_x_log10() + scale_y_log10() +
  geom_smooth(method='lm') +
  geom_text(CCZ,mapping = aes(label=toshow),hjust='inward',vjust='inward',x=-Inf,y=Inf,size=3)

phylopic.url = "http://phylopic.org/assets/images/submissions/"
phylum.uid = c("b8eda501-71d6-4d26-8c9e-731baedd27b2.1024.png", # eutherian
               "f4b6df56-f216-4a4c-9940-4105da8b462e.original.png", # mammals
               "88a07585-846a-405d-9195-c15c010e7443.1024.png", # euarchontoglires
               "865dbf6d-15d6-4940-b0ac-e42dc82e0891.1024.png", # laurasiatheria
               "a2d8a103-daae-4e62-bbf3-c39a4221339b.original.png", # glires.png
               "3717bf88-d959-4bdd-aa69-54e3dccc53b6.original.png", # laurasiatheria.1
               "5248fbfa-3b96-4e34-b962-58008296cec8.original.png", # primates
               "cddf7307-0714-4a16-984a-54a613cb90ae.original.png") %>%  # carnivora
            set_names(phylums)

phylum.png = paste0(phylopic.url,phylum.uid)

# timetree.org (Estimation)
# Eutherian = 98.9MYA
# Boreoeutheria = 94.0MYA
# Euarchotonglires = 87.2MYA
# Laurasiatheria = 81.3MYA
# Glires = 79.6MYA
# Primates = 74.1MYA
# Carnivora = 55.4MYA

# Abundance vs. Evolutionary rate ------------------------------------------

##### Protein abundance in human -----------------------------------------------

hs_paxdb = get.paxdb(tax=9606,abundance = 'integrated',rm.zero=T)
HS_PPM = hs_paxdb %>%
  group_by(protid,id_uniprot) %>%
  summarize( PPM_MIN_ORGAN = min_(ppm_int),
             PPM_MAX_ORGAN = max_(ppm_int),
             PPM_AVG_ORGAN = mean_(ppm_int),
             PPM_MD_ORGAN = median_(ppm_int),
             PPM_geomAVG_ORGAN = geomean(ppm_int)) %>%
  mutate( across(starts_with('PPM_'), log10) ) %>%
  left_join(hs_uni, by=c('id_uniprot'='extid')) %>%
  left_join( pivot_wider(hs_paxdb, id_cols = c(protid,id_uniprot,n_data,n_int),
                         names_from = 'organ', names_prefix = 'PPM_',
                         values_from = 'ppm_int', values_fn = log10) )

HS_ER = left_join(HS_EVO,HS_PPM,by=c('id'='protid'))

ER_COR = cor( HS_ER %>% dplyr::select(starts_with('PPM_')),
              HS_ER %>% dplyr::select(starts_with('r4s_')),
              method= 'spearman',use='pairwise.complete')

library(patchwork)
P1 = (p1 / p2) | (p3 / p4 / p5 / p6)

P2 = (pp1 ) | ( (pA / pB / pC ) | (pX / pY / pZ) )

pER = ggcorrplot::ggcorrplot(ER_COR,method = 'circle',show.diag = F,
                       show.legend = F, outline.color='white',
                       lab = T, lab_col = 'black', lab_size=3, digits = 2)

pR4S
P1
P2
pER
##### plots #####
ggsave(pR4S,filename =  file.path(path_ortho,"correlation-r4s-phylums_to_mammals.pdf"),height=9,width=9)
ggsave(P1, filename = file.path(path_ortho,"scatterplot-r4s-phylums_to_mammals.pdf"),height=12,width=15)
ggsave(pER, filename = file.path(path_ortho,"correlation-r4s_ppm-mammals_phylum.pdf"),height=10,width=15)
ggsave(P2, filename = file.path(path_ortho,"scatterplot-r4s_ppm-mammals_phylum.pdf"),height=10,width=14)

