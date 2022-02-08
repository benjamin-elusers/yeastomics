source(here::here("analysis","function_evorate_fitting.R"))
library(tidyverse)
txt_section_break = repchar("-",50)
tic("Load data")
ABUNDANCE = load.abundance() # EXPRESSION DATA
CLADE = load.clade()
FUNGI = load.fungi.evo()
STRAINS = load.strains.evo()
EVOLUTION = full_join(STRAINS,CLADE,by=c('ORF'='orf')) # EVOLUTION DATA
# ANNOTATION DATA
ANNOTATION=load.annotation()


orthologs = EVOLUTION %>% filter(IS_FUNGI & !is.na(PPM) & !is.na(log10.EVO.FULL))
orf_evolution = EVOLUTION %>% pull(ORF) %>% unique
orf_orthologs = orthologs %>% pull(ORF) %>% unique
toc()

F1.out=make_plot_1A(dat=orthologs,X='PPM',Y='log10.EVO.FULL',add_outliers = 20,ANNOT = ANNOTATION)
x = ggiraph::girafe(ggobj = F1.out)
x <-  ggiraph::girafe_options(x, ggiraph::opts_hover(css = "fill-opacity:1;fill:orange;stroke:red;") )
x

library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)

ppm_out = get_outliers_boundary(orthologs$PPM,n=100)
evo_out = get_outliers_boundary(orthologs$log10.EVO.FULL,n=100)
outliers_ppm = get_extremes(orthologs,'PPM',n=100) %>%
  dplyr::select(ORF,UNIPROT,PPM,log10.EVO.FULL) %>%
  mutate(rk=dense_rank(PPM),
         protexp = cut(PPM,ppm_out,labels=c('  rare',' normal exp','abundant'),include=T),
         evospeed= cut(log10.EVO.FULL,evo_out,labels=c('slow','normal evo','fast'),include=T))

outliers_evo = get_extremes(orthologs,'log10.EVO.FULL',n=100) %>%
  dplyr::select(ORF,UNIPROT,PPM,log10.EVO.FULL) %>%
  mutate(rk=dense_rank(log10.EVO.FULL),
         protexp = cut(PPM,ppm_out,labels=c('  rare',' normal exp','abundant'),include=T),
         evospeed= cut(log10.EVO.FULL,evo_out,labels=c('slow','normal evo','fast'),include=T))

outliers = coalesce_join(outliers_ppm,outliers_evo, by=c('ORF'='ORF'))
table(outliers$protexp)
table(outliers$evospeed)

lo_hi <- compareCluster(ORF~protexp+evospeed, data=outliers,
                        fun="enrichGO",
                        ont           = 'ALL',
                        keyType       = "ORF",
                        OrgDb         = org.Sc.sgd.db,
                        minGSSize	    = 5,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        pool=T)
enrichplot::dotplot(lo_hi)
dotplot(lo_hi, x="protexp",showCategory=1,split='ONTOLOGY',group=T) +  facet_grid(~evospeed)



kk <- enrichKEGG(gene         = outliers$ORF,organism     = 'sce',pvalueCutoff = 0.05)
mkk <- enrichMKEGG(gene = outliers$ORF,
                   organism = 'sce',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1)
enrichWP(outliers$ORF, organism = "Saccharomyces cerevisiae")