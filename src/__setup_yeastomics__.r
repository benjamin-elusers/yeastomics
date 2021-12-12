# Load yeastomics scripts via URL
yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/"
fx=c("annotation","sequence","phylogenetic","analysis","datalocal","datapub")
fct.r=stringr::str_c("function_",fx,".r")
scripts=file.path(yeastomics,"src",c("utils.r",fct.r))
for( script in scripts){ source(script) }
# Load dependencies packages ---------------------------------------------------
library(xfun)
dep.pkg = c("BiocGenerics","Biobase","S4Vectors","parallel","stats4","IRanges","XVector","hutils","AnnotationDbi")
xfun::pkg_attach2(dep.pkg)
# Load visualization packages ----------------------------------------------------
viz.pkg = c("cowplot","ggplot2","ggtext","ggpubr","ggrepel","gridExtra","plotly")
xfun::pkg_attach2(viz.pkg)
# Load graphical packages ----------------------------------------------------
gr.pkg = c("ggsci","ggthemes","hrbrthemes","RColorBrewer")
#"extrafont","extrafontdb"
xfun::pkg_attach2(gr.pkg)
# Load common packages -----------------------------------------------------
main.pkg  = c("tidyverse","tictoc","hablar","here","see","Biostrings")
xfun::pkg_attach2(main.pkg)

# Coloring schemes
violet= RColorBrewer::brewer.pal(name = 'Paired',n=12)[9:10]
purple = colorRampPalette(colors = violet)(5)

# Visualization themes
options(dplyr.summarise.inform = FALSE, dplyr.width=Inf)
theme_set(theme_bw(base_size=16,base_line_size = 0.5) +theme(plot.background = element_blank(),aspect.ratio=1) )

mytheme = theme(plot.background = element_blank(),
                panel.background =element_blank(),
                panel.grid.minor = element_line(color="#AAAAAA"),
                axis.title = element_markdown(family="Helvetica",colour = 'black',size=18),
                axis.text=element_text(color='black',size=10,family="Helvetica"),
                axis.line = element_line(color='black',size=0.5),
                text=element_text(family='Helvetica'))

# source("https://raw.githubusercontent.com/clauswilke/dviz.supp/master/R/dviz.supp.R")
# source("https://raw.githubusercontent.com/clauswilke/dviz.supp/master/R/themes.R")
# empty_theme <- theme_dviz_open(12, rel_small = 1, rel_large = 1) +
#   theme(
#     axis.line = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.ticks.length = grid::unit(0, "pt")
#   )
# theme_set( theme_modern(base_size = 16, axis.text.angle = 45,legend.position = 'none') + theme(aspect.ratio=1) )