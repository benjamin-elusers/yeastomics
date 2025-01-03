# Load yeastomics scripts (either via URL or locally
if(curl::has_internet()){ # If connected to internet
  yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/"
}else{ # else load locally
  yeastomics=here::here("src")
}
fx=c("annotation","alignment","sequence","phylogenetic","analysis","datalocal","datapub")
fct.r=stringr::str_c("function_",fx,".r")
scripts=file.path(yeastomics,"src",c("utils.r",fct.r))
for( script in scripts){ source(script) }
# Load dependencies packages ---------------------------------------------------
library(xfun)
dep.pkg = c("BiocManager","parallel","stats4","hutils","RCurl")
xfun::pkg_attach2(dep.pkg)
# Bioconductor package may be installed once (in particular BiocGenerics seems to reinstall everytime)
bioc.pkg = c("BiocGenerics","Biobase","S4Vectors","XVector","AnnotationDbi","IRanges")
#BiocManager::install(bioc.pkg,update=F)
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
theme_set(
  theme_bw(base_size=16,base_line_size = 0.5) +
    theme(plot.background = element_blank(),aspect.ratio=1)
)

mytheme = theme(plot.background =  element_blank(),
                panel.background = element_blank(),
                panel.grid.minor = element_line(color="#AAAAAA"),
                axis.title = element_markdown(family="Helvetica",colour = 'black',size=18),
                axis.text = element_text(color='black',size=10,family="Helvetica"),
                axis.line = element_line(color='black',linewidth=0.5),
                text = element_text(family='Helvetica'))
# Options
options(dplyr.summarise.inform=FALSE,
        dplyr.width=Inf,
        max.print=2e4,
        timeout = max(600, getOption("timeout")))

# Logger
library(log)
.info  = infoLog()
.error =  errorLog()
.warn  = warningLog()
.succ  = successLog()
.dbg = Logger$new("DEBUG")$
  date()$
  time()$
  hook(crayon::bgWhite)
