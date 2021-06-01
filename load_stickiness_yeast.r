yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/"
source(file.path(yeastomics,"src/utils.r"))
source(file.path(yeastomics,"src/function_annotation.r"))
source(file.path(yeastomics,"src/function_sequence.r"))
source(file.path(yeastomics,"src/function_phylogenetic.r"))
source(file.path(yeastomics,"src/function_analysis.r"))
source(file.path(yeastomics,"src/function_datalocal.r"))
source(file.path(yeastomics,"src/function_datapub.r"))

sgd = load.sgd.proteome()
dub = load.dubreuil2019.data(3)
library(tidyverse)
sc.sti = dub %>%
  group_by(id) %>%
  summarise(sti_full = mean(stickiness),
            sti_iup = mean(stickiness[IUP20==1]),
            sti_dom =  mean(stickiness[domain==1]),
            sti_adom =  mean(stickiness[antidomain==1]))

write_rds(sc.sti,file = 'data/uniprot-stickiness.rds')
