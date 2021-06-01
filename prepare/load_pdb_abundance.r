yeastomics="https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/"
source(file.path(yeastomics,"src/utils.r"))
source(file.path(yeastomics,"src/function_annotation.r"))
source(file.path(yeastomics,"src/function_sequence.r"))
source(file.path(yeastomics,"src/function_phylogenetic.r"))
source(file.path(yeastomics,"src/function_analysis.r"))
source(file.path(yeastomics,"src/function_datalocal.r"))
source(file.path(yeastomics,"src/function_datapub.r"))
library(tidyverse)
library(hablar)

##### 1. UNIPROT-TO-ORF MAPPING
sc_uni = readRDS("data/uniprot-features.rds") %>%
  dplyr::select("UNIPROTKB","SGD","EXISTENCE","SCORE","FAMILIES","PNAME")
sc_sgd = load.sgd.features() %>% dplyr::select("sgdid","type","qual",ORF="name","gname","chr","strand")
uni2sgd = left_join(sc_uni,sc_sgd,by=c('SGD'='sgdid'))

##### 2. PAXDB ABUNDANCE RANKED
sc_pax = get.paxdb(4932,'integrated') %>%
  group_by(taxid) %>%
  mutate(rk = dense_rank(desc(ppm_int)),
         rk.pc = percent_rank(desc(ppm_int)),
         pc= 100 * ppm_int / sum_(ppm_int)) %>%
  arrange(rk.pc)

##### 3. PAXDB MAPPED TO UNIPROT-TO-ORF MAPPING
sc_pax2orf = left_join(sc_pax,uni2sgd, by=c('protid'='ORF')) %>%
  relocate(UNIPROTKB,EXISTENCE,SCORE,FAMILIES,
           ORF,gname,SGD,type,qual,chr,strand)

##### 4. PDB TO ORF MAPPING VIA UNIPROT
sc_complex = load.3dcomplex.yeast()
sc_pdb = get.mapping.3dcomplex.yeast()

##### 5. RANKED PROTEIN ABUNDANCE MAPPED TO PDB
final = left_join(sc_pax2orf,sc_pdb, by=c('UNIPROTKB'='seqid'))
