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
mammals_dir = "/media/WEXAC_data/MAMMALS/"
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

