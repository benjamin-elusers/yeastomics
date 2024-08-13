library(here)
original_data = here("released-dataset","original_data")
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("analysis","function_evorate_fitting.R"))

DATASETS = tibble(source,name,url,accessible,local)

# SGD
SGD = tibble(
  source = c('SGD','SGD','SGD'),
  name = c('proteome','genome_coding','genomic_features'),
  url = c("http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz",
         "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz",
         "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab")
  )


# UNIPROT
# query REST API for single proteins
url_uniprot_rest_query = "https://rest.uniprot.org/uniprotkb/search?query="
UNIPROT =
  tibble(
    source = c('UNIPROT','UNIPROT','UNIPROT','UNIPROT','UNIPROT'),
    name = c('reference_proteome','reference_genome_coding','reference_mapping','full_mapping','localization'),
    url = c("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz",
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292_DNA.fasta.gz",
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.idmapping.gz",
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping.dat.gz",
            "https://rest.uniprot.org/uniprotkb/search?query=")
  )


# ENSEMBL
# query the biomart ensembl database

url_pfam_domain = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes/559292.tsv.gz"
url_pfam_clan = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
url_superfamily = "https://supfam.org/SUPERFAMILY/cgi-bin/save.cgi?var=sc;type=ass"

# Dubreuil 2019 : chemical composition of proteins
url_dubreuil2019 = "https://ndownloader.figshare.com/files/26844182"
# Leueunberger 2017 : Protein stability
url_leuenberger2017 = "https://www.science.org/doi/suppl/10.1126/science.aai7825/suppl_file/aai7825_leuenberger_table-s3.xlsx"
# Jarzab 2020 : Protein stability
url_jarzab2020 = "https://figshare.com/ndownloader/files/21653313"

# Meldal 2019 : protein complexes
url_meldal2019 = "http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv"
# STRING INTERACTIONS
url_string_functional = "https://stringdb-downloads.org/download/protein.links.full.v11.0/4932.protein.links.full.v11.0.txt.gz"
url_string_physical = "https://stringdb-downloads.org/download/protein.physical.links.full.v11.0/4932.protein.physical.links.full.v11.0.txt.gz"

# GENE ONTOLOGY
# use the package org.Sc.sgd.db

# KEGG PATHWAYS
# use the package KEGGREST (keggLink,keggConv)

byrne2005 = load.byrne2005.data()
belle2006 = load.belle2006.data()
pu2008 = load.pu2008.data()
costanzo2010 = load.costanzo2010.data()
barton2010 = load.barton2010.data()
christiano2014 = load.christiano2014.data()
geisberg2014 = load.geisberg2014.data()
dana2014 = load.dana2014.data()
lee2014 = load.lee2014.data()
filleton2015 = load.filleton2015.data()
vanleeuwen2016 = load.vanleeuwen2016.data()
villen2017 = load.villen2017.data()
leuenberger2017 = load.leuenberger2017.data()
mittal2017 = load.mittal2017.data()
meldal2019 = load.meldal.2019.data()
dubreuil2019 = load.dubreuil2019.data(4)
hausser2019 = load.hausser2019.data()
jarzab2020 = load.jarzab2020.data()
szavitsnossan2020 = load.szavitsnossan2020.data()
vanleeuwen2020 = load.vanleeuwen2020.data()
dubreuil2021 = load.dubreuil2021.data(1)

# Abundance
marguerat2012 = load.marguerat2012.data()
ho2018 = load.ho2018.data()
peter2018 = load.peter2018.data(1)


YEASTOMICS_V0 = readRDS(here::here('released-dataset','yeastOmics-290921-v0.rds'))
YEASTOMICS_V1 = readRDS(here::here('released-dataset','yeastOmics-v1.rds'))
dim(YEASTOMICS_V0)
dim(YEASTOMICS_V1)

all_features_v0 = names(YEASTOMICS_V0) %>% str_remove("cat_[^\\.]+\\.")
all_features_v1 = names(YEASTOMICS_V1)

n_distinct(all_features_v0)
n_distinct(all_features_v1)

intersect(all_features_v1,all_features_v0)

matched_features = match_strings(all_features_v1,all_features_v0,max_strings=3,use_soundex = T)
matched_features %>% filter(!is_identical)
