library(AnnotationHub)
source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
ENS_MIRROR='asia'

get_sp = function(full_name){
  spname = str_split_fixed(full_name,"_",n=3)
  first = str_sub(spname[,1],start = 1,end=1) %>% str_to_lower()
  second = spname[,2] %>% str_to_lower()
  return(paste0(first,second))
}

match_ensembl_to_sptree = function(sptree, ens_species, ens_sp, ens_names){
  sptree_names =  sptree$tip.label

  similarity_names = stringdist::stringsimmatrix(sptree_names, ens_species, method='jw',p=0.1)
  maxS= apply(similarity_names,1,max_)
  i_maxS= apply(similarity_names,1,which.max)

  matched = tibble(tree_species=sptree_names, ens_sp = ens_sp[i_maxS],
                   ens_species = ens_species[i_maxS],
                   ens_names= ens_names[i_maxS], sim=maxS)
  return(matched)
}


# 1. get Ensembl vertebrate/mammals species info ---------------------------------------------
library(treeio)
# ens_ftp= get.ensembl.species()  # get information from parsing the html page redirecting to ftp data

# Species with human orthologs data
ens_orthologs = get_ens_filter_ortho()

# Ensembl Vertebrates
ens_vertebrates_tree = get_ensembl_sptree('vertebrates_species-tree_Ensembl')
ens_vertebrates_df = ens_vertebrates_tree %>% as_tibble()
ens_vertebrates_info = get_ensembl_vertebrates()
ens_vertebrates_clades = ape::subtrees(ens_vertebrates_tree, wait=FALSE) %>% set_names(make.unique(ens_vertebrates_tree$node.label))

sum( tolower(ens_vertebrates_tree$tip.label) %pin% tolower(ens_vertebrates_info$species) )
sum( tolower(ens_vertebrates_info$species) %pin% tolower(ens_vertebrates_tree$tip.label) )

stringdist::stringsimmatrix(tolower(ens_vertebrates_tree$tip.label), (ens_vertebrates_info$species),method='jw',p=0.1,useNames='strings')




grep("cricetulus",ens_vertebrates_info$species,v=T)


library(ggtree)
VERTEBRATES = ggtree(ens_vertebrates_tree,ladderize = T,right = T,branch.length = 'none') +
    ggtree::geom_hilight(node=307, alpha=0.2, type="rect",) +
    ggtree::geom_nodelab(size=2.5,geom='label',node = 'internal',) +
    ggtree::geom_tiplab(size=2.5 )
ggsave(VERTEBRATES, filename='~/Desktop/VERTEBRATES_TREE.pdf', scale=2.5)


# Eutherian Mammals (aligned genomes)
ens_mammals_tree = get_ensembl_sptree('43_eutherian_mammals_EPO_default')
#ens_mammals_info = ens_vertebrates_info %>% mutate( similarity = mutate()
ens_mammals_nodes = make.unique(ens_mammals_tree$node.label)
ens_mammals_clades = ape::subtrees(ens_mammals_tree, wait=FALSE) %>%
                     set_names(ens_mammals_nodes)
# Get lineage from root node (common ancestor to eutherian mammals)
ens_mammals_lineages = ape::nodepath(ens_mammals_tree,from = treeio::rootnode(ens_mammals_tree)) %>%
                          set_names(ens_mammals_nodes)

# Find the nodes 4 degrees below root (4 largest divisions)
ens_mammals_4clades = map(ens_mammals_lineages,pluck(4,.default=NA)) %>% unlist() %>% unique %>%
                      set_names( ens_mammals_nodes[.-ens_mammals_tree$Nnode-1] )

# Add column for indicating which of the 4 largest divisions this node belongs to
ens_mammals_df=groupClade(ens_mammals_tree,.node = ens_mammals_4clades) %>% as_tibble()

sp_vertebrates_mammals = ens_vertebrates_clades$Eutheria$tip.label
sp_eutherian_mammals = ens_mammals_tree$tip.label


p = ggtree(ens_mammals_tree,ladderize = T,right = T,branch.length = 'none')  %<+% ens_mammals_df +
  ggtree::geom_nodelab(mapping = aes(x=branch,color=group),nudge_y = 0.7) +
  ggtree::geom_tiplab(aes(color=group),offset=0.5) +
  ggtree::geom_tippoint(aes(size=100*f_orthologs)) +
  scale_color_metro_d() +
  xlim(0,16)
p


spname_similarity = stringdist::stringsimmatrix(ens_mammals_tree$tip.label,
                                                ,
                                                method = 'jw',p=0.1) %>%
                    set_colnames(ens_vertebrates_clades$Eutheria$tip.label) %>%
                    set_rownames(ens_mammals_tree$tip.label)
maxS= apply(spname_similarity,1,max_)
i_maxS= apply(spname_similarity,1,which.max)
setNames( colnames(spname_similarity)[i_maxS], rownames(spname_similarity) )



ens_ftp$taxid = find_ncbi_id(ens_ftp$spname)


# 2. query Ensembl biomart for human orthologuous peptides in mammals ----------


# biomart with human ensembl dataset
library(biomaRt)
ens=useEnsembl('ensembl',mirror = ENS_MIRROR)
hs_ens=useDataset(dataset = 'hsapiens_gene_ensembl', mart = ens)


att_pos = c('start_position','end_position')
att_gene = c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id')
att_uni = c('uniprotswissprot')
att_struct = c('cds_length','transcript_length','transcript_start','transcript_end','ensembl_exon_id','rank','exon_chrom_start','exon_chrom_end','is_constitutive')
att_species = c('ensembl_gene','associated_gene_name','ensembl_peptide','canonical_transcript_protein','subtype',
                'perc_id','perc_id_r1','goc_score','wga_coverage','orthology_confidence')

hs_longest = get_ensembl_hs(verbose=T,longest_transcript=T) %>%
  dplyr::select(-c(start_position,end_position,transcript_start,transcript_end,
                   ensembl_exon_id,rank,is_constitutive,exon_chrom_start,exon_chrom_end,exon_length)) %>%
  distinct()

#count_ens_id(hs_longest)
# hs_longest_paxdb = dplyr::filter(hs_longest, ensembl_peptide_id %in% hs_ppm$protid )

vertebrates_ensembl = get_ensembl_sptree(treename = 'vertebrates_species-tree_Ensembl.nh')
mammals_ensembl = get_ensembl_sptree(treename = '43_eutherian_mammals_EPO_default.nh')

hs_filter_ortho = left_join(get_ensembl_vertebrates(),get_ens_filter_ortho(), by=c('name'='org')) %>%
  mutate( Species = str_to_title(species))

hs_ortho_mammals = match_ensembl_to_sptree(mammals_ensembl,hs_filter_ortho$Species,hs_filter_ortho$sp,hs_filter_ortho$name)
mammals_not_found =  hs_ortho_mammals %>%  dplyr::filter(sim<0.8) %>% pull(tree_species)

HS_QUERY = read_rds("/data/benjamin/Orthologs/results-query-ens-hs_ortho.rds")

#hs_ali = list.files("/data/benjamin/Evolution/HUMAN/fasta", pattern = '\\.mu',full.names = T)
#file.rename(hs_ali, hs_ali %>% str_replace('\\.mu','\\.fasta'))

# Ensembl
library(biomaRt)
# ens_ftp= get.ensembl.species()
# ens_ftp$taxid = find_ncbi_id(ens_ftp$spname)
# ens_ftp$sp = str_replace(ens_ftp$spname,'^([A-Z]).+ ([a-z]+)','\\1\\2') %>% str_to_lower()
ens=useEnsembl('ensembl',mirror = ENS_MIRROR)
hs_ens=useDataset(dataset = 'hsapiens_gene_ensembl', mart = ens)
#hs_ensuni = get_ensembl_hs(verbose=T)
#count_ens_id(hs_ensuni)

att_pos = c('start_position','end_position')
att_gene = c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id')
att_uni = c('uniprotswissprot')
att_struct = c('cds_length','transcript_length','transcript_start','transcript_end','ensembl_exon_id','rank','exon_chrom_start','exon_chrom_end','is_constitutive')
att_species = c('ensembl_gene','associated_gene_name','ensembl_peptide','canonical_transcript_protein','subtype',
                'perc_id','perc_id_r1','goc_score','wga_coverage','orthology_confidence')

hs_longest = get_ensembl_hs(verbose=T,longest_transcript=T) %>%
  dplyr::select(-c(start_position,end_position,transcript_start,transcript_end,
                   ensembl_exon_id,rank,is_constitutive,exon_chrom_start,exon_chrom_end,exon_length)) %>%
  distinct()

#count_ens_id(hs_longest)
# hs_longest_paxdb = dplyr::filter(hs_longest, ensembl_peptide_id %in% hs_ppm$protid )

vertebrates_ensembl = get_ensembl_sptree(treename = 'vertebrates_species-tree_Ensembl.nh')
mammals_ensembl = get_ensembl_sptree(treename = '43_eutherian_mammals_EPO_default.nh')

hs_filter_ortho = left_join(get_ensembl_vertebrates(),get_ens_filter_ortho(), by=c('name'='org')) %>%
  mutate( Species = str_to_title(species))

hs_ortho_mammals = match_ensembl_to_sptree(mammals_ensembl,hs_filter_ortho$Species,hs_filter_ortho$sp,hs_filter_ortho$name)
mammals_not_found =  hs_ortho_mammals %>%  dplyr::filter(sim<0.8) %>% pull(tree_species)

HS_QUERY = read_rds("/data/benjamin/Orthologs/results-query-ens-hs_ortho.rds")

#hs_ali = list.files("/data/benjamin/Evolution/HUMAN/fasta", pattern = '\\.mu',full.names = T)
#file.rename(hs_ali, hs_ali %>% str_replace('\\.mu','\\.fasta'))
# Ensembl
library(biomaRt)
# ens_ftp= get.ensembl.species()
# ens_ftp$taxid = find_ncbi_id(ens_ftp$spname)
# ens_ftp$sp = str_replace(ens_ftp$spname,'^([A-Z]).+ ([a-z]+)','\\1\\2') %>% str_to_lower()
ens=useEnsembl('ensembl',mirror = ENS_MIRROR)
hs_ens=useDataset(dataset = 'hsapiens_gene_ensembl', mart = ens)
#hs_ensuni = get_ensembl_hs(verbose=T)
#count_ens_id(hs_ensuni)

att_pos = c('start_position','end_position')
att_gene = c('ensembl_gene_id','ensembl_transcript_id','ensembl_peptide_id')
att_uni = c('uniprotswissprot')
att_struct = c('cds_length','transcript_length','transcript_start','transcript_end','ensembl_exon_id','rank','exon_chrom_start','exon_chrom_end','is_constitutive')
att_species = c('ensembl_gene','associated_gene_name','ensembl_peptide','canonical_transcript_protein','subtype',
                'perc_id','perc_id_r1','goc_score','wga_coverage','orthology_confidence')

hs_longest = get_ensembl_hs(verbose=T,longest_transcript=T) %>%
  dplyr::select(-c(start_position,end_position,transcript_start,transcript_end,
                   ensembl_exon_id,rank,is_constitutive,exon_chrom_start,exon_chrom_end,exon_length)) %>%
  distinct()

#count_ens_id(hs_longest)
# hs_longest_paxdb = dplyr::filter(hs_longest, ensembl_peptide_id %in% hs_ppm$protid )

vertebrates_ensembl = get_ensembl_sptree(treename = 'vertebrates_species-tree_Ensembl.nh')
mammals_ensembl = get_ensembl_sptree(treename = '43_eutherian_mammals_EPO_default.nh')

hs_filter_ortho = left_join(get_ensembl_vertebrates(),get_ens_filter_ortho(), by=c('name'='org')) %>%
  mutate( Species = str_to_title(species))

hs_ortho_mammals = match_ensembl_to_sptree(mammals_ensembl,hs_filter_ortho$Species,hs_filter_ortho$sp,hs_filter_ortho$name)
mammals_not_found =  hs_ortho_mammals %>%  dplyr::filter(sim<0.8) %>% pull(tree_species)

HS_QUERY = read_rds("/data/benjamin/Orthologs/results-query-ens-hs_ortho.rds")

#hs_ali = list.files("/data/benjamin/Evolution/HUMAN/fasta", pattern = '\\.mu',full.names = T)
#file.rename(hs_ali, hs_ali %>% str_replace('\\.mu','\\.fasta'))


HS_QUERY = list()
for( f in 1:nrow(hs_filter_ortho) ){
  filt_name = hs_filter_ortho$name[f]
  org = hs_filter_ortho$species[f]
  sp =  hs_filter_ortho$sp[f]
  sp_att = sprintf("%s_homolog_%s",sp,att_species)

  QORTHO = tryCatch({
      Q1.0 = query_ens_ortho(species = 'hsapiens',ortho=sp,COUNTER=f)
      #count_ens_id(Q1.0)

      col_ortho_conf = grep('orthology_confidence',sp_att,v=T)
      Q1.1 = Q1.0[ Q1.0[[col_ortho_conf]] == 1,]
      message("--> remove low confidence orthologs")
      #count_ens_id(Q1.1)

      col_ortho_prot = grep('homolog_ensembl_peptide',colnames(Q1.1),v=T)
      sp_long = query_ens_txlen(Sp = sp,ORG=org,COUNTER = f,verbose=F)
      #count_ens_id(sp_long)
      if( is.null(sp_long) ){ next }
      Q1.2 = inner_join(Q1.1, hs_longest, by=att_gene) %>%
             inner_join(sp_long, by=set_names(c("ensembl_gene_id","ensembl_peptide_id"),sp_att[c(1,3)]), suffix=c('',".ortho"))
      count_ens_id(Q1.2)

      message("--> keep orthologs with cds closest in length")
      Q1.3 = Q1.2 %>%
             group_by(ensembl_gene_id,ensembl_transcript_id,ensembl_peptide_id) %>%
             mutate(cds_diff = abs(cds_length - cds_length.ortho)) %>%
              dplyr::filter( cds_diff == min(cds_diff))
      count_ens_id(Q1.3)
      Q1.3
  },
  error=function(cond) {
    cat(sprintf('skipping species: %s!\n',org))
    message(cond)
  }
  )
  HS_QUERY[[sp]] = QORTHO
}
saveRDS(HS_QUERY,"/data/benjamin/Orthologs/results-query-ens-hs_ortho.rds")

