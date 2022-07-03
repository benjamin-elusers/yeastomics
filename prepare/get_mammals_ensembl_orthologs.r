library(AnnotationHub)
source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
ENS_MIRROR='asia'
count_ens_id = function(x,toprint=T){
  nr = nrow(x)
  sp_ensp = str_subset(colnames(x),"homolog_ensembl_peptide")

  ng = n_distinct(x$ensg)
  nt = n_distinct(x$enst)
  np = n_distinct(x$ensp)
  no = n_distinct(x[[sp_ensp]])
  if(!toprint)
    return(c('genes'=ng,'transcript'=nt,'protein'=np,'ortholog'=no))

  message(sprintf('rows = %7s | genes = %6s | transcripts = %6s | proteins = %6s | orthologs = %6s',nr,ng,nt,np,no))
}

get_sp = function(full_name){
  spname = str_split_fixed(full_name,"_",n=3)
  first = str_sub(spname[,1],start = 1,end=1) %>% str_to_lower()
  second = spname[,2] %>% str_to_lower()
  return(paste0(first,second))
}

match_species_names = function(ens_tree, ens_sp){

  library(stringdist)
  if(class(ens_tree) == 'phylo'){
    sp_tree =  ens_tree$tip.label %>% tolower
  }else{
    sp_tree = ens_tree
  }
  sp_ensembl = ens_sp %>% tolower

  sim_names = stringsimmatrix(sp_tree, sp_ensembl, method='jw',p=0.1,useNames='strings')
  maxS= apply(sim_names,1,max_)
  i_maxS= apply(sim_names,1,which.max)

  matched = tibble(sp_tree=sp_tree, sp_ens = sp_ensembl[i_maxS], similarity=maxS)

  return(matched)
}

# 1. get Ensembl vertebrates ---------------------------------------------------
library(treeio)
# Ensembl Vertebrates
ens_vertebrates_tree = get_ensembl_sptree('vertebrates_species-tree_Ensembl')
ens_vertebrates_nodes = make.unique(ens_vertebrates_tree$node.label)
ens_vertebrates_df = ens_vertebrates_tree %>% as_tibble() %>%
                      mutate(is_leaf = label %in% ens_vertebrates_tree$tip.label)

ens_vertebrates_info = get_ensembl_vertebrates() %>%
                       # compute similarity between species tree names and ensembl species
                       left_join(match_species_names(ens_vertebrates_tree,.$species),by=c('species'='sp_ens')) %>%
                       dplyr::rename(vertebrates_tree=sp_tree, vertebrates_similarity=similarity ) %>%
                       relocate(organism,species,vertebrates_tree,vertebrates_similarity)


# Get lineage from root node (common ancestor to vertebrate)
ens_vertebrates_clades = ape::subtrees(ens_vertebrates_tree, wait=FALSE) %>%
                          set_names(ens_vertebrates_nodes)

library(ggtree)
VERTEBRATES = ggtree(ens_vertebrates_tree,ladderize = T,right = T,branch.length = 'none') +
    ggtree::geom_hilight(node=307, alpha=0.2, type="rect") +
    ggtree::geom_nodelab(size=2.0,geom='label',node = 'internal') +
    ggtree::geom_tiplab(size=6,as_ylab = T,align = T)
ggsave(VERTEBRATES, filename='~/Desktop/VERTEBRATES_TREE.pdf', scale=2.5)

# 2. get Ensembl Eutherian Mammals (aligned genomes) ----------------------------
ens_mammals_tree = get_ensembl_sptree('43_eutherian_mammals_EPO_default')
ens_mammals_nodes = make.unique(ens_mammals_tree$node.label)

ens_mammals_info = ens_vertebrates_info %>%
                   # compute similarity between species tree names and ensembl species
                   right_join(match_species_names(ens_mammals_tree,.$species),by=c('species'='sp_ens')) %>%
                   dplyr::rename(mammals_tree=sp_tree, mammals_similarity=similarity) %>%
                   relocate(organism,species,vertebrates_tree,vertebrates_similarity,mammals_tree,mammals_similarity)

ens_mammals_clades = ape::subtrees(ens_mammals_tree, wait=FALSE) %>%
                     set_names(ens_mammals_nodes)

# Get lineage from root node (common ancestor to eutherian mammals)
ens_mammals_lineages = ape::nodepath(ens_mammals_tree, from = treeio::rootnode(ens_mammals_tree)) %>%
                       set_names(c(ens_mammals_nodes,last(ens_mammals_nodes)))

# Find the nodes 4 degrees below root (4 largest divisions)
ens_mammals_4clades = map(ens_mammals_lineages,pluck(4,.default=NA)) %>% unlist() %>% unique %>%
                      set_names( ens_mammals_nodes[.-ens_mammals_tree$Nnode-1] )

# Add column for indicating which of the 4 largest divisions this node belongs to
ens_mammals_df= groupClade(ens_mammals_tree,.node = ens_mammals_4clades) %>% as_tibble() %>%
                mutate(is_leaf = (label  %in% ens_mammals_tree$tip.label),
                       Label = label, label = tolower(Label) ) %>%
                left_join(ens_mammals_info,by=c('label'='mammals_tree')) %>%
                mutate(num_label = paste0(node,".",organism))

ens_mammals_tree$tip.label =  ens_mammals_df$num_label[ens_mammals_df$is_leaf]
MAMMALS = ggtree(ens_mammals_tree,ladderize = T,right = T,branch.length = 'none') %<+% ens_mammals_df  +
  ggtree::geom_nodelab(mapping = aes(x=branch,color=group), size=2.5, geom='label',node = 'internal', nudge_y = 0.7) +
  ggtree::geom_tiplab(aes(color=group),offset=0.1,size=3) +  xlim(0,14) +
  #ggtree::geom_tippoint(aes(size=100*f_orthologs)) +
  scale_color_metro_d() + theme(legend.position = 'none')
ggsave(MAMMALS, filename='~/Desktop/MAMMALS_TREE.pdf', scale=2.5)



# 3. Retrieve reference genome/transcriptome/proteome from ensembl -------------


# Get reference identifiers (Uniprot/ENSP = proteome, ENSG = Genome, ENST=Transcriptome)
hs_uniref = get.uniprot.proteome('9606',DNA = T,fulldesc = T) %>% names
ensp_uniref = hs_uniref %>% str_extract(ENSEMBL.nomenclature()) %>% na.omit() %>% as.vector
n_distinct(ensp_uniref)
ensg_hgnc = load.hgnc(with_protein = T,all_fields = T) %>% drop_na(ensg) %>% pull(ensg)
hs_tx = get_ensembl_tx(ENSG = ensg_hgnc, ENSP = ensp_uniref)


# 4. Query Ensembl biomart for human orthologuous proteins in mammals ----------
library(biomaRt)
ENS_MIRROR='uswest' # if asia fails
ens=useEnsembl('ensembl',mirror = ENS_MIRROR)
hs_ens=useDataset(dataset = 'hsapiens_gene_ensembl', mart = ens)
ens_hs_orthologs = get_ens_filter_ortho() # species with human orthologs data

hs_mammals = left_join(ens_mammals_info,ens_hs_orthologs, by=c('organism'='org')) %>%
             dplyr::rename(filter = name, ens_dataset = sp)

HS_QUERY = read_rds("/data/benjamin/Orthologs/results-query-ens-hs_ortho.rds")

#hs_ali = list.files("/data/benjamin/Evolution/HUMAN/fasta", pattern = '\\.mu',full.names = T)
#file.rename(hs_ali, hs_ali %>% str_replace('\\.mu','\\.fasta'))

att_species = c('ensembl_gene','associated_gene_name','ensembl_peptide',
                'canonical_transcript_protein','subtype',
                'perc_id','perc_id_r1','goc_score',
                'wga_coverage','orthology_confidence')

HS_QUERY = list()
for( f in 1:nrow(hs_mammals) ){
  filt_name = hs_mammals$filter[f]
  org = hs_mammals$organism[f]
  sp =  hs_mammals$ens_dataset[f]
  sp_att = sprintf("%s_homolog_%s",sp,att_species)

  QORTHO = tryCatch({
      Q1.0 = query_ens_ortho(species = 'hsapiens',ortho=sp,COUNTER=f) %>% as_tibble() %>%
             dplyr::rename(ensg=ensembl_gene_id,ensp=ensembl_peptide_id,enst=ensembl_transcript_id) %>%
             filter(!no_orthologs)

      col_ortho_conf = grep('orthology_confidence',sp_att,v=T)
      Q1.1 = Q1.0[ Q1.0[[col_ortho_conf]] == 1,]
      message("--> remove low confidence orthologs")

      col_ortho_prot = grep('homolog_ensembl_peptide',colnames(Q1.1),v=T)
      sp_txlen = query_ens_txlen(Sp = sp,ORG=org,COUNTER = f,verbose=F) %>%
                 filter(ensembl_peptide_id != "" & !is.na(ensembl_peptide_id))
      if( is.null(sp_txlen) ){ next }

      id_hs = c('ensg','enst','ensp')
      id_ortho = set_names(c('ensembl_gene_id','ensembl_peptide_id'),sp_att[c(1,3)])
      Q1.2 = inner_join(Q1.1,hs_tx, by=id_hs) %>%
             inner_join(sp_txlen,by=id_ortho,suffix=c('',".ortho")) %>%
              dplyr::rename(cds_len.ortho=cds_length,transcript_len.ortho=transcript_length) %>%
              as_tibble() %>% distinct()

      message("--> keep orthologs with cds closest in length")
      Q1.3 = Q1.2 %>%
             filter(!is.na(cds_len-cds_len.ortho)) %>%
             group_by(ensg,enst,ensp) %>%
             mutate(cds_diff = abs(cds_len - cds_len.ortho)) %>%
              dplyr::filter( cds_diff == min(cds_diff))
      count_ens_id(Q1.3)
      Q1.3
  },
  error=function(cond) {
    cat(sprintf('skipping species: %s!\n',org))
    message(cond)
  }
  )
  orgname = tolower(str_replace_all(org,'[ -]+','_'))
  saveRDS(QORTHO,here::here('output','ens_hs_ortho',paste0(f,'-',sp,'-',orgname,'.rds')))
  HS_QUERY[[sp]] = QORTHO
}
saveRDS(HS_QUERY,ere::here('output','ens_hs_ortho','hs_mammals_ortho.rds'))

