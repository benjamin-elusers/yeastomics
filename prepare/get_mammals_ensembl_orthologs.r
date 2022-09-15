library(AnnotationHub)
source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
ENS_MIRROR='uswest'
path_ortho = here::here('output','ens_hs_ortho')
library(log)
.info  = infoLog()
.error =  errorLog()
.warn  = warningLog()
.succ  = successLog()

#### FUNCTION ####
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

find_orthologs = function(x,ortho=ens_mammals_df){
  taxid = ortho$tax_id[x]
  org = ortho$organism[x]
  sp =  ortho$ens_dataset[x]
  spname =  ortho$species[x]
  treename = ortho$label[x]
  # phylum.2 = ortho$two[x]
  # phylum.4 = ortho$four[x]
  numlab = ortho$num_label[x]
  col_filter=ortho$ens_filter[x]
  att_species = c('ensembl_gene','associated_gene_name','ensembl_peptide',
                  'canonical_transcript_protein','subtype',
                  'perc_id','perc_id_r1','goc_score',
                  'wga_coverage','orthology_confidence')
  sp_att = sprintf("%s_homolog_%s",sp,att_species)

  if( is.na(sp) ){
    .warn$log(sprintf('Skip %s (%s) for orthologs...',org,spname))
    return(NULL)
  }

  QORTHO = tryCatch({
    id_hs = c('ensg','enst','ensp')
    Q1.0 = query_ens_ortho(species = 'hsapiens',sp_ortho=sp,COUNTER=x) %>% as_tibble() %>%
      dplyr::rename(ensg=ensembl_gene_id,ensp=ensembl_peptide_id,enst=ensembl_transcript_id) %>%
      filter(!no_ortholog)

    col_ortho_conf = grep('orthology_confidence',sp_att,v=T)
    Q1.1 = Q1.0[ Q1.0[[col_ortho_conf]] == 1,]
    .info$log("--> remove low confidence orthologs")

    col_ortho_prot = grep('homolog_ensembl_peptide',colnames(Q1.1),v=T)
    sp_txlen = query_ens_txlen(Sp = sp,ORG=org,COUNTER = x,verbose=F) %>%
      filter(ensembl_peptide_id != "" & !is.na(ensembl_peptide_id)) %>%
      dplyr::rename(ensg=ensembl_gene_id, enst=ensembl_transcript_id, ensp=ensembl_peptide_id) %>%
      dplyr::rename_with(.cols=!ends_with('length'), .fn = Pxx, sp, s='_') %>%
      dplyr::rename(cds_len_ortho=cds_length, tx_len_ortho = transcript_length)
    if( is.null(sp_txlen) ){ return(NULL) }

    id_ortho = set_names(paste0(sp,"_",id_hs[c(1,3)]),sp_att[c(1,3)])
    Q1.2 = inner_join(Q1.1,hs_tx, by=id_hs) %>%
      inner_join(sp_txlen,by=id_ortho) %>%
      as_tibble() %>% distinct()
    .info$log("--> keep orthologs with cds closest in length")

    col_gname = str_subset(colnames(Q1.2),'homolog_associated_gene_name$')
    col_pid = str_subset(colnames(Q1.2),'perc_id$')
    col_canonical = str_subset(colnames(Q1.2),'canonical_transcript_protein')

    Q1.3 = Q1.2 %>%
      dplyr::filter(!is.na(hs_cds_len-cds_len_ortho)) %>%
      group_by(ensg,enst,ensp) %>%
      mutate(cds_diff = abs(hs_cds_len - cds_len_ortho)) %>%
      mutate(tx_diff = abs(hs_transcript_len - tx_len_ortho)) %>%
      dplyr::filter( cds_diff == min(cds_diff)) %>%
      # in case there are >2 orthologs for a human protein, pick the one with the highest pid
      group_by(ensp) %>% dplyr::filter( .data[[col_pid]] == max_(.data[[col_pid]]) ) %>%
      # in case there are >2 orthologs for a human protein, pick the one with the gene name
      dplyr::filter( .data[[col_gname]] != "" ) %>%
      # in case there are >2 orthologs for a human protein, pick the one with closest transcript length
      dplyr::filter( tx_diff == min(tx_diff) ) %>%
      relocate(all_of(hs_cols)) %>%
      dplyr::rename(cds_len=cds_len_ortho,tx_len=tx_len_ortho) %>%
      dplyr::rename_with(.cols = !starts_with(sp) & -any_of(hs_cols), .fn = Pxx, sp, s='_' ) %>%
      ungroup() %>% distinct()

    count_ens_id(Q1.3)
    Q1.3
  },
  error=function(cond) {
    .error$log(sprintf('skipping species: %s!\n',org))
    message(cond)
  }
  )
  orgname = tolower(str_replace_all(org,'[ -]+','_'))
  QORTHO$ens_dataset = sp
  QORTHO$organism = org
  QORTHO$species = spname
  QORTHO$vertebrates_tree = treename
  QORTHO$treename = treename
  QORTHO$id_ortho = QORTHO[[paste0(sp,"_homolog_ensembl_peptide")]]
  QORTHO$pid_ortho = QORTHO[[paste0(sp,"_homolog_perc_id")]]
  QORTHO$cds_delta = QORTHO[[paste0(sp,'_cds_diff')]]
  QORTHO$tx_delta = QORTHO[[paste0(sp,'_tx_diff')]]
  QORTHO$gname_ortho = QORTHO[[paste0(sp,'_homolog_associated_gene_name')]]
  #QORTHO$two = phylum.2
  #QORTHO$four = phylum.4
  QORTHO$num_label =numlab
  QORTHO$filter = col_filter
  saveRDS(QORTHO,file.path(path_ortho,paste0(taxid,'-',sp,'-',orgname,'.rds')))
  return(QORTHO)
}

get_orthologs_clade =function(ens_ortho=df_query, clade_species, clade_name,
                              only_count=T, MIN_N=4, MIN_F=0.7){

  valid_species = intersect(tolower(clade_species),unique(ens_ortho$species))
  clade_size = n_distinct(clade_species)
  clade_nsp = n_distinct(valid_species)
  .info$log(sprintf('clade %s (n=%s/%s)...',clade_name,clade_nsp,clade_size))

  if( clade_nsp < MIN_N ){
    .warn$log(sprintf('Clade too small (less than %s species)!',MIN_N+1))
    cat('\n')
    return(NULL)
  }
  message(sprintf(' -> orthogroup size > %s or %s%% (n=%s) clade',MIN_N,MIN_F*100,round(MIN_F*clade_nsp,0)))

  hs_closest = ens_ortho %>%
    dplyr::filter( treename %in% valid_species) %>%
    mutate( clade_size = clade_nsp ) %>%
    filter(!is.na(id_ortho) & gname_ortho != '') %>%
    group_by(ensp) %>%
    arrange(cds_delta,desc(pid_ortho),tx_delta) %>%
    mutate(rk_pass = row_number(),northo=n_distinct(id_ortho),nsp=n_distinct(treename)) %>%
    ungroup() %>% arrange(ensp) %>%
    mutate(ortho_1to1 =  northo > MIN_N & northo >= MIN_F*max_(nsp) )

  if(only_count){
    clade_ortho = hs_closest %>%
      filter(ortho_1to1) %>%
      mutate(n_orthogroups = n_distinct(ensp),
             f_human = n_orthogroups / n_distinct(ens_ortho$ensp),
             clade=clade_name) %>%
      ungroup() %>%
      dplyr::select(clade,clade_size,species,treename,organism,nsp,n_orthogroups,f_human) %>%
      distinct()
    cat('\n')
    return(clade_ortho)
  }
  cat('\n')
  return(hs_closest)
}

get_mammals_sequence = function(x,ortho=ortho_prefix){

  sp = ortho$ens_dataset[x]
  message(sprintf('%2s %s',x,sp))
  id_col=ortho$col_peptide[x]
  ids = hs_ortho[[id_col]] %>% na.omit() %>% as.vector()

  mart <- useEnsembl("ensembl", dataset=paste0(sp,"_gene_ensembl"),mirror = ENS_MIRROR)
  ortho_seq = getSequence(id=ids, type = "ensembl_peptide_id", seqType = "peptide", mart = mart) %>%
    dplyr::select(ensembl_peptide_id,peptide)
  ortho_fasta = Biostrings::AAStringSet(x=ortho_seq %>% deframe())
  star = Biostrings::subseq(ortho_fasta,start=-1) == '*'
  ortho_fasta_prot = Biostrings::subseq(ortho_fasta,start=1,  end=Biostrings::width(ortho_fasta)-star)
  return(ortho_fasta_prot)
}

prepare_phylodata = function(id_hs, path_group, group_data, sp_seq=all_seq, add_human=T, overwrite=F){

  sptree =  get_ensembl_sptree('43_eutherian_mammals_EPO_default')
  sptree$node.label = NULL
  sptree$tip.label = tolower(sptree$tip.label)

  fasta_file = file.path(path_group,'fasta',paste0(id_hs,'.fa'))
  ali_file = file.path(path_group,'ali',paste0(id_hs,'.mu'))
  tree_file = file.path(path_group,'tree',paste0(id_hs,'.nh'))
  human.fa = sp_seq[ id_hs ]
  human.name = "homo_sapiens_grch38"

  group.id = group_data$id_peptide
  if(add_human){
    group.fa = setNames(c(human.fa,sp_seq[ group.id ]),c(human.name,group_data$mammals_tree))
  }else{
    group.fa = setNames(sp_seq[ group.id ],group_data$mammals_tree)
  }

  # FASTA
  if(!file.exists(fasta_file) || overwrite ){  writeXStringSet(group.fa,fasta_file,format="fasta") }

  # ALIGNMENT
  group.ali = muscle::muscle(group.fa,quiet=T)
  if(!file.exists(ali_file) || overwrite ){  writeXStringSet(as(group.ali, "AAStringSet"),ali_file,format="fasta") }

  # TREE
  group.tree=ape::keep.tip(sptree, tip = rownames(group.ali))
  if(!file.exists(tree_file) || overwrite ){ write.tree(group.tree,file = tree_file) }

}

make_mammals_fasta  = function(irow,ortho=hs_ortho,force.overwrite=F){
#  irow=15

  library(muscle)
  library(ape)
  og = ortho[irow,]
  ENSP = og$ensp

  twos = hs_closest_two %>% filter(ensp == ENSP)
  fours = hs_closest_four %>% filter(ensp == ENSP)

  if(og$n_ortholog<20){ return(NULL) }
  col_ortho_prot = str_subset(colnames(ortho),"homolog_ensembl_peptide")
  ids_mammals = og[,col_ortho_prot] %>% unlist() %>% na.omit() %>% as.vector()
  valid_ids = intersect(ids_mammals,names(all_seq))

  peptide2sp_all = strfind(valid_ids,ortho_prefix$prefix) %>%
    enframe('prefix','id_peptide') %>%
    unnest(id_peptide ) %>%
    left_join(ortho_prefix,by=c('prefix')) %>%
    dplyr::select(id_peptide,ens_dataset,prefix,mammals_tree,two,four) %>%
    filter( id_peptide %in% names(all_seq) )

  ## MAMMALS (all)
  prepare_phylodata(id_hs = ENSP, path_group = eutheria.path, group_data = peptide2sp_all, sp_seq = all_seq,
                    add_human = T, overwrite = T)

  if(og$is_best){

    peptide2sp = strfind(valid_ids,HS_EUTHERIA$prefix) %>%
      enframe('prefix','id_peptide') %>%
      unnest(id_peptide ) %>%
      left_join(HS_EUTHERIA,by=c('prefix')) %>%
      dplyr::select(id_peptide,ens_dataset,prefix,mammals_tree,two,four) %>%
      filter( id_peptide %in% names(all_seq) )

    ## EUTHERIA
    prepare_phylodata(id_hs = ENSP, path_group = eutheria.path, group_data = peptide2sp, sp_seq = all_seq,
                      add_human = T, overwrite = T)

    ## EUARCHONTOGLIRES
    euarcho.df= peptide2sp %>% filter(two=='Euarchontoglires')  %>% arrange(id_peptide)

    prepare_phylodata(id_hs = ENSP, path_group = euarcho.path, group_data = euarcho.df, sp_seq = all_seq,
                      add_human = T, overwrite = T)

    ## RODENTIA
    rodentia.df= peptide2sp %>% filter(four=='Rodentia') %>%
               left_join(fours, by=c('id_peptide'='id_ortho','ens_dataset','mammals_tree','two','four')) %>% arrange(rk_four)

    prepare_phylodata(id_hs = ENSP, path_group = rodentia.path, group_data = rodentia.df, sp_seq = all_seq,
                      add_human = F, overwrite = T)

    ## LAURASIATHERIA
    laura.df= peptide2sp %>% filter(two=='Laurasiatheria') %>%
              left_join(twos, by=c('id_peptide'='id_ortho', 'ens_dataset','mammals_tree','two','four')) %>% arrange(rk_two)

    prepare_phylodata(id_hs = ENSP, path_group = laura.path, group_data = laura.df, sp_seq = all_seq,
                      add_human = F, overwrite = T)

    ## ARTIODACTYLA
    artio.df= peptide2sp %>% filter(four=='Artiodactyla') %>%
               left_join(fours, by=c('id_peptide'='id_ortho','ens_dataset','mammals_tree','two','four')) %>% arrange(rk_four)

    prepare_phylodata(id_hs = ENSP, path_group = artio.path, group_data = artio.df, sp_seq = all_seq,
                      add_human = F, overwrite = T)

    ## PRIMATES
    primates.df= peptide2sp %>% filter(four=='Primates')  %>% arrange(id_peptide)
    prepare_phylodata(id_hs = ENSP, path_group = primates.path, group_data = primates.df, sp_seq = all_seq,
                      add_human = T, overwrite = T)

    ## CARNIVORA
    carnivora.df= peptide2sp %>% filter(four=='Carnivora')  %>%
                  left_join(fours, by=c('id_peptide'='id_ortho','ens_dataset','mammals_tree','two','four')) %>% arrange(rk_four)

    prepare_phylodata(id_hs = ENSP, path_group = carnivora.path, group_data = carnivora.df, sp_seq = all_seq,
                      add_human = F, overwrite = T)
  }
}

ens_hs_orthologs = preload(saved.file = file.path(path_ortho,'ensembl_hsapiens_filter.rds'),
                           { get_ens_filter_ortho() }, # species with human orthologs data
                           'get ensembl homolog filter...')

#### WORKFLOW ####
# 0. retrieve reference genome/transcriptome/proteome from ensembl -------------
# Get reference identifiers (Uniprot/ENSP = proteome, ENSG = Genome, ENST=Transcriptome)
hs_uniref = get.uniprot.proteome(taxid = '9606',DNA = T,fulldesc = T) %>% names
ensp_uniref = hs_uniref %>% str_extract(ENSEMBL.nomenclature()) %>% na.omit() %>% as.vector
n_distinct(ensp_uniref)
ensg_hgnc = load.hgnc(with_protein = T,all_fields = T) %>% drop_na(ensg) %>% pull(ensg)
hs_tx =preload( file.path(path_ortho,'ensembl_hsapiens_tx.rds'),
                { get_ensembl_tx(ENSG = ensg_hgnc, ENSP = ensp_uniref) %>%
                    dplyr::rename_with(.cols = ends_with('_len'), .fn = Pxx, 'hs', s='_') },
                'get ensembl human transcripts...')


# 1. get Ensembl vertebrates ---------------------------------------------------
library(treeio)
# Ensembl Vertebrates
ens_vertebrates_tree = get_ensembl_sptree('vertebrates_species-tree_Ensembl')
ens_vertebrates_tree$node.label = make.unique(ens_vertebrates_tree$node.label)
# Save all unique node labels and replace duplicated node labels by NA in tree
ens_vertebrates_nodes = ens_vertebrates_tree$node.label
is_dup_nodes = str_detect(ens_vertebrates_nodes,".+\\.[0-9]+$")
dup_nodelabels = ens_vertebrates_nodes[is_dup_nodes]
dup_nodenums  = tidytree::nodeid(ens_vertebrates_tree,dup_nodelabels)
dup_nodes = set_names(dup_nodenums,dup_nodelabels)
#ens_vertebrates_tree$node.label[is_dup_nodes] = NA


ens_vertebrates_info = get_ensembl_vertebrates() %>%
  # compute similarity between species tree names and ensembl species
  left_join(match_species_names(ens_vertebrates_tree,.$species),by=c('species'='sp_ens')) %>%
  dplyr::rename(vertebrates_tree=sp_tree, vertebrates_similarity=similarity ) %>%
  left_join(ens_hs_orthologs, by=c('organism'='org')) %>%
  dplyr::rename(ens_filter = name, ens_dataset = sp) %>%
  relocate(organism,species,ens_dataset,ens_filter,vertebrates_tree,vertebrates_similarity)

# Get lineage from root node (common ancestor to vertebrate)
ens_vertebrates_clades = ape::subtrees(ens_vertebrates_tree, wait=FALSE) %>%
  set_names( ens_vertebrates_tree$node.label ) %>%
  purrr::keep(!is_dup_nodes)

ens_vertebrates = ens_vertebrates_tree %>% as_tibble() %>%
                  mutate(depth = ape::node.depth(ens_vertebrates_tree),
                         is_leaf = depth == 1, nodes = label %>% make.unique,
                         label=tolower(label), Label = str_to_title(label) )
ens_vertebrates_df = left_join(ens_vertebrates,ens_vertebrates_info,by=c('label'='vertebrates_tree')) %>%
                     mutate(num_label = ifelse(is_leaf,paste0(node,".",organism),''))

NUM_NODES = ens_vertebrates_df$node[!ens_vertebrates_df$is_leaf]
NUM_TIPS  = ens_vertebrates_df$node[ens_vertebrates_df$is_leaf]

library(ggtree)
VERTEBRATES = ggtree(ens_vertebrates_tree,ladderize = T,right = T,branch.length = 'none')  +
  ggtree::geom_hilight(node=307, alpha=0.2, type="rect") +
  ggtree::geom_nodelab(aes(subset=!(node %in% dup_nodenums)),size=3,geom='label',node = 'internal') +
  ggtree::geom_tiplab(size=3,offset=0.5,align=T) + xlim(0,40)
ggsave(plot=VERTEBRATES, filename=file.path(path_ortho,'VERTEBRATES_TREE.pdf'), height=25, width=25)

# 1.1 query Ensembl biomart for human orthologuous proteins in vertebrates -----
library(biomaRt)
#ENS_MIRROR='uswest' # if asia fails
ens=useEnsembl('ensembl')#,mirror = ENS_MIRROR)
hs_ens=useDataset(dataset = 'hsapiens_gene_ensembl', mart = ens)

hs_cols = c('ensg','enst','ensp',"is_enspref","has_ensp","is_canonical","canonical",
            "hs_gene_len","hs_transcript_len", 'hs_cds_len',
            "has_introns","n_transcripts","n_proteins","n_exons","n_exons_mini")

#ens_ortho = ens_mammals_df %>% filter(is_leaf)
ens_ortho = ens_vertebrates_df %>% filter(is_leaf & !is.na(ens_filter))
# Save all ensembl queries to get human-mammals orthologs
HS_QUERY = preload(file.path(path_ortho,'ensembl_hs_vertebrates_orthologs.rds'),
                   { lapply(X=1:nrow(ens_ortho), FUN=function(x){ find_orthologs(x,ens_ortho) }) %>% compact },
                   'retrieve human to vertebrates orthologs...')

col_ens = c("ens_dataset","species","organism",'treename','num_label','filter',
                    'ensp','id_ortho','gname_ortho','pid_ortho','cds_delta','tx_delta')
df_query =  purrr::map(HS_QUERY, magrittr::extract, col_ens) %>% bind_rows()


# 1.2 Filter small clades  -----------------------------------------------------

# Group species by clades
CLADES_VERTEBRATES = map(names(ens_vertebrates_clades) %>% na.omit %>% as.vector,
                   ~ get_orthologs_clade(clade_species=ens_vertebrates_clades[[.x]]$tip.label, clade_name=.x)) %>%
  purrr::compact() %>%
  bind_rows()

# Count orthogroups and human orthologs per clade
CLADES_COUNT = CLADES_VERTEBRATES %>%
  left_join(ens_vertebrates,by=c('clade'='Label')) %>%
  dplyr::select(node,clade,n_orthogroups,f_human,clade_size,depth) %>%
  distinct  %>%
  arrange(f_human) %>%
  mutate(num_node = paste0(node,".",clade,'\n',n_orthogroups),
         id=node) %>% relocate(id)

p0=ggplot(CLADES_COUNT,aes(y=n_orthogroups, x=reorder(clade,f_human),color=f_human)) +
  geom_point() +
  geom_text(aes(label=paste0(node,'.',clade)), size=3.5, hjust=-0.1,angle=90) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  xlab('') + ylim(0,20000) + scale_x_discrete( expand= expansion(c(.02,.02)))

t0=ggtree(ens_vertebrates_tree,ladderize = T,right = T,branch.length = 'none') %<+% CLADES_COUNT  +
  geom_nodepoint(aes(size = f_human, color = n_orthogroups)) +
  xlim(-10,50)
t0.1 = t0 +   ggtree::geom_nodelab(aes(label=num_node),size=3,geom='text',node = 'internal',hjust=1.1,angle=25)
library(patchwork)
# Show number of orthogroups per taxonomic node
#t0.1 | p0

# Collapse nodes with less species or duplicated names
to_collapse = numeric()
T0 = t0
for(i in NUM_NODES){
  nodename = names(dup_nodes[dup_nodes == i])
  offsprings= treeio::offspring(ens_vertebrates_tree,i,tiponly=F) %>%
             setdiff(NUM_TIPS)
  cat(sprintf('node #%s (%s) : %s offspring...\n',i,nodename,length(offsprings)))

  if(mean(offsprings %in% dup_nodenums) > 0.8 | length(offsprings) < 4){
    to_collapse = c(to_collapse,i)
    T0 = T0 %>% collapse(node = i,mode = "none")
  }
}

tips_to_remove = map(.x = to_collapse,.f=~offspring(ens_vertebrates_tree, .x, tiponly=T)) %>% unlist() %>% unique


T0.0 = t0 + xlim(-10,45) +
  geom_point2(aes(subset=(node %in% to_collapse)), shape=23, size=2, fill='red') +
  ggtree::geom_nodelab(aes(label=node),size=3,geom='text',node = 'internal',hjust=1.4,angle=0) +
  ggeasy::easy_remove_legend()

T0.1 = T0 + xlim(-30,50) + ylim(1,80) +
      geom_point2(aes(subset=(node %in% to_collapse)), shape=23, size=2, fill='red')  +
      ggtree::geom_nodelab(aes(label=num_node),size=3,geom='text',node = 'internal',hjust=1.1,angle=0) +
      theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())

T0.0 | T0.1 | p0


# 2. get Ensembl Eutherian Mammals (aligned genomes) ----------------------------
#ens_mammals_tree = get_ensembl_sptree('43_eutherian_mammals_EPO_default')
ens_mammals_tree = get_ensembl_sptree("91_eutherian_mammals_EPO-Extended_default") # new in release 107
ens_mammals_tree$node.label = make.unique(ens_mammals_tree$node.label)
# Save all unique node labels and replace duplicated node labels by NA in tree
ens_mammals_nodes = ens_mammals_tree$node.label
is_dup_nodes = str_detect(ens_mammals_nodes,".+\\.[0-9]+$")
dup_nodelabels = ens_mammals_nodes[is_dup_nodes]
dup_nodenums  = tidytree::nodeid(ens_mammals_tree,dup_nodelabels)
#ens_mammals_tree$node.label[is_dup_nodes ] = NA

ens_mammals_info = ens_vertebrates_info %>%
  # compute similarity between species tree names and ensembl species
  right_join(match_species_names(ens_mammals_tree,.$species),by=c('species'='sp_ens')) %>%
  dplyr::rename(mammals_tree=sp_tree, mammals_similarity=similarity) %>%
  relocate(organism,species,vertebrates_tree,vertebrates_similarity,mammals_tree,mammals_similarity)

ens_mammals_clades = ape::subtrees(ens_mammals_tree, wait=FALSE) %>%
                     set_names( ens_mammals_tree$node.label ) %>%
                     purrr::keep(!is_dup_nodes)

ens_mammals = ens_mammals_tree %>% as_tibble() %>%
  mutate(depth = ape::node.depth(ens_mammals_tree),
         height = max(depth) - depth + 1,
         is_leaf = depth == 1, nodes = label %>% make.unique,
         label=tolower(label), Label = str_to_title(label) )

CLADES_MAMMALS = map(names(ens_mammals_clades) %>% na.omit %>% as.vector,
                         ~ get_orthologs_clade(clade_species=ens_mammals_clades[[.x]]$tip.label, clade_name=.x)) %>%
  purrr::compact() %>%
  bind_rows()

CLADES_COUNT = CLADES_MAMMALS %>%
  left_join(ens_mammals,by=c('clade'='Label')) %>%
  dplyr::select(node,clade,n_orthogroups,f_human,clade_size,depth,height) %>%
  distinct  %>%
  arrange(f_human) %>%
  mutate(num_node = paste0(clade,'\n',n_orthogroups))

p1=ggplot(CLADES_COUNT,aes(y=n_orthogroups, x=reorder(clade,f_human),color=f_human)) +
  geom_point() +
  geom_text(aes(label=clade), size=3.5, hjust=-0.1,angle=90) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  xlab('') + ylim(0,20000) + scale_x_discrete( expand= expansion(c(.02,.02)))

t1 = ggtree(ens_mammals_tree,ladderize = T,right = T,branch.length = 'none') %<+% CLADES_COUNT  +
  geom_nodepoint(aes(size = f_human, color = n_orthogroups)) +
  xlim(-5,16)

t1.1 = t1 +   ggtree::geom_nodelab(aes(label=num_node),size=3,geom='text',node = 'internal',hjust=1.1,angle=25)
t1.1 | p1

MAMMALS = ggtree(ens_mammals_tree,ladderize = T,right = T,branch.length = 'none') +
  ggtree::geom_nodelab(aes(subset = !(node %in% dup_nodenums)), size=3,geom='label',node = 'internal') +
  ggtree::geom_tiplab(size=3,offset=0.5,align=T) + xlim(0,30)
#ggsave(plot=MAMMALS, filename=file.path(path_ortho,'MAMMALS_TREE.pdf'), height=25, width=25)

# Get lineage from root node (common ancestor to eutherian mammals)
ens_mammals_lineages = ape::nodepath(ens_mammals_tree, from = treeio::rootnode(ens_mammals_tree)) %>%
  set_names(c(ens_mammals_tree$node.label,"out"))

nodenum = tidytree::nodeid(ens_mammals_tree,names(ens_mammals_clades))
#mammals_grp = dendextend::cutree(as.dendrogram(phytools::force.ultrametric(ens_mammals_tree)),k=6)
#table(mammals_grp)
#set_names( ape::node.depth(ens_mammals_tree)[nodenum], names(ens_mammals_clades))
#
# # Find the nodes 4 degrees below root (4 largest divisions)
# mammals_4 = map(ens_mammals_lineages[ names(ens_mammals_clades) ],
#                 pluck(4,.default=NA)) %>% unlist() %>% unique %>%
#             intersect(nodenum) %>%
#             set_names(treeio::nodelab(ens_mammals_tree,.))
#

mammals_4 = c('Carnivora'=121,'Artiodactyla'=103,'Primates'=156,'Rodentia'=135)
mammals_2 = c('Laurasiatheria'=95,'Euarchontoglires'=132)

# mammals_2 = map(ens_mammals_lineages[ names(ens_mammals_clades) ],
#                 pluck(3,.default=NA)) %>% unlist() %>% unique %>%
#             intersect(nodenum) %>%
#             set_names(treeio::nodelab(ens_mammals_tree,.))
#
four= groupClade(ens_mammals,.node = mammals_4) %>% as_tibble() %>% dplyr::rename(four=group)
two = groupClade(ens_mammals,.node = mammals_2) %>% as_tibble() %>% dplyr::rename(two=group)

# Add column for indicating which of the 2/4 largest divisions this node belongs to
ens_mammals_df= left_join(two,four) %>%
                mutate( depth = ape::node.depth(ens_mammals_tree),
                        is_leaf = depth == 1, #is_leaf = (label  %in% ens_mammals_tree$tip.label),
                        Label = str_to_title(label), label = tolower(Label) ) %>%
  left_join(ens_mammals_info,by=c('label'='mammals_tree')) %>%
  mutate(num_label = paste0(node,".",organism))

ens_mammals_tree$tip.label =  ens_mammals_df$num_label[ens_mammals_df$is_leaf]
MAMMALS_2 = ggtree(ens_mammals_tree,ladderize = T,right = T,branch.length = 'none') %<+% ens_mammals_df  +
  ggtree::geom_nodelab(mapping = aes(x=branch,color=two), size=2.5, geom='label',node = 'internal', nudge_y = 0.7) +
  ggtree::geom_tiplab(aes(color=two),offset=0.1,size=3) +  xlim(0,20) +
  #ggtree::geom_tippoint(aes(size=100*f_orthologs)) +
  scale_color_metro_d() + theme(legend.position = 'none')
ggsave(MAMMALS_2, filename=file.path(path_ortho,'MAMMALS_TREE_2phylum.pdf'), scale=2)
MAMMALS_4 = ggtree(ens_mammals_tree,ladderize = T,right = T,branch.length = 'none') %<+% ens_mammals_df  +
  ggtree::geom_nodelab(mapping = aes(x=branch,color=four), size=2.5, geom='label',node = 'internal', nudge_y = 0.7) +
  ggtree::geom_tiplab(aes(color=four),offset=0.1,size=3) +  xlim(0,20) +
  ggtree::geom_tippoint(size=1) +
  scale_color_metro_d() + theme(legend.position = 'none')
ggsave(MAMMALS_4, filename=file.path(path_ortho,'MAMMALS_TREE_4phylum.pdf'), scale=2)

# 2.1 query Ensembl biomart for human orthologuous proteins in mammals ----------
hs_closest_two = df_query %>%
             left_join(two,by=c('species'='label')) %>%
             filter(two != '0' ) %>%
             group_by(two,ensp) %>%
             arrange(cds_delta,desc(pid_ortho),tx_delta) %>%
             mutate(rk_two = row_number()) %>%
             ungroup() %>% arrange(ensp)
hs_closest_four = df_query %>%
                  left_join(four,by=c('species'='label')) %>%
                  filter(four != '0' ) %>%
             group_by(four,ensp) %>%
             arrange(cds_delta,desc(pid_ortho),tx_delta) %>%
             mutate(rk_four = row_number()) %>%
             ungroup() %>% arrange(ensp)

# find prefix
col_group = intersect(colnames(two),colnames(four)) %>% setdiff('label')
ortho_prefix = df_query %>%
               left_join(two,by=c('species'='label')) %>%
               left_join(four,by=c('species'='label',col_group)) %>%
               filter(two != '0' | four != '0' ) %>%
               drop_na() %>%
               mutate(n_human = n_distinct(ensp)) %>%
               group_by(ens_dataset) %>%
               mutate( prefix = hutils::longest_prefix(id_ortho) ) %>%
               #mutate( prefix = str_extract(id_ortho,'^[^0-9]+'), n_human = n_distinct(ensp) ) %>%
               group_by(ens_dataset) %>%
               mutate(n_orthologs = n_distinct(id_ortho),
                      #f_orthologs = n_orthologs/n_protein_coding, n_uniprot = n_swissprot+n_trembl,
                      f_human = n_orthologs/n_human,
                      qq_pid = toString( round(quantile(pid_ortho,c(0.01,0.05,0.25,0.75,0.95,0.99),na.rm=T)) ),
                      md_cds_diff = median_(cds_delta)) %>%
               dplyr::select(-id_ortho,-pid_ortho,-cds_delta,-tx_delta,-gname_ortho,-ensp) %>%
               group_by(two) %>% mutate( num_two= paste0(two, " (n=", n_distinct(species), ")") ) %>%
               group_by(four) %>% mutate( num_four = paste0(four, " (n=", n_distinct(species), ")") ) %>%
               distinct() %>%
               mutate(col_peptide = paste0(ens_dataset,'_homolog_ensembl_peptide')) %>%
               mutate( num_four = factor(num_four, c('out',sort(unique(num_four)))),
                       num_two = factor(num_two, c('out',sort(unique(num_two)))) )

# 3. combine orthologs data of human-to-mammals  -------------------------------
mammals = seq_along(HS_QUERY)
mammals_sp = names(HS_QUERY)
mammals_pep = paste0( names(HS_QUERY), '_homolog_ensembl_peptide')
HS_ORTHO = hs_tx
id_hs = c('ensg','enst','ensp')
for( s in mammals ){ #
  sp = mammals_sp[s]
  cat(sprintf("%3s %-20s\n",s,sp))
  col_ortho_pep = mammals_pep[s]
  count_ens_id( HS_QUERY[[s]] )
  og_col = c('ens_dataset','organsim','species','mammals_tree','two','four','num_label','filter')

  HS_ORTHO = left_join(HS_ORTHO, HS_QUERY[[s]] , by=c(colnames(hs_tx)), suffix=c('','') ) %>%
    dplyr::select(-contains(c('goc_score','wga_coverage','no_ortholog','_homolog_ensembl_gene'))) %>%
    distinct()
  print(dim(HS_ORTHO))
}
#dim(HS_ORTHO)
#colnames(HS_ORTHO)

# 4. get 1-to-1 orthologs human ------------------------------------------------
hs_ortho_1to1 = HS_ORTHO %>%
  dplyr::filter(!is.na(ensg) & !is.dup(ensp) ) %>%
  group_by(ensg,enst,ensp) %>%
  mutate(n_ortholog = sum(!is.na(c_across(ends_with('homolog_ensembl_peptide'))))) %>%
  filter(n_ortholog == max(n_ortholog) & n_ortholog != 0 )
#ggplot(ortho_prefix) +
#  geom_col( aes(y=organism,x=f_human,fill=four)) +
#  geom_vline(xintercept = c(0.5,0.6,0.7,0.75),linetype='dashed')
CUTOFF_FHUMAN = 0.7
HS_EUTHERIA = ortho_prefix %>%
  filter(f_human >CUTOFF_FHUMAN) %>%
  group_by(two) %>% mutate(num_two = sprintf("%s (n=%s)",two,n_distinct(ens_dataset))  ) %>%
  group_by(four) %>% mutate(num_four = sprintf("%s (n=%s)",four,n_distinct(ens_dataset)) ) %>%
  ungroup() %>%
  mutate( num_four = factor(num_four, c('out',sort(unique(num_four)))),
          num_two = factor(num_two, c('out',sort(unique(num_two)))) )

phylums = split(HS_EUTHERIA, HS_EUTHERIA$two) %>%
          append(split(HS_EUTHERIA, HS_EUTHERIA$four)) %>%
          compact

CUTOFF_SP_NA = 5

hs_ortho = hs_ortho_1to1 %>%
  dplyr::select( all_of(colnames(hs_tx)), 'n_ortholog', ends_with('_homolog_ensembl_peptide')) %>%
  rowwise() %>%
  mutate( f_orthogroup =  mean(!is.na(c_across(HS_EUTHERIA$col_peptide))) ) %>%
  mutate(n_artiodactyla = sum.na(c_across(all_of(phylums$Artiodactyla$col_peptide)),notNA=T),
         n_carnivora = sum.na(c_across(all_of(phylums$Carnivora$col_peptide)),notNA=T),
         n_rodentia =sum.na(c_across(all_of(phylums$Rodentia$col_peptide)),notNA=T),
         n_primates = sum.na(c_across(all_of(phylums$Primates$col_peptide)),notNA=T),
         n_laurasiatheria =sum.na(c_across(all_of(phylums$Laurasiatheria$col_peptide)),notNA=T),
         n_euarchontoglires = sum.na(c_across(all_of(phylums$Euarchontoglires$col_peptide)),notNA=T) ) %>%
  # Best orthogroups have less than 4 of the selected species (over 78% orthologs w/r to human)
  mutate( is_best = sum(is.na(c_across(HS_EUTHERIA$col_peptide))) < CUTOFF_SP_NA ) %>%
  ungroup() %>% mutate( p_top = percent_rank(f_orthogroup) )


library(ggplot2)
count_hs_ortho = hs_ortho %>% group_by(n_ortholog) %>% summarize(total=n())
count_hs_min = hs_ortho %>% filter(n_ortholog>20) %>% group_by(n_ortholog) %>% summarize(total=n())
count_hs_best = hs_ortho %>% filter(is_best) %>% group_by(n_ortholog) %>% summarize(total=n())
MIN_SP = min(count_hs_min$n_ortholog)
BEST_SP = min(count_hs_best$n_ortholog)

# hs_ortho_1to1 %>% filter(is_best) %>% ungroup()  %>% summarize(sum(is_best))
# sum(count_hs_best$total)

n_best_OG = sum(count_hs_best$total)
n_min_OG = sum(count_hs_min$total)
n_OG = sum(count_hs_ortho$total[count_hs_ortho$n_ortholog>=MIN_SP])
n_max_OG = sum(count_hs_ortho$total)

library(scales)
pal = palette_pander(10)
#show_col(pal)

MIN_COL = pal[2]
BEST_COL= pal[6]
p1 = ggplot(count_hs_ortho ) +
  geom_vline(xintercept=MIN_SP, linetype='dotted', color=MIN_COL) +
  geom_vline(xintercept=BEST_SP, linetype='dotted', color=BEST_COL) +
  geom_col(aes(y=total,x=n_ortholog)) +
  geom_col(data=count_hs_best,aes(y=total,x=n_ortholog),fill=BEST_COL) +
  xlab('# mammals species with orthologs') +
  ylab('# orthogroups') +
  geom_text(data = data.frame(),aes(x=MIN_SP,y=Inf,label=sprintf('>%s orthologs\n(n=%s)',MIN_SP,n_OG)),
            size=4,color=MIN_COL,vjust=2.5,hjust=1.1) +
  geom_text(data = data.frame(),aes(x=BEST_SP,y=Inf,label=sprintf('Best\nspecies\n(n=%s)',n_best_OG)),
            size=4,color=BEST_COL,vjust=1.5,hjust=-0.1) +
  geom_text(data = data.frame(),aes(x=-Inf,y=Inf, label=sprintf('\nOrthogroups=%s\n',n_max_OG)),
            hjust=-.1,vjust='inward',size=4)

# 5. Download human fasta sequences --------------------------------------------
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",mirror = ENS_MIRROR)
hs_seq = getSequence(id=hs_ortho$ensp, type="ensembl_peptide_id", seqType="peptide", mart=mart) %>%
         dplyr::select(ensembl_peptide_id,peptide) %>% deframe()
hs_fasta = Biostrings::AAStringSet(x=hs_seq)
star = Biostrings::subseq(hs_fasta,start=-1) == '*'
hs_fasta_prot = Biostrings::subseq(hs_fasta,start=1,  end=Biostrings::width(hs_fasta)-star)

# 6. Download ortholog fasta sequences -----------------------------------------

mammals_seq =preload(file.path(path_ortho,'ensembl_mammals_seq.rds'),
                     { lapply(1:nrow(ortho_prefix), get_mammals_sequence) },
                     'get mammals ortholog sequences...')

p2 = ggplot(ortho_prefix, aes(y=reorder(organism,f_human),x=f_human, fill=num_two)) +
  geom_point(shape=21,size=2) + ylab('Mammals Orthologuous species') +
  theme( axis.text.y = element_text(size=9), axis.title = element_text(size=10))+
  geom_vline(xintercept=0.78, linetype='dashed') +
  scale_color_metro_d() + scale_fill_metro_d()
p3 = ggplot(HS_EUTHERIA, aes(y=reorder(organism,f_human),x=f_human, fill=num_two)) +
  geom_point(shape=21,size=3,color=BEST_COL) + ylab('Mammals Orthologuous species') +
  theme( axis.text.y = element_text(size=9), axis.title = element_text(size=10)) +
  scale_color_metro_d() + scale_fill_metro_d()


p2.1 = ggplot(ortho_prefix, aes(y=reorder(organism,f_human),x=f_human, fill=num_four)) +
  geom_point(shape=21,size=2) + ylab('Mammals Orthologuous species') +
  theme( axis.text.y = element_text(size=9), axis.title = element_text(size=10))+
  geom_vline(xintercept=0.78, linetype='dashed') +
  scale_color_metro_d() + scale_fill_metro_d()
p3.1 = ggplot(HS_EUTHERIA, aes(y=reorder(organism,f_human),x=f_human, fill=num_four)) +
  geom_point(shape=21,size=3,color=BEST_COL) + ylab('Mammals Orthologuous species') +
  theme( axis.text.y = element_text(size=9), axis.title = element_text(size=10)) +
  scale_color_metro_d() + scale_fill_metro_d()

library(patchwork)
graphics.off()
MAMMALS_2 = MAMMALS_2 +
  geom_tippoint(data = . %>% filter(label %in% HS_EUTHERIA$num_label), color=BEST_COL,size=3)

MAMMALS_4 = MAMMALS_4 +
  geom_tippoint(data = . %>% filter(label %in% HS_EUTHERIA$num_label), color=BEST_COL,size=3)

PP = ((p1) | (MAMMALS_2)) / (p2 | p3)
ggplot2::ggsave(PP, height=12,width=18, filename = file.path(path_ortho,'Mammals-orthologs-2_phylum.pdf'))
PP.1 = ((p1) | (MAMMALS_4)) / (p2.1 | p3.1)
ggplot2::ggsave(PP.1, height=12,width=18, filename = file.path(path_ortho,'Mammals-orthologs-4_phylum.pdf'))

# 7. Write the fasta for mammals and subphylum ---------------------------------
hs_ortho_1to1.rds=file.path(path_ortho,'hs_ortho_1to1.rds')
hs_ortho.rds=file.path(path_ortho,'hs_orthologs.rds')
ortho_prefix.rds = file.path(path_ortho,'ensembl_species.rds')
HS_EUTHERIA.rds=file.path(path_ortho,'ensembl_mammals.rds')
hs_fasta.rds =file.path(path_ortho,'ensembl_human_seq.rds')
hs_closest.rds = file.path(path_ortho,'hs_closest_ortho.rds')

#saveRDS(hs_ortho_1to1, hs_ortho_1to1.rds)
#saveRDS(hs_ortho, hs_ortho.rds)
#saveRDS(ortho_prefix, ortho_prefix.rds)
#saveRDS(HS_EUTHERIA, HS_EUTHERIA.rds)
#saveRDS(hs_fasta_prot, hs_fasta.rds)
#saveRDS(list(two=hs_closest_two,four=hs_closest_four),hs_closest.rds)

#hs_ortho_1to1=readRDS(hs_ortho_1to1.rds)
#HS_EUTHERIA=readRDS(HS_EUTHERIA.rds)
#ortho_prefix=readRDS(ortho_prefix.rds)
#hs_closest_two = readRDS(hs_closest.rds)$two
#hs_closest_four = readRDS(hs_closest.rds)$four

# get all sequences
mammals_seq=readRDS(file.path(path_ortho,'ensembl_mammals_seq.rds'))
all_seq=readRDS(hs_fasta.rds)
for( i in 1:length(mammals_seq) ){ all_seq = AAStringSet(c(all_seq,mammals_seq[[i]])) }

phylum.path = sapply( phylums, function(x){ file.path(path_ortho,paste0("hs_",x)) })
eutheria.path  = file.path(path_ortho,'hs_eutheria')
boreo.path   = file.path(path_ortho,'hs_boreoeutheria')
rodentia.path    = file.path(path_ortho,'hs_rodentia')
artio.path     = file.path(path_ortho,'hs_artiodactyla')
primates.path  = file.path(path_ortho,'hs_primates')
carnivora.path = file.path(path_ortho,'hs_carnivora')
euarcho.path  = file.path(path_ortho,'hs_euarchontoglires')
laura.path = file.path(path_ortho,'hs_laurasiatheria')


all_paths = c(mammals0.path,mammals.path,glires.path,laura1.path,primates.path,carnivora.path,euarcho.path,laura.path)
all_folders = c(file.path(all_paths,'fasta'),file.path(all_paths,'ali'),file.path(all_paths,'tree'))
sapply(all_folders, dir.create, showWarnings = F,recursive = T)

library(muscle)
library(pbmcapply)
STEPS=2000
tictoc::tic('align orthologs...')
for( S in 1:10 ){
  ii = seq(STEPS*(S-1)+1,len=STEPS)
  iseq = ii[ ii<=nrow(hs_ortho_1to1)]
  print(range(iseq))
  make_fasta <- pbmclapply(X=iseq, FUN=function(x){ make_mammals_fasta(x) }, mc.cores=14, mc.cleanup=T)
}
tictoc::toc()

check_fasta_to_ali = function(path){

  fasta = list.files(file.path(path,'ali'),pattern='.mu$') %>% str_replace("\\.mu",'')
  ali = list.files(file.path(path,'fasta'),pattern='.fa$') %>% str_replace("\\.fa",'')
  trees = list.files(file.path(path,'tree'),pattern='.nh$') %>% str_replace("\\.nh",'')

  n_ali = n_distinct(ali)
  n_fa = n_distinct(fasta)
  n_tr = n_distinct(trees)
  message(sprintf("alignments : %10s",n_fa))
  message(sprintf("fasta      : %10s",n_ali))
  message(sprintf("trees      : %10s",n_tr))

  f2a = setdiff(fasta,ali)
  a2f = setdiff(ali,fasta)
  if( n_ali < n_fa){ message( paste0(f2a,"\n") ) ; return(f2a) }
  if( n_ali > n_fa){ message( paste0(a2f,"\n") ); return(a2f) }

}

mapply(check_fasta_to_ali,all_paths)
best = hs_ortho %>% filter( is_best)
table(best$n_glires)
table(best$n_carnivora)
table(best$n_laura.1)
table(best$n_primates)
table(best$n_laurasiatheria)
table(best$n_euarchontoglires)
