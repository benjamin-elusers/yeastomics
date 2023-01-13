library(AnnotationHub)
source(here::here("src","__setup_yeastomics__.r"))
library(tidyverse)
library(here)
ENS_MIRROR='uswest'

#### FUNCTION ####
count_ens_id = function(x,toprint=T){
  nr = nrow(x)

  id_sp = str_subset(colnames(x),"homolog_ensembl_peptide")
  if( length(id_sp) == 0L ){ id_sp="id_ortho" }

  ng = n_distinct(x$ensg)
  nt = n_distinct(x$enst)
  np = n_distinct(x$ensp)
  no = n_distinct(x[[id_sp]])
  if(!toprint)
    return(c('genes'=ng,'transcript'=nt,'protein'=np,'ortholog'=no))

  message(sprintf('rows = %7s | genes = %6s | transcripts = %6s | proteins = %6s | orthologs = %6s',nr,ng,nt,np,no))
}

find_orthologs = function(ortho,BM=sc_biomart,use_cache=T){

  current_org = ortho$ens_org
  current_sp = ortho$ens_sp
  current_dataset = paste0(current_sp,"_eg_gene") # "_gene_ensembl" => for vertebrates
  orgname = tolower(str_replace_all(current_org,'[ -]+','_'))
  cache_file = file.path(path_ortho,paste0(ortho$tax_id,'-',current_sp,'-',orgname,'.rds'))

  if(use_cache & file.exists(cache_file)){
    cache_date = file.info(cache_file)[,'ctime'] %>% format(format="%d-%m-%Y %H:%M:%S")
    .info$log(sprintf("Using cached data for '%s' created on %s",current_org,cache_date))
    return(readRDS(cache_file))
  }

  sc_cols = c('ensg','enst','ensp',"is_enspref","has_ensp","is_canonical","canonical",
              "sc_gene_len","sc_transcript_len", 'sc_cds_len',
              "has_introns","n_transcripts","n_proteins","n_exons","n_exons_mini")

  att_species = c('ensembl_gene','ensembl_peptide',
                  'canonical_transcript_protein','associated_gene_name','subtype',
                  'perc_id','perc_id_r1','goc_score',
                  'wga_coverage','orthology_confidence')
  sp_att = sprintf("%s_eg_homolog_%s",current_sp,att_species) # _homolog_ for vertebrates


  ens_id = tibble(seqtype = c('gene','peptide','transcript'),
                  all = c('ensg','ensp','enst'),
                  sc = sprintf("ensembl_%s_id",seqtype),
                  ortho = c(sp_att[1:2], paste0(current_sp,'_homolog_ensembl_transcript')))

  bm_ortho = tryCatch({
    get_ensembl_dataset("fungi_mart", current_dataset) # ENSEMBL_MART_ENSEMBL for vertebrates
  },error=function(cond) {
    .error$log(sprintf('skipping species: %s!\n',current_org))
    message(cond)
    return(NULL)
  })

  if( is.na(current_sp) | is.null(bm_ortho) ){
    .warn$log(sprintf('Skip %s (%s) for orthologs...',current_org,current_sp))
    return(NULL)
  }

  QORTHO = tryCatch({
    Q1.0 = query_ens_ortho(sp_ortho=current_sp, species='scerevisiae', BIOMART = BM,COUNTER=1)

    if( is.null(Q1.0) ){ return(NULL) }

    col_ortho_conf = grep('orthology_confidence',sp_att,v=T)
    Q1.1 = Q1.0[ Q1.0[[col_ortho_conf]] == 1,] %>%
           dplyr::rename(ens_id %>% pull(sc,all)) %>%
           filter(!no_ortholog)
    .info$log("--> remove low confidence orthologs")

    col_ortho_prot = grep('homolog_ensembl_peptide',colnames(Q1.1),v=T)

    id_ortho = ens_id %>% pull(all,ortho)
    sp_txlen = query_ens_txlen(ORG=current_org, BIOMART=bm_ortho, COUNTER = 1, verbose=F, debug=T) %>%
      dplyr::rename(ens_id %>% pull(sc,all)) %>%
      filter(ensp != "" & !is.na(ensp)) %>%
      dplyr::rename(cds_len_ortho=cds_length, tx_len_ortho = transcript_length) %>%
      dplyr::rename( id_ortho )

    if( is.null(sp_txlen) ){ return(NULL) }

    #id_ortho = set_names(paste0(current_sp,"_",id_hs[c(1,3)]),sp_att[c(1,3)])

    Q1.2 = inner_join(Q1.1,sc_tx, by=ens_id$all ) %>%
      inner_join(sp_txlen,by=names(id_ortho)[1:2]) %>%
      as_tibble() %>% distinct() %>% type_convert()
    .info$log("--> keep orthologs with cds closest in length")

    col_gname = str_subset(colnames(Q1.2),'homolog_associated_gene_name$')
    col_id = str_subset(colnames(Q1.2),'ensembl_peptide$')
    col_pid = str_subset(colnames(Q1.2),'perc_id$')
    col_canonical = str_subset(colnames(Q1.2),'canonical_transcript_protein')

    ortho_cols= c('id_ortho','pid_ortho','cds_delta','tx_delta','gname_ortho')
    Q1.3 = Q1.2 %>%
      dplyr::filter(!is.na(sc_cds_len-cds_len_ortho)) %>%
      group_by(ensg,enst,ensp) %>%
      mutate(cds_delta = abs(sc_cds_len - cds_len_ortho)) %>%
      mutate(tx_delta = abs(sc_transcript_len - tx_len_ortho)) %>%
      dplyr::filter( cds_delta == min(cds_delta)) %>%
      # in case there are >2 orthologs for a yeast protein, pick the one with the highest pid
      group_by(ensp) %>% dplyr::filter( .data[[col_pid]] == max_(.data[[col_pid]]) ) %>%
      # in case there are >2 orthologs for a yeast protein, pick the one with the gene name
      #dplyr::filter( .data[[col_gname]] != "" ) %>%
      # in case there are >2 orthologs for a yeast protein, pick the one with closest transcript length
      dplyr::filter( tx_delta == min(tx_delta) ) %>%
      relocate(all_of(sc_cols)) %>%
      dplyr::rename(cds_len=cds_len_ortho,tx_len=tx_len_ortho,
                    id_ortho=!!col_id, pid_ortho=!!col_pid, gname_ortho = !!col_gname) %>%
      dplyr::rename_with(.cols = !starts_with(current_sp) & -any_of(c(sc_cols,ortho_cols)), .fn = Pxx, current_sp, s='_' ) %>%
      ungroup() %>% distinct()

    count_ens_id(Q1.3)
    Q1.3
  },
  error=function(cond) {
    .error$log(sprintf('skipping species: %s!\n',current_org))
    message(cond)
    return(NULL)
  }
  )

  if(!is.null(QORTHO)){
    QRES = bind_cols(ortho,QORTHO) %>% mutate(ens_dataset=current_sp)
    saveRDS(QRES,cache_file)
  }
  return(QRES)
}

get_orthologs_clade =function(ens_ortho=df_query, clade_species, clade_name,
                              only_count=T, MIN_N=4, MIN_F=0.7){

  valid_species = intersect(tolower(clade_species),unique(ens_ortho$label))
  clade_size = n_distinct(clade_species)
  clade_nsp = n_distinct(valid_species)
  .info$log(sprintf('clade %s (n=%s/%s)...',clade_name,clade_nsp,clade_size))

  if( clade_nsp < MIN_N ){
    .warn$log(sprintf('Clade too small (less than %s species)!',MIN_N+1))
    cat('\n')
    return(NULL)
  }
  message(sprintf(' -> orthogroup size > %s or %s%% (n=%s) clade',MIN_N,MIN_F*100,round(MIN_F*clade_nsp,0)))

  sc_closest = ens_ortho %>%
    dplyr::filter( label %in% valid_species) %>%
    mutate( clade_size = clade_nsp ) %>%
    filter(!is.na(id_ortho) & gname_ortho != '') %>%
    group_by(ensp) %>%
    arrange(cds_delta,desc(pid_ortho),tx_delta) %>%
    mutate(rk_pass = row_number(),northo=n_distinct(id_ortho),nsp=n_distinct(label)) %>%
    ungroup() %>% arrange(ensp) %>%
    mutate(ortho_1to1 =  northo > MIN_N & northo >= MIN_F*max_(nsp) )

  if(only_count){
    clade_ortho = sc_closest %>%
      filter(ortho_1to1) %>%
      mutate(n_orthogroups = n_distinct(ensp),
             f_yeast = n_orthogroups / n_distinct(ens_ortho$ensp),
             clade=clade_name) %>%
      ungroup() %>%
      dplyr::select(clade,clade_size,species,label,ens_org,nsp,n_orthogroups,f_yeast) %>%
      distinct()
    cat('\n')
    return(clade_ortho)
  }
  cat('\n')
  return(sc_closest)
}


get_fungi_sequence = function(x,ortho=ortho_prefix){

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
  prepare_phylodata(id_hs = ENSP, path_group = mammals0.path, group_data = peptide2sp_all, sp_seq = all_seq,
                    add_human = T, overwrite = T)

  if(og$is_best){

    # peptide2sp = strfind(valid_ids,HS_BOREOEUTHERIA$prefix) %>%
    #   enframe('prefix','id_peptide') %>%
    #   unnest(id_peptide ) %>%
    #   left_join(HS_BOREOEUTHERIA,by=c('prefix')) %>%
    #   dplyr::select(id_peptide,ens_dataset,prefix,mammals_tree,two,four) %>%
    #   filter( id_peptide %in% names(all_seq) )
  }
}


#### WORKFLOW ####
path_ortho = here::here('output','ens_sc_ortho','release-55')
sc_biomart = get_ensembl_dataset('fungi_mart','cerevisiae')
#sp_biomart = get_ensembl_dataset('fungi_mart','pombe')

ens_sc_orthologs = preload(saved.file = file.path(path_ortho,'ensembl_scerevisiae_filter.rds'),
                           { get_ens_filter_ortho(BIOMART = sc_biomart) }, # species with yeast orthologs data
                           'get ensembl homolog filter...')

# 0. retrieve reference genome/transcriptome/proteome from ensembl -------------
# Get reference identifiers (Uniprot/ENSP = proteome, ENSG = Genome, ENST=Transcriptome)
#find.uniprot_refprot(c('cerevisiae','yeast','s288c'))
orfref = load.sgd.proteome() %>% names
sc_uniref = get.uniprot.proteome(taxid = '559292',DNA = T,fulldesc = F) %>% names
ensp_uniref = get.uniprot.mapping(taxid = 559292, targetdb='EnsemblGenome_PRO') %>% pull(extid)
ensg_uniref = get.uniprot.mapping(taxid = 559292, targetdb='EnsemblGenome') %>% pull(extid)

sgd_gene = load.sgd.features() %>%
           mutate( gname=ifelse(gname=="",NA,gname),
                   gene = coalesce(gname,name)) %>%
           dplyr::pull(gene)

sc_tx =preload( file.path(path_ortho,'ensembl_scerevisiae_tx.rds'),
                { get_ensembl_tx(ENSG=orfref, ENSP=orfref, BIOMART=sc_biomart) %>%
                    dplyr::rename_with(.cols = ends_with('_len'), .fn = Pxx, 'sc', s='_') },
                'get ensembl yeast transcripts...')


# 1. get Ensembl Fungi ---------------------------------------------------
library(treeio)
# Ensembl Fungi
fu_tree = get_ensembl_sptree(treename = 'fungi_protein-trees_default.nh',
                             URL_SPTREE = 'http://ftp.ebi.ac.uk/ensemblgenomes/pub/fungi/current/compara/species_trees/')

fu_node = tibble( node = fu_tree$node.label,
                  node_ = make.unique(node), is_duplicated = duplicated(node))

ens_fu = get_ensembl_fungi() %>% mutate( sp = get_ensembl_sp(species) ) %>%
         filter(peptide_compara == 'Y' & sp %in% ens_sc_orthologs$sp)

# based on ncbi taxid and synonyms
fu_names = tribble(~common_name,       ~sci_name,
                   "ashbya_gossypii",                "eremothecium_gossypii_atcc_10895",
                   "aspergillus_fumigatus",          "aspergillus_fumigatus_af293",
                   "aspergillus_fumigatusa1163",     "aspergillus fumigatus a1163",
                   "candida_auris",                  "_auris_strain_strain_b8441",
                   "candida_duobushaemulonis",       "_duobushaemulonis_strain_strain_b09383",
                   "candida_glabrata",               "_glabrata_cbs_138",
                   "candida_haemuloni",              "_haemuloni_strain_b11899",
                   "candida_pseudohaemulonis",       "_pseudohaemulonii_strain_b12108",
                   "colletotrichum_gloeosporioides", "colletotrichum_fructicola_nara_gc5",
                   "fusarium_solani",                "fusarium_vanettenii_77-13-4_strain_mpvi_77-13-4",
                   "gaeumannomyces_graminis",        "gaeumannomyces_tritici_r3-111a-1",
                   "komagataella_pastoris",          "komagataella_phaffii_gs115",
                   "magnaporthe_oryzae",             "pyricularia_oryzae_70-15",
                   "magnaporthe_poae",               "magnaporthiopsis_poae_atcc_64411",
                   "melampsora_laricipopulina",      "melampsora_larici-populina_98ag31",
                   "neosartorya_fischeri",           "aspergillus_fischeri_nrrl_181",
                   "phaeosphaeria_nodorum",          "parastagonospora_nodorum_sn15",
                   "puccinia_graminis",              "puccinia_graminis_f._sp._tritici_crl_75-36-700-3",
                   "puccinia_graminisug99",          "puccinia_graminis_f._sp._tritici_04ken156/4",
                   "pyrenophora_triticirepentis",    "pyrenophora_tritici-repentis_pt-1c-bfp",
                   "pyricularia_oryzae",             "pyricularia_oryzae_strain_br32",
                   "verticillium_dahliae",           "verticillium_dahliae_vdls.17",
                   "verticillium_dahliaejr2",        "verticillium_dahliae_jr2",
)


fu_species = match_strings(ens_fu$species,fu_tree$tip.label,use_soundex = F) %>%
             bind_rows(fu_names %>% dplyr::rename(s1='common_name')) %>%
             dplyr::filter( is_matched | !is.na(sci_name)) %>%
             group_by(s1) %>% add_count(name='n1') %>%
             mutate( s2 = ifelse(is.na(s2), sci_name, s2),
                     is_matched = n1 == 1,
                     verified = ifelse( is.na(sci_name), 'by_similarity', 'taxid_synonym') ) %>%
             dplyr::select(-sci_name)



fu_info = inner_join(ens_sc_orthologs,ens_fu, by=c('Org'='organism','sp')) %>%
          left_join(fu_species , by = c('species'='s1')) %>%
          dplyr::rename(fu_tree=s2,ens_org=org,ens_filter=name,ens_sp=sp) %>%
          filter(!is.na(fu_tree)) %>%
          relocate(tax_id, Org, ens_org, ens_sp,species,fu_tree,ens_filter,
                   osa:soundex, is_identical:is_matched, n1, verified,
                   species_id:coverage,division,cored_db,
                   assembly_version:genebuild,variation:other_alignments
          )


# Get lineage from root node (common ancestor to vertebrate)
fu_clades = ape::subtrees(fu_tree, wait=FALSE) %>%
            set_names( fu_tree$node.label ) %>%
            purrr::discard(fu_node$is_duplicated)

fu_data = fu_tree %>% as_tibble() %>%
          mutate(depth = ape::node.depth(fu_tree),
          is_leaf = depth == 1, nodes = label %>% make.unique,
          label=tolower(label), Label = str_to_title(label) )

fu_df = left_join(fu_data,fu_info,by=c('label'='fu_tree')) %>%
  mutate(num_label = ifelse(is_leaf,paste0(node,".",ens_org),''))

NUM_NODES = fu_df$node[!fu_df$is_leaf]
NUM_TIPS  = fu_df$node[fu_df$is_leaf]

library(ggtree)
dup_nodes = fu_node$node[fu_node$is_duplicated]
sc_ortho_node = fu_data$node[ fu_data$is_leaf & fu_data$label %in% fu_species$s2]

FUNGI = ggtree(fu_tree,ladderize = T,right = T,branch.length = 'none')  +
  ggtree::geom_tippoint(aes(subset=(node %in% sc_ortho_node)),size=1, color='red') +
  ggtree::geom_nodelab(aes(subset=!(node %in% dup_nodes)),size=3,geom='label',node = 'internal') +
  ggtree::geom_tiplab(size=3,offset=0.5,align=T) + xlim(0,30) + ylim(0,400)
ggsave(plot=FUNGI, filename=file.path(path_ortho,'FUNGI_TREE.pdf'), height=25, width=25)

# 1.1 query Ensembl biomart for yeast orthologuous proteins in fungi -----
ens_ortho = fu_df %>% filter(is_leaf & !is.na(ens_filter)) %>%
            mutate(ens_dataset = paste0(ens_sp,'_eg_gene'))

# Save all ensembl queries to get human-mammals orthologs
SC_QUERY = preload(file.path(path_ortho,'ensembl_sc_fungi_orthologs.rds'),
                   { lapply(X=1:nrow(ens_ortho), FUN=function(x){ find_orthologs(ortho = ens_ortho[x,],BM = sc_biomart) }) %>% compact },
                   'retrieve yeast to fungi orthologs...')

# Get lineage from root node (common ancestor to fungi)
#fu_lineages = ape::nodepath(fu_tree, from = treeio::rootnode(fu_tree)) %>%
#              set_names(fu_tree$tip.label)

fu_lineages = find.lineage(fu_tree,full = F)

col_ens = c(colnames(ens_ortho),'ensp','id_ortho','gname_ortho','pid_ortho','cds_delta','tx_delta')
df_query =  purrr::map(SC_QUERY, magrittr::extract, col_ens) %>% bind_rows()

CLADES_FUNGI = map(names(fu_clades) %>% na.omit %>% as.vector,
                   ~ get_orthologs_clade(clade_species=fu_clades[[.x]]$tip.label,
                                         clade_name=.x,MIN_N = 1,MIN_F = 0.1)) %>%
  purrr::compact() %>%
  bind_rows()

CLADES_COUNT = CLADES_FUNGI %>%
  left_join(fu_df,by=c('clade'='Label')) %>%
  dplyr::select(clade,n_orthogroups,f_yeast,clade_size,depth) %>%
  distinct  %>%
  arrange(f_yeast) %>%
  mutate(num_node = paste0(clade,'\n',n_orthogroups))

p0=ggplot(CLADES_COUNT,aes(y=n_orthogroups, x=reorder(clade,f_yeast),color=f_yeast)) +
  geom_point() +
  geom_text(aes(label=clade), size=3.5, hjust=-0.1,angle=90) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  xlab('') + ylim(0,6000) + scale_x_discrete( expand= expansion(c(.02,.02)))

t0=ggtree(fu_tree,ladderize = T,right = T,branch.length = 'none') %<+% CLADES_COUNT  +
  geom_nodepoint(aes(size = f_yeast, color = n_orthogroups)) +
  xlim(-10,50)
t0.1 = t0 +   ggtree::geom_nodelab(aes(label=num_node),size=3,geom='text',node = 'internal',hjust=1.1,angle=25)
library(patchwork)
# Show number of orthogroups per taxonomic node
#t0.1 | p0

T0=t0
for(i in dup_nodes){
  T0 = T0 %>% collapse(node = i,mode = "none")
}
T0.1 = T0 + ylim(1,55) + xlim(-20,50) +
  geom_point2(aes(subset=(node %in% dup_nodes)), shape=23, size=2, fill='red')  +
  ggtree::geom_nodelab(aes(label=num_node),size=3,geom='text',node = 'internal',hjust=1.1,angle=0) +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())

#T0.1 | p0

#test = read.tree( file.path(path_ortho,'ens_vertebrates_species_time_in_MYA.nwk'))
#test %>% as_tibble
#ggtree(test) + geom_tiplab() + geom_nodelab(geom='shadowtext') + xlim(0,1500)



# Find the nodes 4 degrees below root (4 largest divisions)
mammals_4 = map(ens_mammals_lineages[ names(ens_mammals_clades) ],
                pluck(4,.default=NA)) %>% unlist() %>% unique %>%
  intersect(nodenum) %>%
  set_names(treeio::nodelab(ens_mammals_tree,.))

mammals_2 = map(ens_mammals_lineages[ names(ens_mammals_clades) ],
                pluck(3,.default=NA)) %>% unlist() %>% unique %>%
  intersect(nodenum) %>%
  set_names(treeio::nodelab(ens_mammals_tree,.))

four= groupClade(ens_mammals_tree,.node = mammals_4) %>% as_tibble() %>% dplyr::rename(four=group)
two = groupClade(ens_mammals_tree,.node = mammals_2) %>% as_tibble() %>% dplyr::rename(two=group)

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
  ggtree::geom_tiplab(aes(color=two),offset=0.1,size=3) +  xlim(0,15) +
  #ggtree::geom_tippoint(aes(size=100*f_orthologs)) +
  scale_color_metro_d() + theme(legend.position = 'none')
ggsave(MAMMALS_2, filename=file.path(path_ortho,'MAMMALS_TREE_2phylum.pdf'), scale=2)
MAMMALS_4 = ggtree(ens_mammals_tree,ladderize = T,right = T,branch.length = 'none') %<+% ens_mammals_df  +
  ggtree::geom_nodelab(mapping = aes(x=branch,color=four), size=2.5, geom='label',node = 'internal', nudge_y = 0.7) +
  ggtree::geom_tiplab(aes(color=four),offset=0.1,size=3) +  xlim(0,15) +
  #ggtree::geom_tippoint(aes(size=100*f_orthologs)) +
  scale_color_metro_d() + theme(legend.position = 'none')
ggsave(MAMMALS_4, filename=file.path(path_ortho,'MAMMALS_TREE_4phylum.pdf'), scale=2)


# 2.1 query Ensembl biomart for human orthologuous proteins in mammals ----------

#nodeid(ens_vertebrates_tree,unique(df_query$treename) %>% str_to_title()) %>%



sc_closest_two = df_query %>%
  filter(two != '0' ) %>%
  group_by(two,ensp) %>%
  arrange(cds_delta,desc(pid_ortho),tx_delta) %>%
  mutate(rk_two = row_number()) %>%
  ungroup() %>% arrange(ensp)
sc_closest_four = df_query %>%
  filter(two != '0' ) %>%
  group_by(four,ensp) %>%
  arrange(cds_delta,desc(pid_ortho),tx_delta) %>%
  mutate(rk_four = row_number()) %>%
  ungroup() %>% arrange(ensp)

# find prefix
ortho_prefix = df_query %>%
  drop_na() %>%
  mutate(n_human = n_distinct(ensp)) %>%
  group_by(ens_dataset) %>%
  mutate( prefix = hutils::longest_prefix(id_ortho) ) %>%
  #mutate( prefix = str_extract(id_ortho,'^[^0-9]+'), n_human = n_distinct(ensp) ) %>%
  group_by(ens_dataset) %>%
  mutate(n_orthologs = n_distinct(id_ortho),
         #f_orthologs = n_orthologs/n_protein_coding, n_uniprot = n_swissprot+n_trembl,
         f_yeast = n_orthologs/n_yeast,
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

# 4. get 1-to-1 orthologs yeast ------------------------------------------------
HS_BOREOEUTHERIA = ortho_prefix %>%
  filter(f_human >0.78) %>%
  group_by(two) %>% mutate(num_two = sprintf("%s (n=%s)",two,n_distinct(ens_dataset))  ) %>%
  group_by(four) %>% mutate(num_four = sprintf("%s (n=%s)",four,n_distinct(ens_dataset)) ) %>%
  ungroup() %>%
  mutate( num_four = factor(num_four, c('out',sort(unique(num_four)))),
          num_two = factor(num_two, c('out',sort(unique(num_two)))) )

phylums = split(HS_BOREOEUTHERIA, HS_BOREOEUTHERIA$two) %>% append(split(HS_BOREOEUTHERIA, HS_BOREOEUTHERIA$four)) %>% compact

hs_ortho_1to1 = HS_ORTHO %>%
  dplyr::filter(!is.na(ensg) & !is.dup(ensp) ) %>%
  group_by(ensg,enst,ensp) %>%
  mutate(n_ortholog = sum(!is.na(c_across(ends_with('homolog_ensembl_peptide'))))) %>%
  filter(n_ortholog == max(n_ortholog) & n_ortholog != 0 )

hs_ortho = hs_ortho_1to1 %>%
  dplyr::select( all_of(colnames(hs_tx)), 'n_ortholog', ends_with('_homolog_ensembl_peptide')) %>%
  rowwise() %>%
  mutate( f_orthogroup =  mean(!is.na(c_across(HS_BOREOEUTHERIA$col_peptide))) ) %>%
  mutate(n_glires = sum.na(c_across(all_of(phylums$Glires$col_peptide)),notNA=T),
         n_carnivora = sum.na(c_across(all_of(phylums$Carnivora$col_peptide)),notNA=T),
         n_laura.1 =sum.na(c_across(all_of(phylums$Laurasiatheria.1$col_peptide)),notNA=T),
         n_primates = sum.na(c_across(all_of(phylums$Primates$col_peptide)),notNA=T),
         n_laurasiatheria =sum.na(c_across(all_of(phylums$Laurasiatheria$col_peptide)),notNA=T),
         n_euarchontoglires = sum.na(c_across(all_of(phylums$Euarchontoglires$col_peptide)),notNA=T) ) %>%
  # Best orthogroups have less than 4 of the selected species (over 78% orthologs w/r to human)
  mutate( is_best = sum(is.na(c_across(HS_BOREOEUTHERIA$col_peptide))) < 3 ) %>%
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
p3 = ggplot(HS_BOREOEUTHERIA, aes(y=reorder(organism,f_human),x=f_human, fill=num_two)) +
  geom_point(shape=21,size=3,color=BEST_COL) + ylab('Mammals Orthologuous species') +
  theme( axis.text.y = element_text(size=9), axis.title = element_text(size=10)) +
  scale_color_metro_d() + scale_fill_metro_d()


p2.1 = ggplot(ortho_prefix, aes(y=reorder(organism,f_human),x=f_human, fill=num_four)) +
  geom_point(shape=21,size=2) + ylab('Mammals Orthologuous species') +
  theme( axis.text.y = element_text(size=9), axis.title = element_text(size=10))+
  geom_vline(xintercept=0.78, linetype='dashed') +
  scale_color_metro_d() + scale_fill_metro_d()
p3.1 = ggplot(HS_BOREOEUTHERIA, aes(y=reorder(organism,f_human),x=f_human, fill=num_four)) +
  geom_point(shape=21,size=3,color=BEST_COL) + ylab('Mammals Orthologuous species') +
  theme( axis.text.y = element_text(size=9), axis.title = element_text(size=10)) +
  scale_color_metro_d() + scale_fill_metro_d()

library(patchwork)
graphics.off()
MAMMALS_2 = MAMMALS_2 +
  geom_tippoint(data = . %>% filter(label %in% HS_BOREOEUTHERIA$num_label), color=BEST_COL,size=3)

MAMMALS_4 = MAMMALS_4 +
  geom_tippoint(data = . %>% filter(label %in% HS_BOREOEUTHERIA$num_label), color=BEST_COL,size=3)

PP = ((p1) | (MAMMALS_2)) / (p2 | p3)
ggplot2::ggsave(PP, height=12,width=18, filename = file.path(path_ortho,'Mammals-orthologs-2_phylum.pdf'))
PP.1 = ((p1) | (MAMMALS_4)) / (p2.1 | p3.1)
ggplot2::ggsave(PP.1, height=12,width=18, filename = file.path(path_ortho,'Mammals-orthologs-4_phylum.pdf'))

# 7. Write the fasta for mammals and subphylum ---------------------------------
hs_ortho_1to1.rds=file.path(path_ortho,'hs_ortho_1to1.rds')
hs_ortho.rds=file.path(path_ortho,'hs_orthologs.rds')
ortho_prefix.rds = file.path(path_ortho,'ensembl_species.rds')
HS_BOREOEUTHERIA.rds=file.path(path_ortho,'ensembl_mammals.rds')
hs_fasta.rds =file.path(path_ortho,'ensembl_human_seq.rds')
hs_closest.rds = file.path(path_ortho,'hs_closest_ortho.rds')

#saveRDS(hs_ortho_1to1, hs_ortho_1to1.rds)
#saveRDS(hs_ortho, hs_ortho.rds)
#saveRDS(ortho_prefix, ortho_prefix.rds)
#saveRDS(HS_BOREOEUTHERIA, HS_BOREOEUTHERIA.rds)
#saveRDS(hs_fasta_prot, hs_fasta.rds)
#saveRDS(list(two=hs_closest_two,four=hs_closest_four),hs_closest.rds)

#hs_ortho_1to1=readRDS(hs_ortho_1to1.rds)
#HS_BOREOEUTHERIA=readRDS(HS_BOREOEUTHERIA.rds)
#ortho_prefix=readRDS(ortho_prefix.rds)
#hs_closest_two = readRDS(hs_closest.rds)$two
#hs_closest_four = readRDS(hs_closest.rds)$four

# get all sequences
mammals_seq=readRDS(file.path(path_ortho,'ensembl_mammals_seq.rds'))
all_seq=readRDS(hs_fasta.rds)
for( i in 1:length(mammals_seq) ){ all_seq = AAStringSet(c(all_seq,mammals_seq[[i]])) }

mammals0.path  = file.path(path_ortho,'hs_mammals_all')
mammals.path   = file.path(path_ortho,'hs_mammals')
glires.path    = file.path(path_ortho,'hs_glires')
laura1.path     = file.path(path_ortho,'hs_laurasiatheria.1')
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





