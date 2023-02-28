library(tidyverse)
library(here)
source(here::here("src","__setup_yeastomics__.r"))
library(treeio)
library(tidytree)
ncbi_dir = here::here("data","ncbi")
#eggnog_fasta = Biostrings::readAAStringSet("http://eggnog5.embl.de/download/latest/e5.proteomes.faa")

rename.tips <- function(phy, old_names, new_names) {
  mpos <- match(old_names,phy$tip.label)
  phy$tip.label[mpos] <- new_names
  return(phy)
}

reduce_tree = function(treein, fasta_1to1){
  ids = names(fasta_1to1)
  taxons = ids %>% str_extract("^[0-9]+")
  reduced_tree = keep.tip(treein,taxons)
  return(rename.tips(reduced_tree, taxons, ids))
}

remove_outliers = function(BS,max_fX=0.1){
  BS_normalized = normalize_sequence(BS)
  df_aa =  count_aa(BS_normalized) %>%
           ungroup() %>%
           mutate( is_outlier = (f_noAA >= max_fX) ) %>%
           filter( !is_outlier )
  n_seq = n_distinct(names(BS))
  n_outliers = n_seq - n_distinct(df_aa$ids)
  cat(sprintf("removed %s outliers out of %s sequences\n",n_outliers,n_seq))
  return(BS_normalized[df_aa$ids])
}

reduce_orthologs = function(fastaseq, id_ref, id_orthogroup){

  nref = length(id_ref)
  fastaseq_norm = normalize_sequence(fastaseq)
  # VALID ORTHOLOG = present in the fasta files
  id_orthologs = intersect(names(fastaseq_norm), id_orthogroup)

  reference_seq = fastaseq_norm[id_ref] %>% as.character()
  orthogroup_seq = fastaseq_norm[id_orthologs] %>% as.character()

  df_seq = tibble(ids=id_orthologs) %>%
           separate(col=ids, into=c('taxon','string'), sep = "\\.", extra ='merge', remove = F) %>%
           group_by(taxon) %>%
           add_count(name='np')

  is_one2one = df_seq$np==1
  id_one2one = df_seq$ids[is_one2one]

  is_dup = df_seq$np>1
  id_dup = df_seq$ids[is_dup]


  DF_CLOSEST = map_dfr(df_seq$ids,
                    ~score_ali(p1=orthogroup_seq[.x], s1=.x,p2=reference_seq[1], s2=id_ref[1])) %>%
      mutate(taxon1=str_split_fixed(ID1,'\\.',2)[,1] ) %>%
      mutate(delta_ol = abs(100-OL1), delta_pid = abs(100-PID1) ) %>%
      group_by(taxon1) %>%
      arrange( taxon1, delta_ol, delta_pid, desc(SCORE.BLOSUM62),desc(S), .by_group = TRUE )

  CLOSEST= DF_CLOSEST %>%
      dplyr::slice(1,.preserve =T) %>%
      ungroup() %>%
      mutate( has_dup = ID1 %in% id_dup, is_ref = ID1 == ID2) %>%
      arrange( desc(is_ref), delta_ol, delta_pid )

  if(nref==2){
    DF_CLOSEST_2 = map_dfr(df_seq$ids,
                         ~score_ali(p1=orthogroup_seq[.x], s1=.x,p2=reference_seq[2], s2=id_ref[2])) %>%
      mutate(taxon1=str_split_fixed(ID1,'\\.',2)[,1] ) %>%
      mutate(delta_ol = abs(100-OL1), delta_pid = abs(100-PID1) ) %>%
      group_by(taxon1) %>%
      arrange( taxon1, delta_ol, delta_pid, desc(SCORE.BLOSUM62),desc(S), .by_group = TRUE )

    CLOSEST2= DF_CLOSEST_2 %>%
          anti_join(CLOSEST %>% filter(has_dup), by=c('ID1')) %>%  # NOT REUSING THE 1ST DUPLICATES
          dplyr::slice(1,.preserve =T) %>%
          ungroup() %>%
          mutate( has_dup = ID1 %in% id_dup, is_ref = ID1 == ID2) %>%
          arrange( desc(is_ref), delta_ol, delta_pid )

    return(list(fastaseq_norm[CLOSEST$ID1],fastaseq_norm[CLOSEST2$ID1]))
  }
  return(list(fastaseq_norm[CLOSEST$ID1]))
}

open_fasta = function(fastafile, fastaurl,debug=F, max_fX = 0.1){
  if(debug){
    .dbg$log(sprintf('get fasta sequences for the full orthogroup (%s proteins)...',NPROT))
  }

  if( file.exists(fastafile) ){
    return( Biostrings::readAAStringSet(fastafile) )
  }else{
    # Get fasta sequence of the full orthogroup
    temp = tempfile()
    cat('\rdownloading fasta from url...             ')
    download.file(fastaurl,destfile = temp,quiet = T)
    Sys.sleep(3)
    return( Biostrings::readAAStringSet(temp) )
  }
}

check_ortho = function(outfasta,outtree,debug=F){

  has_fasta = file.exists( file.path(outfasta) )
  has_tree = file.exists( file.path(outtree) )
  SEQ=if(has_fasta) Biostrings::readAAStringSet(file.path(outfasta)) else NULL
  TREE=if(has_tree) ape::read.tree(file.path(outtree)) else NULL
  return( list(fasta=SEQ,tree=TREE) )
}

write_r4s_input = function(tree,fasta,outfasta,outtree,debug=F){

  if(debug){
    .dbg$log('...writing the tree and alignments for 1-to-1 orthologs...')
  }

  fastaname = basename(outfasta)
  treename = basename(outtree)

  np_tree  = n_distinct(tree$tip.label)
  np_fasta = n_distinct(names(fasta))

  if(np_tree != np_fasta){
    .error$log(sprintf('Not the same number of 1-to-1 orthologs in tree (%s) and alignment (%s)',np_tree,np_fasta))
  }else{
    if(np_tree<3){
      .error$log(sprintf("not enough 1-to-1 orthologs found (%s)...",np_tree))
    }else{

      if(!file.exists(outfasta)){
        .dbg$log(sprintf("writing fasta file %s",fastaname))
        writeXStringSet(fasta,outfasta,append = F,format = 'fasta')
      }

      #if(!file.exists(outtree)){
      if(debug){ .dbg$log(sprintf("writing tree file %s",treename)) }
      #tree$edge.length=NULL
      tree$node.label=NULL
      #tree$root.edge=NULL
      #print(tree)
      write.tree(tree,outtree)
    }
  }
}

filter_orthogroups = function( orthologs_count, ref_sp = 4932, ref_tree,
                               debug=F, force=F, eggnog_tree=F){

  if(missing(orthologs_count)){
    .error$log('Must input a dataframe with the orthogroups from a taxonomic level across valid subclades...')
    return(NULL)
    #node = find_eggnog_node()
    #clades = get_eggnog_taxonomy(node$id,.print = T) %>% mutate(clade_desc = sprintf("%s_%s (n=%s)",clade_id,clade_name,clade_size))
    #selected = select.list(clades$clade_desc,multiple = T)
    #orthologs_count = count_taxons_eggnog_node(node$id, clades$clade_id[ clades$clade_desc %in% selected] ) %>% bind_rows()
  }

  .info$log("remove orthogroups without reference species...")
  og_ref = orthologs_count
  nstep = nrow(og_ref)
  pc1 = as.integer(nstep/ min(nstep,1000))

  # 3NTXF.2
  #i = which(og_ref$OG == "3NVWA")#3NTWU # 3NU1J") # test outlier

  for( i in 1:nrow(og_ref) ){

    og_data = og_ref[i,]
    perc = 100*(i/nstep)
    if(i %% pc1 == 0){ cat(sprintf("[%.1f%%] ---> %6d/%6d      \r",perc,i,nstep)) }
    OG = og_data$OG
    NPROT = og_data$clade_np

    NODE_ID  = og_data$node_id
    NODE_NAME = og_data$node_name
    NODE  = og_data$node_orthologs[[1]]
    FASTA_DIR = here::here('data','eggnog','fasta',paste0(NODE_ID,'-',NODE_NAME) %>% tolower)

    REF_ID = NODE$id[ NODE$taxid == ref_sp ]
    nref = length(REF_ID)

    if(nref==0){
      .error$log(sprintf('(i=%s) orthogroup %s has no reference sequence (taxid=%s)',i,OG,ref_sp))
      next
    }else if( nref>2 ){
      if(debug){
        .warn$log(sprintf('(i=%s) orthogroup %s has multiple reference sequences (taxid=%s n=%s)',i,OG,ref_sp,nref))
      }
      next
    }

    CLADE_FNAME = og_data$clade_dirname
    CLADE_ID = og_data$clade_id
    CLADE_SIZE = og_data$clade_size
    CLADE_N = og_data$clade_n
    CLADE_NS = og_data$clade_ns
    CLADE_FS = og_data$clade_f
    is_one2one = og_data$clade_one2one

    CLADE = og_data$clade_orthologs[[1]] %>%
             mutate(ref = (taxid==ref_sp), is_dup = is.dup(taxid) ) %>%
             arrange(desc(ref),desc(is_dup))
    CLADE_NOREF = CLADE |> filter(!ref)

    NODE_DESC =  og_data$node_desc
    CLADE_DESC = og_data$clade_desc
    CLADE_DIR = here::here('data','eggnog',NODE_DESC,CLADE_DESC)
    dir.create(CLADE_DIR,showWarnings = F,recursive = T)
    dir.create(file.path(CLADE_DIR,'fasta'),showWarnings = F,recursive = T)
    dir.create(file.path(CLADE_DIR,'tree'),showWarnings = F,recursive = T)

    og_fastafile = file.path(FASTA_DIR,paste0(OG,'.fasta'))
    fasta_in = open_fasta(og_fastafile, og_data$url_fasta) %>%
               remove_outliers(max_fX = 0.1)

    OG_unique = rep(OG, times=nref) |> makeunique::make_unique(sep='.',wrap_in_brackets = F)
    fasta_out = sprintf("%s_%s_%s_1to1-orthologs.fa",CLADE_FNAME,OG_unique,REF_ID)
    tree_out =  str_replace(fasta_out,'\\.fa$','.nwk')
    if(eggnog_tree){
      tree_out = str_replace(fasta_out,'\\.fa$','.nwk_eggnog') # FOR SPECIES TREE
    }
    fastapath = file.path(CLADE_DIR,"fasta",fasta_out)
    treepath = file.path(CLADE_DIR,"tree",tree_out)
    has_output_fasta = file.exists(fastapath)
    has_output_tree = file.exists(treepath)

    if(nref==2){
      ogs = paste0(OG_unique,collapse=" ")
      if(debug){
        .dbg$log(sprintf('(i=%s) orthogroup %s has %s reference orthologs (to be splitted into: %s)',i,OG,nref,ogs))
      }
    }

    if(!force & all(has_output_fasta) & all(has_output_tree)){
      next
    }

    fasta_1to1 = lapply(fastapath[has_output_fasta],Biostrings::readAAStringSet)
    tree_1to1 = lapply(treepath[has_output_tree],ape::read.tree)
    if(force | (length(fasta_1to1) != nref) ){
      fasta_1to1 = reduce_orthologs(fastaseq = fasta_in,
                                    id_ref = REF_ID,
                                    id_orthogroup = CLADE$id)
    }

    if(force | (length(tree_1to1) != nref) ){
      tree_1to1 = lapply(fasta_1to1,function(fa){ reduce_tree(ref_tree,fa) })
    }

    for(R in 1:nref){
      write_r4s_input(tree = tree_1to1[[R]],  fasta = fasta_1to1[[R]],
                      outfasta = fastapath[R], outtree=treepath[R])
    }
  }
}

fungi_clades = c('4890'='ascomycota', '5204'='basidiomycota')
ascomycota_clades = c('451866'='taphrimycotina','4891'='saccharomycetes','147541'='dothideomycetes',
                      '147545'='eurotiomycetes','147550'='sordariomycetes','147548'='leotiomycetes')
#### FUNGI ####
# Fungi = 4751
# |--Basidiomycota = 5204
# |--Ascomycota = 4890
# |--|-----Taphrimycotina = 451866 ----> S. pombe 4896
# |---|----Saccharomycetes = 4891 ---> S. cerevisiae 4932
# |----|--Dothideomycetes = 147541
# |----|---Eurotiomycetes = 147545
# |-----|--Sordariomycetes = 147550
# |-----|--Leotiomycetes = 147548

fu_dir = here::here("data","eggnog","4751_Fungi")
fu_node    = find_eggnog_node(node = 4751,.print = F) %>% dplyr::rename_with(~Pxx(.,px='node',s="_"))
fuNOG      = get_eggnog_node(node = 4751)
fu_tax     = get_eggnog_taxonomy(4751)
fu_species = get_eggnog_species(node = 4751)

####_download NCBI common tree from list of species (manual) ####
# write_lines(fu_species$taxid, file.path(ncbi_dir,'4751-fungi-179taxids.txt'))
# write_lines(fu_species$taxon, file.path(ncbi_dir,'4751-fungi-179taxons.txt'))
fungi = treeio::read.tree(file.path(ncbi_dir,"ncbi-fungi.phy"))
fungi$tip.label = str_remove_all(fungi$tip.label,"['\\[\\]]") #%>% str_replace_all(" ","_")

fu_ncbi = preload( saved.file = file.path(ncbi_dir,'ncbi-to-eggnog-fungi-179species.rds'),
                    { match_strings(fungi$tip.label, SP2=fu_species$taxon, use_soundex = F, manual = T) },
                    'match ncbi species tree to eggnog fungal species...')
fu_ncbi_eggnog = fu_ncbi %>%
   left_join( fu_species %>% mutate(Taxon=taxon, taxon=tolower(taxon)), c('s2'='taxon')) %>%
   dplyr::rename(ncbi_name=s1,eggnog_name=s2) %>%
   dplyr::select(-c(is_identical:is_substring,osa:n1))

####_get NCBI tree from taxon identifiers (not necessarily find all ids) ####
#fu_ncbi = get_ncbi_tree(fu_species$taxid)
#fungi_eggnog$tip.label = fu_ncbi_eggnog$taxid

library(taxize)
fungi$tip.label = fu_ncbi_eggnog$taxid[match(tolower(fungi$tip.label),fu_ncbi_eggnog$ncbi_name)]
# fungi_19$taxids = fu_ncbi_eggnog$taxid[match(tolower(fungi_19$tip.label),fu_ncbi_eggnog$taxid)]
# fungi_19$taxons = fu_ncbi_eggnog$Taxon[match(tolower(fungi_19$tip.label),fu_ncbi_eggnog$taxid)]
# fungi_19$tip.label = fungi_19$taxons

fu_clades = subtrees(fungi) %>%
            set_names(fungi$node.label %>% str_remove_all("['\\[\\]]"))  %>%
            map_dfr( ~ tibble(clade_sp = list(.x$tip.label %>% sort %>% factor(.,levels=fu_species$taxid)), clade_size= n_distinct(unlist(clade_sp))), .id = "clade_name")  %>%
            filter(clade_size > 2 & clade_name != "") %>%
            bind_cols( fu_node ) %>%
            mutate( clade_id = taxize::get_ids(clade_name, db='ncbi',verbose = F)$ncbi %>% as.character(),
                    clade_dirname = clade_name %>% str_replace_all("[^a-zA-Z0-9_\\-]+","."),
                    node_dirname = node_name %>% str_replace_all("[^a-zA-Z0-9_\\-]+","."),
                    is_clade = T, is_eggnog=F,
                    is_subnode = clade_size < node_size | clade_size == node_id,
                    clade_desc =sprintf('%s_%s_%ssp',clade_id,clade_dirname,clade_size),
                    node_desc =sprintf('%s_%s',node_id,node_dirname),
                    clade_has_4896.4932 = map(clade_sp, ~ any(.x %in% c(4932,4896))) %>% unlist) %>%
            relocate(node_id,node_name,node_dirname,node_desc,node_size,
                     clade_id,clade_name,clade_dirname,clade_desc,clade_size,
                     is_clade,is_subnode,clade_sp) %>%
            arrange(desc(clade_size))

fu_clades_toprocess = fu_clades %>% filter(clade_size > 20 | clade_has_4896.4932)
fu_yeast   = eggnog_annotations_species(node = 4751, species = c(4932,4896))
fu_og  = preload(saved.file = here::here("data/eggnog/4751_Fungi-orthogroups.rds"),
                 loading.call = { count_eggnog_orthologs(4751) },
                 doing = "fungi othogroups...")

fu_orthogroups  = preload(saved.file = here::here("data/eggnog/4751_Fungi-orthogroups-clades.rds"),
                          loading.call = { count_clade_orthologs(fu_og, fu_clades_toprocess) },
                          doing = "orthogroups by fungi subclades...")

fu_ref = fu_orthogroups %>%
  bind_rows() %>%
  mutate(ref_sp = 4932,
         # remove orthogroup from fungi node present in other subnodes
         node_only = !( clade_id == 4751 & OG %in% OG[clade_id != 4751] ),
         node_has_ref = map(node_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist,
         clade_has_ref = map(clade_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist,
         node_has_4896 = map(node_orthologs, ~ sum(.x$taxid %in% 4896)) %>% unlist,
         clade_has_4896 = map(clade_orthologs, ~ sum(.x$taxid %in% 4896)) %>% unlist)


fu_toprocess = fu_ref %>%
  filter(node_has_ref %in% c(1,2) & (clade_ns > 20 | clade_has_ref %in% c(0,1,2) | clade_has_4896 %in% c(1,2)) & clade_f > 0.5 ) %>%
  group_by(node_id,clade_id) %>% add_count(name='clade_n')

dim(fu_orthogroups)
dim(fu_ref)
dim(fu_toprocess)
table(fu_toprocess$clade_name,fu_toprocess$node_has_ref)
janitor::tabyl(fu_toprocess,clade_name,clade_ns)

####_filter orthogroups fasta to keep at most 179 species (ortholog closest to yeast) ####

#chunk_size=100
#chunks=seq(1,nrow(fu_toprocess),by=chunk_size)
#list_fu_toprocess = split(fu_toprocess, fu_toprocess$clade_desc)
#chunks=seq(1,length(list_fu_toprocess),by=1)
# pbmcapply::pbmclapply(chunks,FUN = function(og){
#   filter_orthogroups(
#     orthologs_count = list_fu_toprocess[[og]],
#     ref_tree = fungi_19, eggnog_tree = T,
#     ref_sp = 4932,
#     debug = F,
#     force = F) },mc.cores = 10)

fu_sptree=here::here(fu_dir,"fungi_sptree","cog_100-alg_concat_default-raxml_default")
fungi_19 = treeio::read.tree(file.path(fu_sptree,"fungi179-19orthogroups-seqs-sorted.fa.gz.final_tree.nw"))


fu_toprocess_left = fu_toprocess %>% filter( clade_name %in% c('Fungi','Ascomycota','Dikarya') )

list_fu_toprocess =  split( fu_toprocess_left , cut_width(1:nrow(fu_toprocess_left), width=1,boundary = 0) )
chunks=seq(1,length(list_fu_toprocess),by=1)

#tmp = fungi_outliers %>% bind_rows() %>% filter(f_noAA > 0.1) %>% arrange(desc(f_noAA))
#OG_DATA = fu_toprocess %>% filter(clade_name == "Pezizomycotina" )

library(furrr)
plan(multisession, workers = 14)
furrr::future_map(list_fu_toprocess, function(og){
  filter_orthogroups( orthologs_count = og,
                      ref_tree = fungi_19, eggnog_tree = T, ref_sp = 4932,
                      debug = F, force = F) },.progress = T)
#
# pbmcapply::pbmclapply(chunks,FUN = function(irow){
#                      filter_orthogroups(
#                           orthologs_count = fu_toprocess[irow:(irow+chunk_size),],
#                           ref_tree = fungi_19,
#                           ref_sp = 4932,
#                           debug = F,
#                           force = F) },mc.cores = 14)

#### _retrieveing orthogroups to build species tree ####
# Select orthogroups for building the fungi species tree
# ---> orthogroups must be at the taxonomic level of fungi (id=4751)
# ---> all species must be present (179 fungal species)
# ---> at least 95% of orthologs should be 1-to-1
# ---> s. cerevisiae ortholog must be single (not duplicated)
fu_sptree = fu_ref %>%
            filter(clade_id == node_id & clade_f1to1>0.975 & clade_f == 1 & node_has_ref == 1) %>%
            mutate( ref_id = map_chr(node_orthogroup,~str_subset(.x, "^4932\\.")) )

fu_sptree_fasta = with(fu_sptree,
                       file.path(fu_dir,clade_desc,"fasta",
                                 paste0(clade_name,"_",OG,"_",ref_id,"_","1to1-orthologs.fa")))
fu_sptree_tree = with(fu_sptree,
                      file.path(fu_dir,clade_desc,"tree",
                                paste0(clade_name,"_",OG,"_",ref_id,"_","1to1-orthologs.nwk")))

fu_sptree_dir = here::here("data","eggnog","4751_Fungi","FUNGI-SPECIES-TREE")
#file.exists(fu_sptree_fasta)
#file.exists(fu_sptree_tree)
file.copy(fu_sptree_fasta,fu_sptree_dir)
file.copy(fu_sptree_tree,fu_sptree_dir)

orthogroup_fasta = Biostrings::readAAStringSet(fu_sptree_fasta)
names(orthogroup_fasta) = names(orthogroup_fasta)  %>% str_replace("\\.","|")
orthogroup_fasta.sorted = orthogroup_fasta[gtools::mixedsort(names(orthogroup_fasta))]

orthologs_ids =map(fu_sptree_fasta,
                   ~Biostrings::readAAStringSet(.x) %>%
                     names() %>%
                     str_replace("\\.","|") %>%
                     gtools::mixedsort() )
Biostrings::writeXStringSet(orthogroup_fasta.sorted,
                            filepath=file.path(mz_sptree_dir,"fungi179-19orthogroups-seqs-sorted.fa.gz"), compress = T, format = 'fasta')

write_lines(unlist(lapply(orthologs_ids, paste, collapse="\t")),
            file = file.path(mz_sptree_dir,"fungi179-19orthogroups-seqs-sorted.ids"))

#### METAZOA ####
# Metazoa = 33208
# Mammalia = 40674

mz_dir = here::here('data','eggnog','33208_Metazoa')
mz_node    = find_eggnog_node(node = 33208,.print = F) %>% dplyr::rename_with(~Pxx(.,px='node',s="_"))
mzNOG      = get_eggnog_node(node = 33208)
mz_tax     = get_eggnog_taxonomy(33208)
mz_species = get_eggnog_species(node = 33208)

metazoa_clades = c('7742'='vertebrata','40674'='mammalia','50557'='insecta','6231'='nematoda')

####_download NCBI common tree from list of species (manual) ####
# write_lines(mz_species$taxid, here::here('data','ncbi','33208-metazoa-161taxids.txt'))
# write_lines(mz_species$taxon, here::here('data','ncbi','33208-metazoa-161taxons.txt'))
metazoa = treeio::read.tree(here::here('data','ncbi','ncbi-metazoan.phy'))
metazoa$tip.label = str_remove_all(metazoa$tip.label,"['\\[\\]]") #%>% str_replace_all(" ","_")

mz_ncbi = preload( saved.file = file.path(ncbi_dir,'ncbi-to-eggnog-metazoa-161species.rds'),
                   { match_strings(metazoa$tip.label, SP2=mz_species$taxon, use_soundex = F, manual = T) },
                   'match ncbi species tree to eggnog fungal species...')
mz_ncbi_eggnog = mz_ncbi %>%
  left_join( mz_species %>% mutate(Taxon=taxon, taxon=tolower(taxon)), c('s2'='taxon')) %>%
  dplyr::rename(ncbi_name=s1,eggnog_name=s2) %>%
  dplyr::select(-c(is_identical:is_substring,osa:n1)) %>%
  mutate(ncbi_Name = str_to_sentence(ncbi_name))

####_get NCBI tree from taxon identifiers (not necessarily find all ids) ####
#mz_ncbi = get_ncbi_tree(mz_species$taxid)
#metazoa$tip.label = mz_ncbi_eggnog$taxid

library(taxize)
metazoa$tip.label = mz_ncbi_eggnog$taxid[match(tolower(metazoa$tip.label),mz_ncbi_eggnog$ncbi_name)]
mz_clades = subtrees(metazoa) %>%
  set_names(metazoa$node.label %>% str_remove_all("['\\[\\]]"))  %>%
  map_dfr( ~ tibble(clade_sp = list(.x$tip.label %>% sort), clade_size= n_distinct(unlist(clade_sp))), .id = "clade_name")  %>%
  filter(clade_size > 2, clade_name != "") %>%
  bind_cols( mz_node ) %>%
  mutate( clade_id = taxize::get_ids(clade_name, db='ncbi',verbose = F)$ncbi %>% as.character(),
          clade_dirname = clade_name %>% str_replace_all("[^a-zA-Z0-9_\\-]+","."),
          node_dirname = node_name %>% str_replace_all("[^a-zA-Z0-9_\\-]+","."),
          is_clade = T, is_eggnog=F,
          is_subnode = clade_size < node_size | clade_size == node_id,
          clade_desc =sprintf('%s_%s_%ssp',clade_id,clade_dirname,clade_size),
          node_desc =sprintf('%s_%s',node_id,node_dirname),
          clade_has_9606 = map(clade_sp, ~ any(.x %in% c(9606))) %>% unlist) %>%
  relocate(node_id,node_name,node_dirname,node_desc,node_size,
           clade_id,clade_name,clade_dirname,clade_desc,clade_size,
           is_clade,is_subnode,clade_sp) %>%
  arrange(desc(clade_size))

mz_clades_keep=c(Metazoa = 33208, Bilateria = 33213,
                 Protostomia = 33317,  Deuterostomia = 33511,
                 Endopterygota = 33392, Nematoda = 6231,
                 Diptera = 7147, Drosophila = 32281,
                 Neopterygii = 41665, Tetrapoda = 32523,
                 Sauria = 32561, Mammalia = 40674,
                 Eutheria = 9347, Laurasiatheria = 314145,
                 Boreoeutheria = 1437010, Euarchontoglires = 314146,
                 Afrotheria = 311790, Chiroptera = 9397,
                 Carnivora = 33554, Artiodactyla = 91561,
                 Primates = 9443, Hominoidea = 314295, Glires = 314147)

mz_clades_toprocess = mz_clades %>% filter(clade_size > 5 | clade_has_9606) %>%
  filter(clade_id %in% mz_clades_keep)

mz_human   = eggnog_annotations_species(node = 33208, species = c(9606))

mz_og  = preload(saved.file = here::here("data/eggnog/33208_Metazoa-orthogroups.rds"),
                 loading.call = { count_eggnog_orthologs(33208) },
                 doing = "metazoa orthogroups...")

mz_orthogroups  = preload(saved.file = here::here("data/eggnog/33208_Metazoa-orthogroups-clades_keep.rds"),
                          loading.call = { count_clade_orthologs(mz_og, mz_clades_toprocess) },
                          doing = "orthogroups by metazoa subclades...")

mz_ref = mz_orthogroups %>%
  bind_rows() %>%
  mutate(ref_sp = 9606,
         # remove orthogroup from metazoa node present in other subnodes
         node_only = !( clade_id == 33208 & OG %in% OG[clade_id != 33208] ),
         node_has_ref = map(node_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist,
         clade_has_ref = map(clade_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist)

mz_toprocess = mz_ref %>%
  filter(node_has_ref %in% c(1,2) & (clade_ns > 5 | clade_has_ref %in% c(1,2)) & clade_f > 0.5 ) %>%
  group_by(node_id,clade_id) %>% add_count(name='clade_n')

dim(mz_orthogroups)
dim(mz_ref)
dim(mz_toprocess)
table(mz_toprocess$clade_name,mz_toprocess$node_has_ref)
janitor::tabyl(mz_toprocess,clade_name,clade_ns)

####_filter orthogroups fasta to keep at most 161 species (ortholog closest to human) ####
#chunk_size=100
#chunks=seq(1,nrow(mz_toprocess),by=chunk_size)
#mz_toprocess.2 = mz_toprocess %>% filter(clade_ns < 20 & clade_has_ref == 0)
#list_mz_toprocess = split(mz_toprocess, mz_toprocess$clade_desc)

mz_sptree=here::here(mz_dir,"metazoa_sptree","cog_100-alg_concat_default-raxml_default")
metazoa_15 = treeio::read.tree(file.path(mz_sptree,"metazoa161-15orthogroups-seqs-sorted.fa.gz.final_tree.nw"))

piority_clade = c("Mammalia","Deuterostomia","Tetrapoda","Bilateria","Metazoa")
mz_process_priority = mz_toprocess %>% filter( clade_id %in% mz_clades_keep[piority_clade])
list_mz_toprocess =  split( mz_process_priority , cut_width(1:nrow(mz_process_priority), width=1,boundary = 0) )
chunks=seq(1,length(list_mz_toprocess),by=1)

# ## TREE MUST HAVE TAXON ID AS LEAVES
# pbmcapply::pbmclapply(chunks,FUN = function(og){
#   filter_orthogroups(
#     orthologs_count = list_mz_toprocess[[og]],
#     ref_tree = metazoa_15, eggnog_tree = T,
#     ref_sp = 9606,
#     debug = F,
#     force = F) },mc.cores = 10)

library(furrr)
og = mz_toprocess %>% filter(OG == "3B94D" & clade_id == mz_clades_keep['Boreoeutheria'])
plan(multisession, workers = 14)
furrr::future_map(list_mz_toprocess, function(og){
  filter_orthogroups(
    orthologs_count = og,
    ref_tree = metazoa_15, eggnog_tree = T,ref_sp = 9606,
    debug = F,force = T)},
.progress = T)



#### _retrieveing orthogroups to build species tree ####
# Select orthogroups for building the metazoa species tree
# ---> orthogroups must be at the taxonomic level of fungi (id=33208)
# ---> all species must be present (161 metazoan species)
# ---> at least 95% of orthologs should be 1-to-1
# ---> h. sapiens ortholog must be single (not duplicated)
mz_sptree = mz_ref %>%
  filter(clade_id == mz_node$node_id & clade_f1to1>0.95 & clade_f == 1 & node_has_ref == 1) %>%
  mutate( ref_id = map_chr(node_orthogroup,~str_subset(.x, "^9606\\.")) )

mz_sptree_fasta = with(mz_sptree,file.path(mz_dir,clade_desc,"fasta",
                         paste0(clade_name,"_",OG,"_",ref_id,"_","1to1-orthologs.fa")))
mz_sptree_tree = with(mz_sptree,file.path(mz_dir,clade_desc,"tree",
                                           paste0(clade_name,"_",OG,"_",ref_id,"_","1to1-orthologs.nwk")))
mz_sptree_dir = here::here("data","eggnog","33208_Metazoa_speciestree")
#file.exists(mz_sptree_fasta)
#file.exists(mz_sptree_tree)
file.copy(mz_sptree_fasta,mz_sptree_dir)
file.copy(mz_sptree_tree,mz_sptree_dir)

orthogroup_fasta = Biostrings::readAAStringSet(mz_sptree_fasta)
names(orthogroup_fasta) = names(orthogroup_fasta)  %>% str_replace("\\.","|")
orthogroup_fasta.sorted = orthogroup_fasta[gtools::mixedsort(names(orthogroup_fasta))]


orthologs_ids =map(mz_sptree_fasta,
                   ~Biostrings::readAAStringSet(.x) %>%
                     names() %>%
                     str_replace("\\.","|") %>%
                     gtools::mixedsort() )
Biostrings::writeXStringSet(orthogroup_fasta.sorted,
                            filepath=file.path(mz_sptree_dir,"metazoa161-15orthogroups-seqs-sorted.fa.gz"), compress = T, format = 'fasta')

write_lines(unlist(lapply(orthologs_ids, paste, collapse="\t")),
            file = file.path(mz_sptree_dir,"metazoa161-15orthogroups.ids"))

# laurasiatheria = extract.clade(mammalia_eggnog, 'Laurasiatheria', root.edge = 0, collapse.singles = TRUE)
# laurasiatheria_clades=c('cetartiodactyla'=91561,'carnivora'=33554,'chiroptera'=9397)
# id_laurasiatheria = laurasiatheria$tip.label %>% as.numeric %>% sum
#
#pbmcapply::pbmclapply(mzNOG$url_fasta, safe_download, mc.cores=14, path=mz_dir, ext='.fasta')

#### EGGNOG V6 ####
#"http://eggnog6.embl.de/download/eggnog_6.0/e6.og2level.tsv"
#"http://eggnog6.embl.de/download/eggnog_6.0/e6.og2seqs_and_species.tsv"
# e6_dir = here::here('data','eggnog','eggnog_v6')
# ncbi_dir = here::here('data','ncbi')
# e6.tax = readr::read_delim(file.path(e6_dir,"e6.og2level.tsv"),
#                            delim = '\t',col_names = c("node","OG"))
#
# e6.species = readr::read_delim(file.path(e6_dir,"e6.og2seqs_and_species.tsv"),
#                                delim = '\t',col_types = 'cciicc',
#                                col_names = c("node","OG","nsp","nprot","taxids","string_ids"))
#
# e6.ogs = readr::read_delim(file.path(e6_dir,"e6.seq2ogs.tsv"), delim = '\t',
#                            col_types = 'cc',
#                            col_names = c("string_id","OGs"))
#
# e6.hierarchy = readr::read_delim(file.path(e6_dir,"e6.og2parents_and_children.new.tsv"),
#                                  delim = '\t', col_types = 'ciicc',
#                                  col_names = c("OG","nparent","nchild","OG_parent","OG_child"))
#

#### NCBI
# ncbi_names = readr::read_delim(file.path(ncbi_dir,'names.dmp'),
#                                col_names=c('taxid','name','uname','class'),
#                                delim='\t|\t', col_types='cccc', col_select = 1:4) %>%
#   mutate(class = str_replace(class,'\t\\|',""))
#
# ncbi_sciname = ncbi_names %>% filter(class == 'scientific name')
# library(hutils)
# ncbi_sciname %>% filter( name %pin% fungi_eggnog$node.label)
#
#
# node_pnames = ncbi_sciname %>% filter( name %pin% fungi_eggnog$node.label) %>% pull(name)
# tmp = match_strings(fungi_eggnog$node.label, node_pnames, manual = T)

# ncbi_node = readr::read_delim(file.path(ncbi_dir,'nodes.dmp'),
#                               col_names=c('taxid','taxid_parent','rank','embl_code'),
#                               delim='\t|\t', col_types='cccccccccc', col_select=1:13)

# ncbi_lineage = readr::read_delim(file.path(ncbi_dir,'rankedlineage.dmp'),
#                                  col_names=c(), delim='\t|\t', col_types='cccc', col_select=1:5)
