library(tidyverse)
library(here)
source(here::here("src","__setup_yeastomics__.r"))
library(treeio)
library(tidytree)
ncbi_dir = here::here("data","ncbi")
#eggnog_fasta = Biostrings::readAAStringSet("http://eggnog5.embl.de/download/latest/e5.proteomes.faa")

reduce_orthologs = function(fastaseq, id_ref, id_orthogroup){

  reference_seq = fastaseq[id_ref] %>% as.character() %>% chartr("U","X",x=.)
  orthogroup_seq = fastaseq[id_orthogroup] %>% as.character() %>% chartr("U","X",x=.)
  df_seq = tibble(ids = id_orthogroup) %>%
           separate(col=ids, into=c('taxon','string'), sep = "\\.", extra ='merge', remove = F) %>%
           group_by(taxon) %>%
           add_count(name='np')

  is_one2one = df_seq$np==1
  id_one2one = df_seq$ids[is_one2one]

  is_dup = df_seq$np>1
  id_dup = df_seq$ids[is_dup]

  CLOSEST = map_dfr(seq_along(id_orthogroup), ~score_ali(p1=orthogroup_seq[.x], s1=id_orthogroup[.x],
                                        p2=reference_seq, s2=id_ref)) %>%
    mutate(taxon1=str_split_fixed(ID1,'\\.',2)[,1] ) %>%
    mutate(delta_ol = abs(100-OL1), delta_pid = abs(100-PID1) ) %>%
    group_by(taxon1) %>%
    arrange( delta_ol, delta_pid ) %>%
    dplyr::slice_head(n=1) %>%
    ungroup() %>%
    arrange( delta_ol, delta_pid )

    return(fastaseq[CLOSEST$ID1])
}

open_fasta = function(fastafile, fastaurl,debug=F){
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
    Sys.sleep(1)
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
                               debug=F, force=F){

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

    NODE_DESC =  og_data$node_desc
    CLADE_DESC = og_data$clade_desc
    CLADE_DIR = here::here('data','eggnog',NODE_DESC,CLADE_DESC)
    dir.create(CLADE_DIR,showWarnings = F,recursive = T)
    dir.create(file.path(CLADE_DIR,'fasta'),showWarnings = F,recursive = T)
    dir.create(file.path(CLADE_DIR,'tree'),showWarnings = F,recursive = T)

    og_fastafile = file.path(FASTA_DIR,paste0(OG,'.fasta'))
    fasta_in = open_fasta(og_fastafile, og_data$url_fasta)

    OG_unique = rep(OG, times=nref) |> makeunique::make_unique(sep='.',wrap_in_brackets = F)
    fasta_out = sprintf("%s_%s_%s_1to1-orthologs.fa",CLADE_FNAME,OG_unique,REF_ID)
    tree_out =  str_replace(fasta_out,'\\.fa$','.nwk') #%>% paste0("_eggnog")
    fastapath = file.path(CLADE_DIR,"fasta",fasta_out)
    treepath = file.path(CLADE_DIR,"tree",tree_out)
    has_output = file.exists(fastapath)

    if(nref==2){
      ogs = paste0(OG_unique,collapse=" ")
      if(debug){
        .dbg$log(sprintf('(i=%s) orthogroup %s has %s reference orthologs (to be splitted into: %s)',i,OG,nref,ogs))
      }
    }

    for(R in 1:nref){

      refid = REF_ID[R]
      refog = OG_unique[R]

      outputs=check_ortho(fastapath[R],treepath[R])
      fasta_1to1 = outputs$fasta
      tree_1to1  = outputs$tree

      if( !is.null(fasta_1to1) && !is.null(tree_1to1) ){
        if(debug){  .dbg$log(sprintf('(i=%s) orthogroup %s already processed',i,refog)) }
        break
      }else{
        if(is.null(fasta_1to1)){
          fasta_1to1 = reduce_orthologs(fasta_in, refid, CLADE$id)
          if(is_one2one){ .dbg$log(sprintf('(i=%s) orthogroup %s has 1-to-1 ortholog',i,refog)) }
          if(debug){
              .dbg$log(sprintf('(i=%s - %s) orthogroup %s has %s taxons for %s duplicated orthologs',i,CLADE_FNAME,OG,n_distinct(CLADE$taxid),sum(CLADE$is_dup)))
          }
        }

        if(is.null(tree_1to1)){
          # Final set of sequences with 1-to-1 orthologs for current clade
          DF_CLADE_1TO1 = CLADE %>% filter( id %in% names(fasta_1to1) ) %>% arrange(match(id, names(fasta_1to1)))

          ref_tree_orthogroup = inner_join(ref_tree %>% as_tibble(), DF_CLADE_1TO1,by=c('label'='taxid'))
          tree_1to1 = keep.tip(ref_tree,ref_tree_orthogroup$label)
          tree_1to1$tip.label = ref_tree_orthogroup$id
        }

        write_r4s_input(tree = tree_1to1,  fasta = fasta_1to1,
                        outfasta = fastapath[R], outtree=treepath[R])
      }
    }
  }
  #return(og_ref)
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

sptree=here::here(fu_dir,"basic_sptree","cog_100-alg_concat_default-raxml_default")
fungi_19 = treeio::read.tree(file.path(sptree,"fungi179-19orthogroups-seqs-sorted.fa.gz.final_tree.nw"))

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
#fungi_19$tip.label = fu_ncbi_eggnog$taxid[match(tolower(fungi_19$tip.label),fu_ncbi_eggnog$taxid)]

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
                          loading.call = { count_clade_orthologs(fu_og, fu_clades_toprocess) %>% bind_rows() },
                          doing = "orthogroups by fungi subclades...")

fu_ref = fu_orthogroups %>%
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
chunk_size=100
chunks=seq(1,nrow(fu_toprocess),by=chunk_size)
pbmcapply::pbmclapply(chunks,FUN = function(irow){
                     filter_orthogroups(
                          orthologs_count = fu_toprocess[irow:(irow+chunk_size),],
                          ref_tree = fungi_19,
                          ref_sp = 4932,
                          debug = F,
                          force = F) },mc.cores = 14)

# Select orthogroups for building the fungi species tree
# ---> orthogroups must be at the taxonomic level of fungi (id=4751)
# ---> all species must be present (179 fungal species)
# ---> at least 95% of orthologs should be 1-to-1
# ---> s. cerevisiae ortholog must be single (not duplicated)
fu_sptree = fu_ref %>%
            filter(clade_id == node_id & clade_f1to1>0.975 & clade_f == 1 & node_has_ref == 1)

#fu_processed = fu_ref %>% filter(node_has_ref %in% c(1,2) & clade_has_ref %in% c(0,1,2))


#### METAZOA ####
# Metazoa = 33208
# Mammalia = 40674

mz_dir = here::here('data','eggnog','33208_Metazoa')
mz_node    = find_eggnog_node(node = 33208,.print = F) %>% dplyr::rename_with(~Pxx(.,px='node',s="_"))
mzNOG      = get_eggnog_node(node = 33208)
mz_tax     = get_eggnog_taxonomy(33208)
mz_species = get_eggnog_species(node = 33208)

####_download NCBI common tree from list of species (manual) ####
# write_lines(mz_species$taxid, here::here('data','ncbi','33208-metazoa-161taxids.txt'))
# write_lines(mz_species$taxon, here::here('data','ncbi','33208-metazoa-161taxons.txt'))
metazoa = treeio::read.tree(here::here('data','ncbi','ncbi-metazoan.phy'))
metazoa$tip.label = str_remove_all(metazoa$tip.label,"['\\[\\]]") #%>% str_replace_all(" ","_")
metazoa_clades = c('7742'='vertebrata','40674'='mammalia','50557'='insecta','6231'='nematoda')
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

mz_clades_toprocess = mz_clades %>% filter(clade_size > 5 | clade_has_9606)

library(ggtree)
metazoa = treeio::read.tree(here::here('data','ncbi','ncbi-metazoan.phy'))
df_metazoa = metazoa %>% as_tibble() %>%
             mutate( label = stringr::str_remove_all(label,pattern="'"),
                     depth = ape::node.depth(metazoa),
                     is_leaf = depth == 1) %>%
             left_join(mz_ncbi_eggnog %>% dplyr::select(ncbi_name,ncbi_Name,eggnog_name,taxid), by=c('label'='ncbi_Name')) %>%
             mutate( to_process = label %in% mz_clades_toprocess$clade_name) %>%
             left_join(mz_clades, by=c('label'='clade_name'))

df_clade = df_metazoa %>% filter( !is_leaf )

ggtree(metazoa,layout = 'dendrogram',branch.length = 'none') %<+% df_metazoa +
  geom_nodepoint(aes(subset=node %in% df_clade$node, col = to_process),size=2) +
  geom_text_repel(aes(subset=node %in% df_clade$node, label=clade_desc, col = to_process),
                  fontface='bold', size=3)

mz_human   = eggnog_annotations_species(node = 33208, species = c(9606))

mz_og  = preload(saved.file = here::here("data/eggnog/33208_Metazoa-orthogroups.rds"),
                 loading.call = { count_eggnog_orthologs(33208) },
                 doing = "metazoa orthogroups...")

mz_orthogroups  = preload(saved.file = here::here("data/eggnog/33208_Metazoa-orthogroups-clades_over_5sp.rds"),
                          loading.call = { count_clade_orthologs(mz_og, mz_clades_toprocess) %>% bind_rows() },
                          doing = "orthogroups by metazoa subclades...")

mz_ref = mz_orthogroups %>%
  mutate(ref_sp = 9606,
         # remove orthogroup from fungi node present in other subnodes
         node_only = !( clade_id == 33208 & OG %in% OG[clade_id != 33208] ),
         node_has_ref = map(node_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist,
         clade_has_ref = map(clade_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist)

mz_toprocess = mz_ref %>%
  filter(node_has_ref %in% c(1,2) & (clade_ns > 20 | clade_has_ref %in% c(1,2)) & clade_f > 0.5 ) %>%
  group_by(node_id,clade_id) %>% add_count(name='clade_n')

dim(mz_orthogroups)
dim(mz_ref)
dim(mz_toprocess)
table(mz_toprocess$clade_name,mz_toprocess$node_has_ref)
janitor::tabyl(mz_toprocess,clade_name,clade_ns)

####_filter orthogroups fasta to keep at most 179 species (ortholog closest to yeast) ####
chunk_size=100
chunks=seq(1,nrow(mz_toprocess),by=chunk_size)
pbmcapply::pbmclapply(chunks,FUN = function(irow){
  filter_orthogroups(
    orthologs_count = mz_toprocess[irow:(irow+chunk_size),],
    ref_tree = metazoa,
    ref_sp = 9606,
    debug = F,
    force = F) },mc.cores = 6)


mz_ref %>% filter(node_has_ref == 2 & clade_has_ref == 2)


# Filter orthogroups fasta to keep at most 179 species (ortholog closest to yeast)
mz_processed = filter_orthogroups(orthologs_count = mz_ref,
                                  ref_tree = metazoa, ref_sp = 9606,
                                  debug = F, force = F)


#### MAMMALIA ####
# mammalia = treeio::read.tree(here::here('data','ncbi','ncbi-mammalia.phy'))
# mammalia$tip.label = str_remove_all(mammalia$tip.label,"['\\[\\]]") #%>% str_replace_all(" ","_")
# mammalia_clades = c('40674'='mammalia','314146'='euarchontoglires','9443'='primates',
#                     '9989'='rodentia','33554'='carnivora','91561'='cetartiodactyla',
#                     '9397'='chiroptera','311790'='Afrotheria')
#
#
# ma_dir = here::here('data','eggnog','40674_Mammalia')
# maNOG = get_eggnog_node(40674)
# ma_tax=get_eggnog_taxonomy(40674)
#
# ma_species = get_eggnog_species(node = 40674)
# write_lines(ma_species$taxid, here::here('data','ncbi','40674-mammalia-70taxids.txt'))
# write_lines(ma_species$taxon, here::here('data','ncbi','40674-mammalia-70taxons.txt'))
#
# ma_ncbi = preload( saved.file = here('data','ncbi','ncbi-to-eggnog-mammalia-70species.rds'),
#                    { match_strings(SP1=mammalia$tip.label, SP2=ma_species$taxon, use_soundex = F, manual = T) },
#                    'match ncbi species tree to eggnog mammalian species...')
#
# mammalia_eggnog = treeio::read.tree(here::here('data','ncbi','ncbi-mammalia.phy'))
# mammalia_eggnog$tip.label = ma_ncbi_eggnog$taxid
#
# laurasiatheria = extract.clade(mammalia_eggnog, 'Laurasiatheria', root.edge = 0, collapse.singles = TRUE)
# laurasiatheria_clades=c('cetartiodactyla'=91561,'carnivora'=33554,'chiroptera'=9397)
# id_laurasiatheria = laurasiatheria$tip.label %>% as.numeric %>% sum
#
#pbmcapply::pbmclapply(mzNOG$url_fasta, safe_download, mc.cores=14, path=mz_dir, ext='.fasta')
#mz_ali = get_eggnog_alignment(node=33208, use_trimmed = F)
#mz_clades = map(metazoa_clades, ~get_eggnog_taxonomy(.x))



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
