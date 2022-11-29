library(tidyverse)
library(here)
source(here::here("src","__setup_yeastomics__.r"))
library(treeio)
library(tidytree)
ncbi_dir = here::here("data","ncbi")
#eggnog_fasta = Biostrings::readAAStringSet("http://eggnog5.embl.de/download/latest/e5.proteomes.faa")

filter_orthogroups = function( orthologs_count, ref_sp = 4932,ref_tree,
                               debug=F, force=F, min_species=5 ){

  if(missing(orthologs_count)){
    .warn$log('Must input a taxonomic level and some valid clades...')
    node = find_eggnog_node()
    clades = get_eggnog_taxonomy(node$id,.print = T) %>% mutate(clade_desc = sprintf("%s_%s (n=%s)",clade_id,clade_name,clade_size))
    selected = select.list(clades$clade_desc,multiple = T)
    orthologs_count = count_taxons_eggnog_node(node$id, clades$clade_id[ clades$clade_desc %in% selected] ) %>% bind_rows()
  }

  .info$log("remove orthogroups without reference species...")
  df_og = orthologs_count %>%
    mutate( node_has_ref = map(node_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist,
            clade_has_ref = map(clade_orthologs, ~ sum(.x$taxid %in% ref_sp)) %>% unlist)

  og_count =  df_og %>% ungroup() %>% mutate(nog = n_distinct(OG)) %>%
    group_by(node_id,node_name,node_has_ref) %>% mutate( node_nog =n_distinct(OG) ) %>%
    group_by(node_id,node_name,clade_id,clade_name,clade_has_ref) %>% add_count(name='clade_nog') %>%
    dplyr::select(node_id,node_name,node_size, nog, node_has_ref, node_nog,
                  clade_has_ref,clade_id,clade_name,clade_has_ref,clade_nog) %>%
    distinct()
  og_ref = df_og %>% filter(node_has_ref > 0  & clade_size >  min_species & clade_ns >  min_species ) %>%
    group_by(node_id,clade_id) %>% add_count(name='clade_n')

  #df_og %>% janitor::tabyl(node_has_ref,clade_name)
  #df_og %>% janitor::tabyl(clade_has_ref,clade_name)
  N0 = nrow(df_og)
  n0 = n_distinct(df_og$OG)
  has_ref= og_ref$node_has_ref>0
  n1 = n_distinct(og_ref$OG[has_ref])
  N1 = sum(has_ref)
  .succ$log(sprintf("number of orthogroups with reference species  = %s/%s (total=%s/%s)",n1,n0,N1,N0))
  has_many_ref= og_ref$node_has_ref>1
  n2.0 = n_distinct(og_ref$OG[has_ref & has_many_ref])
  N2.0 = sum(has_many_ref)
  .succ$log(sprintf("   ---> single reference species  = %s/%s (total=%s/%s)",n2.0,n1,N2.0,N1))
  has_single_ref= og_ref$node_has_ref==1
  n2 = n_distinct(og_ref$OG[has_single_ref])
  N2 = sum(has_single_ref)
  .succ$log(sprintf("   ---> single reference species  = %s/%s (total=%s/%s)",n2,n1,N2,N1))
  is_1to1 = og_ref$clade_one2one==T
  n3 = n_distinct(og_ref$OG[has_single_ref & is_1to1])
  N3 = sum(has_single_ref & is_1to1)
  .succ$log(sprintf("   -------> one-to-one orthologs  = %s/%s (total=%s/%s)",n3,n2,N3,N2))
  n3.0 = n_distinct(og_ref$OG[has_single_ref & !is_1to1])
  N3.0 = sum(has_single_ref & !is_1to1)
  .succ$log(sprintf("   -------> one-to-many orthologs = %s/%s (total=%s/%s)",n3.0,n2,N3.0,N2))

  for( i in 1:nrow(og_ref) ){

    #cat(sprintf("%-5s/r",i))
    OG = og_ref$OG[i]
    NPROT = og_ref$clade_np[i]
    is_one2one = og_ref$clade_one2one[i]

    NODE_ID  = og_ref$node_id[i]
    NODE_NAME = og_ref$node_name[i]
    NODE  = og_ref$node_orthologs[[i]]
    FASTA_DIR = here::here('data','eggnog','fasta',paste0(NODE_ID,'-',NODE_NAME))

    REF_ID = NODE$id[ NODE$taxid == ref_sp ]
    nref = length(REF_ID)

    CLADE_NAME = og_ref$clade_name[i]
    CLADE_ID = og_ref$clade_id[i]
    CLADE_N = og_ref$clade_n[i]
    CLADE = og_ref$clade_orthologs[[i]] %>%
      mutate(ref = (taxid==ref_sp), is_dup = is.dup(taxid) ) %>%
      arrange(desc(ref),desc(is_dup))
    ntaxdup = n_distinct(CLADE$taxid[CLADE$is_dup])
    ndup = sum(CLADE$is_dup)
    ID_DUP = CLADE$id[CLADE$is_dup]
    RANGE_DUP = seq_along(ID_DUP)

    NODE_DESC = paste0(NODE_NAME,"_",NODE_ID)
    CLADE_DESC = paste0(CLADE_NAME,"_",CLADE_ID)
    CLADE_DIR = here::here('data','eggnog',NODE_DESC,CLADE_DESC)
    dir.create(CLADE_DIR,showWarnings = F,recursive = T)

    CLADE_NS = og_ref$clade_ns[i]

    if(nref==0){
      .error$log(sprintf('(i=%s) orthogroup %s has no reference sequence (taxid=%s)',i,OG,ref_sp))
    }else if(nref==1){

      fasta_filename = sprintf("%s-%s_sp-%s-%s-1to1_orthologs.fa",CLADE_NAME,CLADE_NS,OG,REF_ID)
      tree_filename =  str_replace(fasta_filename,'\\.fa$','.nwk')
      has_file = file.exists( file.path(CLADE_DIR,fasta_filename) )

      if( has_file ){
        .dbg$log(sprintf('(i=%s) orthogroup %s already processed',i,OG))
      }else{
        if(debug){
          .dbg$log(sprintf('get fasta sequences for the full orthogroup (%s proteins)...',NPROT))
        }

        og_fasta = file.path(FASTA_DIR,paste0(OG,'.fasta'))
        if( !file.exists(og_fasta) ){
          # Get fasta sequence of the full orthogroup
          temp = tempfile()
          download.file(og_ref$url_fasta[i],destfile = temp,quiet = T)
          fasta = Biostrings::readAAStringSet(temp)
        }else{
          fasta = Biostrings::readAAStringSet(og_fasta)
        }

        # Get fasta sequence of the single reference sequence
        ref_seq = fasta[[REF_ID]] %>%  as.character() %>% chartr("U","X",x = .) %>% setNames(REF_ID)

        if( is_one2one ){
          .succ$log(sprintf('(i=%s) orthogroup %s has 1-to-1 ortholog',i,OG))
          fasta_1to1 = fasta[CLADE$id]

          ref_tree_orthogroup = inner_join(ref_tree %>% as_tibble(), CLADE,by=c('label'='taxid'))
          tree_1to1 = keep.tip(ref_tree,ref_tree_orthogroup$label)

        }else{

          # Count and get fasta sequences of duplicated orthologs from the current clade
          if(debug){
            .dbg$log(sprintf('(i=%s - %s) orthogroup %s has %s taxons for %s duplicated orthologs',i,CLADE_NAME,OG,ntaxdup,ndup))
          }
          ortho2many_seq = fasta[ID_DUP] %>% as.character() %>% chartr("U","X",x=.)

          # Align each duplicated ortholog to the reference sequence and select best representative per taxon
          #.dbg$log('select most representative orthologs compared to reference sequence...')
          CLOSEST = map_dfr(RANGE_DUP, ~score_ali(p1=ortho2many_seq[.x], s1=ID_DUP[.x], p2=ref_seq, s2=REF_ID)) %>%
            mutate(taxon1=str_split_fixed(ID1,'\\.',2)[,1] ) %>%
            mutate(delta_ol = abs(100-OL1), delta_pid = abs(100-PID1) ) %>%
            group_by(taxon1) %>%
            arrange( delta_ol, delta_pid ) %>%
            #filter( nearest(100,OL1) & nearest(PID1,100)) %>%
            dplyr::slice_head(n=1)

          ONE2ONE =  CLADE %>% filter( !is.dup(taxid) )
          DUP_REPRE = CLADE %>% filter( id %in% CLOSEST$ID1 )
          DF_CLADE_1TO1 = bind_rows(ONE2ONE,DUP_REPRE)

          # Final set of sequences with 1-to-1 orthologs for current clade
          fasta_1to1 = fasta[DF_CLADE_1TO1$id]
          NP_1to1 = n_distinct(DF_CLADE_1TO1$taxid)
          NS_1to1 = n_distinct(DF_CLADE_1TO1$id)
          if(debug){
            .dbg$log(sprintf("ns = %s np = %s",NP_1to1,NS_1to1))
          }

          ref_tree_orthogroup = inner_join(ref_tree %>% as_tibble(), DF_CLADE_1TO1,by=c('label'='taxid'))
          tree_1to1 = keep.tip(ref_tree,ref_tree_orthogroup$label)


          fasta_path = file.path(CLADE_DIR,fasta_filename)
          writeXStringSet(fasta_1to1,fasta_path,append = F,format = 'fasta')

          tree_path = file.path(CLADE_DIR,tree_filename)
          write.tree(tree_1to1,tree_path)
        }
      }
    }else if(nref> 1){
      .warn$log(sprintf('(i=%s) orthogroup %s has multiple reference sequences (taxid=%s n=%s)',i,OG,ref_sp,nref))
      #print(og_ref[i,])
      #readline('continue?')
    }
  }
  return(og_ref)
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

fu_dir = here::here("data","eggnog","Fungi_4751")
fuNOG      = get_eggnog_node(node = 4751)
fu_tax     = get_eggnog_taxonomy(4751)
fu_species = get_eggnog_species(node = 4751)

####_NCBI tree from taxons/taxids (manual) ####
# write_lines(fu_species$taxid, here::here('data','ncbi','4751-fungi-179taxids.txt'))
# write_lines(fu_species$taxon, here::here('data','ncbi','4751-fungi-179taxons.txt'))
# fungi=treeio::read.tree("/home/benjamin/Desktop/GitHub/yeastomics/data/ncbi/ncbi-fungi.phy")
# fungi$tip.label = str_remove_all(fungi$tip.label,"['\\[\\]]") #%>% str_replace_all(" ","_")
# fu_ncbi = preload( saved.file = file.path(ncbi_dir,'ncbi-to-eggnog-fungi-179species.rds'),
#                    { match_strings(fungi$tip.label, SP2=fu_species$taxon, use_soundex = F, manual = T) },
#                    'match ncbi species tree to eggnog fungal species...')
# fu_ncbi_eggnog = fu_ncbi %>%
#   left_join( fu_species %>% mutate(Taxon=taxon, taxon=tolower(taxon)), c('s2'='taxon')) %>%
#   dplyr::rename(ncbi_name=s1,eggnog_name=s2) %>%
#   dplyr::select(-c(is_identical:is_substring,osa:n1))

fu_ncbi = get_ncbi_tree(fu_species$taxid)
#fungi_eggnog$tip.label = fu_ncbi_eggnog$taxid


library(taxize)
fu_ncbi  = subtrees(fungi_eggnog) %>%
           set_names(fungi_eggnog$node.label %>% str_remove_all("['\\[\\]]"))  %>%
           map_df(~.x$tip.label %>% n_distinct %>% as_tibble, .id = "clades")  %>%
           filter(value>4 & clades != "") %>%
           mutate( ncbi_id = taxize::get_ids(clades, db='ncbi',verbose = F)$ncbi )

fu_yeast   = eggnog_annotations_species(node = 4751, species = c(4932,4896))
fu_og  = count_eggnog_orthologs(4751)

fu_taxons  = count_clade_orthologs(4751, subnode=c(4751,4890,5204,451866,4891,147541,147545,147550)  )

#fu_og_
#subnode=c(4751,4890,5204,451866,4891,147541,147545,147550)


fu_taxons  = count_taxons_eggnog_node(4751, subnode=c(4751,4890,5204,451866,4891,147541,147545,147550)  )

# remove orthogroup from fungi node present in other subnodes
fu_orthogroups = bind_rows(fu_taxons) %>%
  mutate( node_only = !( clade_id == 4751 & OG %in% OG[clade_id != 4751] ),
          has_4932 = map(clade_orthologs,~ sum(.x$taxid == 4932)))
#og_4932 = fu_yeast %>% filter(taxid == 4932) %>% pull(OG,string)

# Filter orthogroups fasta to keep at most 179 species (ortholog closest to yeast)
fu_og = filter_orthogroups(orthologs_count = fu_orthogroups,
                           ref_tree = fungi_eggnog, ref_sp = 4932,
                           debug = F, force = F, min_sp = 5)

fu_ = fu_orthogroups %>% filter(clade_id == node_id & clade_ns == node_size & has_4932==1)
n_distinct(fu_179sp_with4932)

FASTA1TO1 = list.files(fu_dir,  recursive = T, pattern = '.fa', full.names = T,include.dirs = F,no.. = F)
ortho1to1 = str_split_fixed(basename(FASTA1TO1),"-",n=4) %>%
  as_tibble() %>%
  set_names(c('clade','nspecies','OG','refid')) %>%
  mutate(nspecies = str_remove(nspecies,'_sp$'),
         refid = str_remove(refid,"-1to1_orthologs.fa"),
         fastafile=FASTA1TO1) %>%
  separate(refid, into=c('taxid','string'), sep='\\.' ) %>%
  type_convert()

fungi_orthoseq = read.sequences(ortho1to1$fastafile,ncores = 14,type='AA')
ortho1to1$NSP = sapply(fungi_orthoseq,length)
sum(ortho1to1$nspecies != ortho1to1$NSP)
#unlink(ortho1to1$fastafile[ortho1to1$nspecies != NSP])

fungi_179sp = ortho1to1 %>%
  left_join(fu_orthogroups, by=c('clade'='clade_name','OG')) %>%
  filter(clade_id == 4751 )


#### METAZOA ####
# Metazoa = 33208
# Mammalia = 40674
metazoa = treeio::read.tree(here::here('data','ncbi','ncbi-metazoan.phy'))
metazoa$tip.label = str_remove_all(metazoa$tip.label,"['\\[\\]]") #%>% str_replace_all(" ","_")
metazoa_clades = c('7742'='vertebrata','40674'='mammalia','50557'='insecta','6231'='nematoda')

mz_species = get_eggnog_species(node = 33208)
write_lines(mz_species$taxid, here::here('data','ncbi','33208-metazoa-161taxids.txt'))
write_lines(mz_species$taxon, here::here('data','ncbi','33208-metazoa-161taxons.txt'))


#### MAMMALIA ####
mammalia = treeio::read.tree(here::here('data','ncbi','ncbi-mammalia.phy'))
mammalia$tip.label = str_remove_all(mammalia$tip.label,"['\\[\\]]") #%>% str_replace_all(" ","_")
mammalia_clades = c('40674'='mammalia','314146'='euarchontoglires','9443'='primates',
                    '9989'='rodentia','33554'='carnivora','91561'='cetartiodactyla',
                    '9397'='chiroptera','311790'='Afrotheria')


ma_dir = here::here('data','eggnog','40674_Mammalia')
maNOG = get_eggnog_node(40674)
ma_tax=get_eggnog_taxonomy(40674)

ma_species = get_eggnog_species(node = 40674)
write_lines(ma_species$taxid, here::here('data','ncbi','40674-mammalia-70taxids.txt'))
write_lines(ma_species$taxon, here::here('data','ncbi','40674-mammalia-70taxons.txt'))

ma_ncbi = preload( saved.file = here('data','ncbi','ncbi-to-eggnog-mammalia-70species.rds'),
                   { match_strings(SP1=mammalia$tip.label, SP2=ma_species$taxon, use_soundex = F, manual = T) },
                   'match ncbi species tree to eggnog mammalian species...')

ma_ncbi_eggnog = ma_ncbi %>%
  left_join( ma_species %>% mutate(Taxon=taxon, taxon=tolower(taxon)), c('s2'='taxon')) %>%
  dplyr::rename(ncbi_name=s1,eggnog_name=s2) %>%
  dplyr::select(-c(is_identical:is_substring,osa:n1))


mammalia_eggnog = treeio::read.tree(here::here('data','ncbi','ncbi-mammalia.phy'))
mammalia_eggnog$tip.label = ma_ncbi_eggnog$taxid

ma_human = eggnog_annotations_species(node = 40674, species = c(9606))
ma_taxons = count_taxons_eggnog_node(node=40674, subnode=c(40674,314146,9443,9989,33554,314294,9397,311790))

ma_taxons$`40674_Mammalia (n=70)`

laurasiatheria = extract.clade(mammalia_eggnog, 'Laurasiatheria', root.edge = 0, collapse.singles = TRUE)
laurasiatheria_clades=c('cetartiodactyla'=91561,'carnivora'=33554,'chiroptera'=9397)
id_laurasiatheria = laurasiatheria$tip.label %>% as.numeric %>% sum

ma


ma_orthogroups = bind_rows(ma_taxons) %>%
  mutate( node_only = !( clade_id == 40674 & OG %in% OG[clade_id != 40674] ),
          has_9606 = map(clade_orthologs,~ sum(.x$taxid == 9060)))




ma_og = filter_orthogroups(orthologs_count = ma_orthogroups,

                           ref_tree=mammalia_eggnog, ref_sp = 9606,
                           debug = T, force = T, min_sp = 1)

FASTA1TO1 = list.files(ma_dir,  recursive = T, pattern = '.fa', full.names = T,include.dirs = F,no.. = F)
ortho1to1 = str_split_fixed(basename(FASTA1TO1),"-",n=4) %>%
  as_tibble() %>%
  set_names(c('clade','nspecies','OG','refid')) %>%
  mutate(nspecies = str_remove(nspecies,'_sp$'),
         refid = str_remove(refid,"-1to1_orthologs.fa"),
         fastafile=FASTA1TO1) %>%
  separate(refid, into=c('taxid','string'), sep='\\.' ) %>%
  type_convert()

fungi_orthoseq = read.sequences(ortho1to1$fastafile,ncores = 14,type='AA')
ortho1to1$NSP = sapply(fungi_orthoseq,length)
sum(ortho1to1$nspecies != ortho1to1$NSP)
#unlink(ortho1to1$fastafile[ortho1to1$nspecies != NSP])


#pbmcapply::pbmclapply(mzNOG$url_fasta, safe_download, mc.cores=14, path=mz_dir, ext='.fasta')
#mz_ali = get_eggnog_alignment(node=33208, use_trimmed = F)
#mz_clades = map(metazoa_clades, ~get_eggnog_taxonomy(.x))


#### EGGNOG V6 ####
#"http://eggnog6.embl.de/download/eggnog_6.0/e6.og2level.tsv"
#"http://eggnog6.embl.de/download/eggnog_6.0/e6.og2seqs_and_species.tsv"
e6_dir = here::here('data','eggnog','eggnog_v6')
ncbi_dir = here::here('data','ncbi')
e6.tax = readr::read_delim(file.path(e6_dir,"e6.og2level.tsv"),
                           delim = '\t',col_names = c("node","OG"))

e6.species = readr::read_delim(file.path(e6_dir,"e6.og2seqs_and_species.tsv"),
                               delim = '\t',col_types = 'cciicc',
                               col_names = c("node","OG","nsp","nprot","taxids","string_ids"))

e6.ogs = readr::read_delim(file.path(e6_dir,"e6.seq2ogs.tsv"), delim = '\t',
                           col_types = 'cc',
                           col_names = c("string_id","OGs"))

e6.hierarchy = readr::read_delim(file.path(e6_dir,"e6.og2parents_and_children.new.tsv"),
                                 delim = '\t', col_types = 'ciicc',
                                 col_names = c("OG","nparent","nchild","OG_parent","OG_child"))


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