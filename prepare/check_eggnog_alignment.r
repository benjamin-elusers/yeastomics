library(tidyverse)
count.fasta = function(fastafile){
  return(sum(grepl("^>",readLines(fastafile))))
}

count.file = function(directory, fileext=""){
  if(length(fileext)>1){
    resfile = list.files(directory, full.names=T, recursive=T)
    return( map(fileext, ~str_subset(resfile, pattern = paste0("\\.",.x,"$")) %>% length ) )
  }else{
    resfile = list.files(directory, pattern=fileext, full.names=T, recursive=T)
    return(length(resfile))
  }
}

count.trimal = function(directory){
  cat(sprintf("count trimal results %s...\n",basename(directory)))
  params = paste0("trimal_gap",c(10,20,50,80,'pyout'))
  map(params, ~count.file(file.path(directory,.x), paste0("\\.",.x,"$"))) %>%
    set_names(nm = paste0("n_",params)) %>% as_tibble()
}

count.r4s = function(directory){
  cat(sprintf("count rate4site results %s...\n",basename(directory)))
  params = c('notrim', paste0("trimal_gap",c(10,20,50,80,"pyout")))
  r4s_dir = c('r4s_muscle', str_replace(params[-1],"trimal","r4s"))
  r4s_ext="\\.eggnog_sptree\\.r4s_raw"
  map(seq(params), ~count.file(file.path(directory,r4s_dir[.x]), paste0("\\.",params[.x],r4s_ext,"$"))) %>%
    set_names(nm = paste0("n_",params)) %>% as_tibble()
  #count.file(file.path(directory,"r4s"), paste0(params,r4s_ext)) %>%
  #  set_names(nm = paste0("n_r4s_",params)) %>% as_tibble()
}

wexac_fu="/media/WEXAC/EGGNOG/4751_Fungi"
clade_regex="[0-9]+_[a-zA-Z]+_[0-9]+sp"
col_clade = paste0("clade_",c("id","name","ns"))

clades = list.dirs(wexac_fu,recursive = F) %>% basename %>% str_subset(clade_regex)
df_clade = tibble(clade_desc=clades) %>%
           separate(col='clade_desc', into=col_clade, sep = "_", remove=F) %>%
           rowwise() %>%
           mutate( resdir = file.path(wexac_fu,clade_desc) )
df_clade = df_clade %>% mutate( n_muscle = count.file(file.path(resdir,"muscle"), "\\.mu$") )
df_clade = df_clade %>% mutate( n_trimal = count.trimal( resdir ) )
df_clade = df_clade %>% mutate( n_r4s = count.r4s( resdir ) )


df_clade$n_muscle - df_clade$n_trimal
df_clade$n_muscle - df_clade$n_r4s

muscle_files = list.files(wexac_fu, recursive = T,  pattern="\\.mu$", full.names=T)
#table(str_extract(dirname(muscle_files),clade_regex))
trimal_files = list.files(wexac_fu, recursive = T,  pattern="\\.trimal_gap(10|20|50|80|pyout)$", full.names=T)
#table(str_extract(dirname(trimal_files),clade_regex))
r4s_files  = list.files(wexac_fu, recursive = T,  pattern="\\.eggnog_sptree\\.r4s_raw$", full.names=T)
#table(str_extract(dirname(r4s_files),clade_regex))
#rowSums(df_clade$n_r4s) %>% setNames(df_clade$clade_desc)

clade_desc = dirname(muscle_files) %>% str_extract(clade_regex)
orthogroup = basename(muscle_files) %>% str_split_fixed(string = ., "_",n=4)
orthoname = apply(orthogroup[,1:2],1,paste0,collapse="_")
source(here::here("src","function_alignment.r"))
muscle_gaps = pbmcapply::pbmcmapply(count_gaps_fromfile, muscle_files,  mc.cores = 14)
names(muscle_gaps) = orthoname
muscle_gapstat = tibble(ortho=orthoname,
                        clade=clade_desc,
                        gaps=muscle_gaps) %>%
               rowwise() %>%
               mutate( agap = mean(gaps,na.rm=T),
                       sgap = sd(gaps,na.rm=T),
                       mgap = median(gaps,na.rm=T),
                       gap5 = quantile(gaps,probs=0.05),
                       gap20 = quantile(gaps,probs=0.2),
                       gap80 = quantile(gaps,probs=0.8),
                       gap95 = quantile(gaps,probs=0.95)) %>%
              mutate(og=str_split_fixed(ortho,"_",2)[,2]) %>%
              group_by(clade) %>% add_count(name='nclade')

clade_desc = dirname(trimal_files) %>% str_extract(clade_regex)
orthogroup = basename(trimal_files) %>% str_split_fixed(string = ., "_",n=4)
orthoname = apply(orthogroup[,1:2],1,paste0,collapse="_")
trimal_param = dirname(trimal_files) %>% basename()

trim_gaps = pbmcapply::pbmcmapply(trimal_files, FUN=count_gaps_fromfile, mc.cores = 14)
names(trim_gaps) = orthoname
trim_gapstat = tibble(ortho=orthoname,
                      clade=clade_desc,
                      params=trimal_param,
                      gaps=trim_gaps) %>%
  rowwise() %>%
  mutate( agap = mean(gaps,na.rm=T),
          sgap = sd(gaps,na.rm=T),
          mgap = median(gaps,na.rm=T),
          gap5 = quantile(gaps,probs=0.05),
          gap20 = quantile(gaps,probs=0.2),
          gap80 = quantile(gaps,probs=0.8),
          gap95 = quantile(gaps,probs=0.95)) %>%
  mutate(og=str_split_fixed(ortho,"_",2)[,2]) %>%
  group_by(clade,params) %>% add_count(name='nclade')

library(tidyverse)
library(ggplot2)
library(ggthemes)
pg=ggplot(muscle_gapstat) +
  geom_density(aes(x=agap)) +
  geom_text(aes(label=paste0('n=',nclade)),x=Inf,y=Inf,hjust='inward',vjust='inward',check_overlap = T,size=3) +
  geom_density(data=trim_gapstat,aes(x=agap,color=clade_desc)) +
  facet_wrap(~params,scales='free_y',nrow = 4) +
  xlab('Average % gap in alignment') +
  theme_clean()
pg

ggsave(path=here::here("plots"),
       filename="4751_Fungi-gaps_frequency_eggnog.pdf",
       plot=pg,
       scale=1.5)

###
source(here::here("src","__setup_yeastomics__.r"))
fudir = here::here("data","eggnog","4751_Fungi")
r4s_res=pbmcapply::pbmcmapply(r4s_files,FUN = read.R4S, mc.cores = 14, SIMPLIFY=F) # DO NOT SIMPLIFY TO KEEP AS A LIST
saveRDS(r4s_res,file.path(fudir,"rate4site-4751_fungi-clades.rds"))



df_r4s = tibble(r4sfile  = names(r4s_res),
                clade = str_extract(dirname(r4sfile),clade_regex),
                OG = str_split_fixed(basename(r4sfile),'_',n=3)[,2],
                r4s_res = r4s_res,
                param = str_extract(basename(r4sfile),"\\.[^\\.]*trim[^\\.]*\\.") ) %>%
         unnest(r4s_res)

df_rate = df_r4s %>% group_by(clade,OG) %>%
          mutate(ER=mean_(SCORE)) %>%
          dplyr::select(OG,clade,param,ER) %>%
          distinct() %>%
          pivot_wider(id_cols=OG, names_from=c(clade,param), names_glue = "{clade}{param}", values_from = ER, names_prefix = 'ER.')

COR_ER_clade=cor(df_rate[,-1],use='pairwise.complete',met='spearman')

fungi=ape::read.tree(file=here::here("data/ncbi/ncbi-fungi.phy"))
fu_ncbi = preload( saved.file = file.path(here::here("data/ncbi/"),'ncbi-to-eggnog-fungi-179species.rds'),
                   { match_strings(fungi$tip.label, SP2=fu_species$taxon, use_soundex = F, manual = T) },
                   'match ncbi species tree to eggnog fungal species...')
fu_ncbi_eggnog = fu_ncbi %>%
  left_join( get_eggnog_species(4751) %>% mutate(Taxon=taxon, taxon=tolower(taxon)), c('s2'='taxon')) %>%
  dplyr::rename(ncbi_name=s1,eggnog_name=s2) %>%
  dplyr::select(-c(is_identical:is_substring,osa:n1)) %>%
  mutate(ncbi_Name = str_to_sentence(ncbi_name))
library(tidytree)
library(ggtree)
df_fungi = fungi %>% as_tibble() %>%
           mutate( label = stringr::str_remove_all(label,pattern="'"),
                   depth = ape::node.depth(fungi),
                   is_leaf = depth == 1) %>%
           left_join(fu_ncbi_eggnog, by=c('label'='ncbi_Name')) %>%
           left_join(CLADE, by=c('label'='clade_name'))
df_yeast = df_fungi %>% filter(is_leaf & label %in% c('Saccharomyces cerevisiae','Schizosaccharomyces pombe'))

df_clade = df_fungi %>% filter(label %in% clade_name)

p1.1=ggtree(fungi,ladderize = T,branch.length = 'none') %<+% df_fungi +
  geom_text2(aes(subset=node %in% df_clade$node, label=clade_desc),size=4,col='blue',hjust=1,vjust=-1) +
  geom_nodepoint(aes(subset=node %in% df_clade$node, label=clade_desc),size=3,col='blue') +
  geom_tiplab(aes(subset=node %in% df_yeast$node),size=3) +
  geom_tippoint(aes(subset=node %in% df_yeast$node), col=c('dodgerblue','orange')) +
  xlim(-5,16)
p1.1

p1.2=ggcorrplot::ggcorrplot(COR_ER_clade,
                            lab=T,
                            method = 'circle',type = 'lower')

library(patchwork)

p1 = p1.1 | p1.2

p1.3=GGally::ggpairs(df_rate[,-1])

ggsave(path=here::here("plots"),
       filename="4751_Fungi-evorate-clades.pdf",
       plot=p1.2,
       scale=3)
ggsave(path=here::here("plots"),
       filename="4751_Fungi-evorate-clades-scatterplot.pdf",
       plot=p1.3,
       scale=1.8)

