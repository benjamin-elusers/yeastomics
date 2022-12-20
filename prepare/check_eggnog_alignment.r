library(tidyverse)
count.fasta = function(fastafile){
  return(sum(grepl("^>",readLines(fastafile))))
}

count.file = function(directory, fileext=""){
  resfile = list.files(directory, pattern=fileext, full.names=T, recursive=T)
  if(length(fileext)>1){
    return( map(fileext, ~str_subset(resfile, pattern = paste0("\\.",.x,"$")) %>% length ) )
  }
  return(length(resfile))
}

count.trimal = function(directory){
  cat(sprintf("count trimal results %s...\n",basename(directory)))
  params = paste0("trimal_gap",c(50,80,90,'pyout'))
  map(params, ~count.file(file.path(directory,.x), paste0("\\.",.x,"$"))) %>%
    set_names(nm = paste0("n_",params)) %>% as_tibble()
}

count.r4s = function(directory){
  cat(sprintf("count rate4site results %s...\n",basename(directory)))
  params = paste0("trimal_gap",c(50,80,90,'pyout'))
  r4s_ext="\\.eggnog_sptree\\.r4s_raw"
  count.file(file.path(directory,"r4s"), paste0(params,r4s_ext)) %>%
    set_names(nm = paste0("n_r4s_",params)) %>% as_tibble()
}

wexac_fu="/media/WEXAC/EGGNOG/4751_Fungi"
clade_regex="[0-9]+_[a-zA-Z]+_[0-9]+sp"
col_clade = paste0("clade_",c("id","name","ns"))

clades = list.dirs(wexac_fu,recursive = F) %>% basename %>% str_subset(clade_regex)
df_clade = tibble(clade_desc=clades) %>%
           separate(col='clade_desc', into=col_clade, sep = "_", remove=F) %>%
           rowwise() %>%
           mutate( resdir = file.path(wexac_fu,clade_desc),
                   n_muscle = count.file(file.path(resdir,"muscle"), "\\.mu$"),
                   n_trimal = count.trimal( resdir ),
                   n_r4s = count.r4s( resdir )
           )

# with(df_count,
#     cat(sprintf("
#     [%s]
#     --> muscle = %5s
#     --> trimal = %5s (g50=%4s | g80=%4s | g90=%4s | gout=%4s)
#     --> r4s    = %5s (g50=%4s | g80=%4s | g90=%4s | gout=%4s)
#     ",
#     dir_clade,
#     n_muscle,
#     n_trimal,n_trimal_gap50,n_trimal_gap80,n_trimal_gap90,n_trimal_gappyout,
#     n_r4s,n_r4s_gap50,n_r4s_gap80,n_r4s_gap90,n_r4s_gappyout))
#   )
# }

muscle_files = list.files(wexac_fu, recursive = T,  pattern="\\.mu$", full.names=T)
trimal_files = list.files(wexac_fu, recursive = T,  pattern="\\.trimal_gap(50|80|90|pyout)$", full.names=T)
#old_trimal_files = list.files(wexac_fu, recursive = T,  pattern="\\.mu.trimal_gap(10|20)$", full.names=T)
#unlink(old_trimal_files)

table(str_extract(dirname(muscle_files),"[0-9]+_[a-zA-Z]+_[0-9]+sp"))
table(str_extract(dirname(trim_files),"[0-9]+_[a-zA-Z]+_[0-9]+sp"))

r4s_files  = list.files(wexac_fu, recursive = T,  pattern="\\.r4s_raw_eggnog$", full.names=T)




CLADE_DESC = str_extract(dirname(r4s_files),"[0-9]+_[a-zA-Z]+_[0-9]+sp")
CLADE = tibble( clade_desc = CLADE_DESC %>% unique,
                clade_id   = str_split_fixed(clade_desc,"_",3)[,1],
                clade_name = str_split_fixed(clade_desc,"_",3)[,2],
                clade_ns   = str_split_fixed(clade_desc,"_",3)[,3])


clade_desc = basename(dirname(muscle_files))
orthogroup = basename(muscle_files) %>% str_split_fixed(string = ., "_",n=4)
orthoname = apply(orthogroup[,1:2],1,paste0,collapse="_")

muscle_gaps = pbmcapply::pbmcmapply(muscle_files, FUN=count_gaps, mc.cores = 14)
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

trim_gaps = pbmcapply::pbmcmapply(trimmed_files, FUN=count_gaps, mc.cores = 14)
names(trim_gaps) = orthoname
trim_gapstat = tibble(ortho=orthoname,
                        clade=clade_desc,
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
  group_by(clade) %>% add_count(name='nclade')

library(tidyverse)
library(ggplot2)

p=ggplot(muscle_gapstat) +
  geom_density(aes(x=agap)) +
  geom_text(aes(label=paste0('n=',nclade)),x=Inf,y=Inf,hjust='inward',vjust='inward',check_overlap = T,size=3) +
  facet_wrap(~clade,scales='free_y',nrow = 2,) +
  xlab('Average % gap in alignment') +
  geom_density(data=trim_gapstat,aes(x=agap),col='red')+
  theme_clean()
p

ggsave(path=here::here("plots"),
       filename="4751_Fungi-gaps_frequency_eggnog.pdf",
       plot=p,
       scale=1.5)

###
source(here::here("src","__setup_yeastomics__.r"))
fudir = here::here("data","eggnog","4751_Fungi")
r4s_res=pbmcapply::pbmcmapply(r4s_files,FUN = read.R4S, mc.cores = 14, SIMPLIFY=F) # DO NOT SIMPLIFY TO KEEP AS A LIST
saveRDS(r4s_res,file.path(fudir,"rate4site-fungi-clades.rds"))

df_r4s = bind_rows(r4s_res) %>%
         mutate(clade = str_extract(dirname(r4sfile),"[0-9]+_[a-zA-Z]+_[0-9]+sp"),
                r4sfile  = names(r4s_res) %>% basename,
                OG = str_split_fixed(r4sfile,'_',n=3)[,2] )

df_rate = df_r4s %>% group_by(clade,OG) %>%
          mutate(ER=mean_(SCORE)) %>%
          dplyr::select(OG,clade,ER) %>%
          distinct() %>%
          pivot_wider(id_cols=OG, names_from=clade, values_from = ER, names_prefix = 'ER.')

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
       plot=p1,
       scale=1.8)
ggsave(path=here::here("plots"),
       filename="4751_Fungi-evorate-clades-scatterplot.pdf",
       plot=p1.3,
       scale=1.8)

