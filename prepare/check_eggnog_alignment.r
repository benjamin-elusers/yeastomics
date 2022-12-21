library(tidyverse)
library(tidytree)
library(ggtree)
library(ggplot2)
library(ggthemes)
library(patchwork)
source(here::here("src","__setup_yeastomics__.r"))
fudir = here::here("data","eggnog","4751_Fungi")

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

count.muscle = function(directory){
  cat(sprintf("count muscle alignments %s...\n",basename(directory)))
  return( count.file(file.path(directory,"muscle"),"\\.mu$") )
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
}

wexac_fu="/media/WEXAC/EGGNOG/4751_Fungi"
clade_regex="[0-9]+_[a-zA-Z]+_[0-9]+sp"
col_clade = paste0("clade_",c("id","name","ns"))

clades = list.dirs(wexac_fu,recursive = F) %>% basename %>% str_subset(clade_regex)
df_clade = tibble(clade_desc=clades) %>%
           separate(col='clade_desc', into=col_clade, sep = "_", remove=F) %>%
           rowwise() %>%
           mutate( resdir = file.path(wexac_fu,clade_desc) )

df_clade = df_clade %>% mutate( n_muscle = count.muscle(resdir ) )
df_clade = df_clade %>% mutate( n_trimal = count.trimal( resdir ) )
df_clade = df_clade %>% mutate( n_r4s = count.r4s( resdir ) )


df_clade$n_muscle - df_clade$n_trimal
df_clade$n_muscle - df_clade$n_r4s

muscle_files = list.files(wexac_fu, recursive = T,  pattern="\\.mu$", full.names=T)
#table(str_extract(dirname(muscle_files),clade_regex))
#df_clade$n_muscle %>% setNames(df_clade$clade_desc)
trimal_files = list.files(wexac_fu, recursive = T,  pattern="\\.trimal_gap(10|20|50|80|pyout)$", full.names=T)
#table(str_extract(dirname(trimal_files),clade_regex))
#rowSums(df_clade$n_trimal) %>% setNames(df_clade$clade_desc)
r4s_files  = list.files(wexac_fu, recursive = T,  pattern="\\.eggnog_sptree\\.r4s_raw$", full.names=T)
#table(str_extract(dirname(r4s_files),clade_regex))
#rowSums(df_clade$n_r4s) %>% setNames(df_clade$clade_desc)

save.image(here::here("data/eggnog","4751_Fungi-eggnog-results.rdata"))


clade_desc = dirname(muscle_files) %>% str_extract(clade_regex)
orthogroup = basename(muscle_files) %>% str_split_fixed(string = ., "_",n=4)
orthoname = apply(orthogroup[,1:2],1,paste0,collapse="_")
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
r4s_res=pbmcapply::pbmcmapply(r4s_files,FUN = read.R4S, mc.cores = 14, SIMPLIFY=F) # DO NOT SIMPLIFY TO KEEP AS A LIST
saveRDS(r4s_res,file.path(fudir,"rate4site-4751_fungi-clades.rds"))


#### LOAD EVOLUTIONARY DATA ####
load(here::here("data/eggnog","4751_Fungi-eggnog-results.rdata"))
r4s_res = readRDS(file.path(fudir,"rate4site-4751_fungi-clades.rds"))
R4S = map(r4s_res, ~separate(data=.x, col="MSA", remove=F, into=c('nmsa','nseq'), sep="/"))

df_r4s = tibble(r4sfile  = names(r4s_res),
                clade = str_extract(dirname(r4sfile),clade_regex),
                OG = str_split_fixed(basename(r4sfile),'_',n=3)[,2],
                r4s_res = r4s_res,
                param = str_extract(basename(r4sfile),"\\.[^\\.]*trim[^\\.]*\\.") ) %>%
         unnest(r4s_res)


tmp = df_r4s %>%
    separate(col="MSA", remove=F, into=c('nmsa','nseq'), sep="/")
head(tmp)

         #mutate(f_MSA = nmsa/nseq)
df_rate = df_r4s %>% group_by(clade,OG,ID,param) %>%
          mutate(ER=mean_(SCORE),
                 f_gap = mean( !(SEQ %in% Biostrings::AA_STANDARD) )) %>%
          dplyr::select(OG,ID,clade,param,ER,f_MSA) %>%
          distinct() %>%
          pivot_wider(id_cols=c(OG,ID,param), names_from=clade, values_from = ER, names_prefix = 'ER.')

parse_fraction = function(x){ eval(parse(text=as.character(x))) }

tmp = df_r4s %>%
  rowwise() %>%
  mutate( f_aligned = parse_fraction(MSA),
          bin_aligned = cut(f_aligned, seq(0,1,len=11)))


df_rate_param= split( df_rate, df_rate$param)

dim(df_er_ppm)


#%>%
  group_by(clade,param) %>%
  summarise( notAA = mean( !(SEQ %in% Biostrings::AA_STANDARD) ))

ggplot(df_gap) +
  geom_col(aes(x=param,y=notAA,fill=param)) +
  facet_wrap(~clade) +
  theme_clean()


#### CALCULATE ER CORRELATION BETWEEN CLADES ####
COR_CLADE = map(df_rate_param, ~cor(.x[,-c(1:3)], method = "spearman",use='pairwise.complete'))
#N_CLADE = df_rate %>% select(starts_with("ER")) %>% pair_n

fungi=ape::read.tree(file=here::here("data/ncbi/ncbi-fungi.phy"))
df_fungi=readRDS(file.path(fudir,'4751_Fungi-tree_data.rds'))
#saveRDS(df_fungi, file.path(fudir,'4751_Fungi-tree_data.rds'))

df_yeast= df_fungi %>% filter(taxid %in% c(4932, 4896))
p1.1=ggtree(fungi,ladderize = T,branch.length = 'none') %<+% df_fungi +
  geom_text2(aes(subset=node %in% df_fungi$node, label=clade_desc),size=3,col='3CB371	',hjust=1,vjust=-1) +
  geom_nodepoint(aes(subset=node %in% df_fungi$node, label=clade_desc),size=3,col='3CB371	') +
  geom_tiplab(aes(subset=node %in% df_fungi$node),size=3) +
  geom_tippoint(aes(subset=node %in% df_yeast$node), col=c('dodgerblue','orange')) +
  xlim(-5,16)
p1.1

library(corrgram)
pcor=map(seq(COR_CLADE), ~
           ggcorrplot::ggcorrplot(COR_CLADE[[.x]],lab=T,method='circle',type='lower', tl.cex=6,
                                  show.legend = F,  ggtheme=theme_clean(), lab_size=3, title=names(COR_CLADE)[.x] ))

(pcor[[1]] + pcor[[2]] + pcor[[3]]) /
  (pcor[[4]] + pcor[[5]] + pcor[[6]])

align_patches(pcor) + plot_layout(ncol=3, nrow=2)
#p1.2=ggcorrplot::ggcorrplot(COR_ER_clade,
 #                           lab=T,
  #                          method = 'circle',type = 'lower')


p1 = p1.1 | p1.2

p1.3=GGally::ggpairs(df_rate[,-1])

ggsave(path=here::here("plots"),
       filename="4751_Fungi-evorate-clades.pdf",
       plot=p1,
       scale=3)
ggsave(path=here::here("plots"),
       filename="4751_Fungi-evorate-clades-scatterplot.pdf",
       plot=p1.3,
       scale=1.8)


#### CALCULATE CLADE ER vs. ABUNDANCE CORRELATION ####

# integrated paxdb datasets
ho_et_al_2018 = load.abundance() %>%
  set_names(c("orf","ho2018_MPC","ho2018_MDPC","ho2018_gfp","ho2018_ms")) %>%
  ungroup()
sc_abundance = get.paxdb(tax=4932,abundance = 'integrated',rm.zero=T)  %>%
              full_join(ho_et_al_2018,by = c('protid'='orf'))
df_er_ppm = left_join(df_rate, sc_abundance, by=c('ID'='protid'))


