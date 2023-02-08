library(tidyverse)
library(tidytree)
library(ggtree)
library(ggplot2)
library(ggthemes)
library(patchwork)
source(here::here("src","__setup_yeastomics__.r"))
fudir = here::here("data","eggnog","4751_Fungi")
mzdir = here::here("data","eggnog","33208_Metazoa")
fu_wexac="/media/WEXAC/EGGNOG/4751_Fungi"
mz_wexac="/media/WEXAC/EGGNOG/33208_Metazoa"
outdir="/data/benjamin/Evolution/EGGNOG/"
fu_out=file.path(outdir,"4751_Fungi")
mz_out=file.path(outdir,"33208_Metazoa")

count.fasta = function(fastafile){
  return(sum(grepl("^>",readLines(fastafile))))
}

find.file = function(directory, fileext=""){
  if(length(fileext)>1){
    resfile = list.files(directory, full.names=T, recursive=T)
    return( map(fileext, ~str_subset(resfile, pattern = paste0("\\.",.x,"$"))) )
  }else{
    resfile = list.files(directory, pattern=fileext, full.names=T, recursive=T)
    return(resfile)
  }
}

count.file = function(directory, fileext=""){
  return( map(find.file(directory, fileext),length) )
}

find.muscle = function(directory){
  cat(sprintf("get muscle alignments files %s...\n",basename(directory)))
  return( find.file(file.path(directory,"muscle"),"\\.mu$") )
}

find.trimal = function(directory,params=c(10,20,50,80,'pyout')){
  cat(sprintf("get trimal result files %s...\n",basename(directory)))
  params = paste0("trimal_gap",params)
  files = map(params, ~find.file(file.path(directory,.x), paste0("\\.",.x,"$"))) %>%
          set_names(nm = paste0("files_",params))
  n=map(files,n_distinct) %>% unique
  if( length(n) == 1){ return(files %>% as_tibble())  }
  return(files)
}

find.r4s = function(directory,params=c(10,20,50,80,'pyout')){
  cat(sprintf("get rate4site result files %s...\n",basename(directory)))
  trimal_params = paste0("trimal_gap",params) %>% str_replace("trimal_gap0","notrim")
  r4s_dir = str_replace(trimal_params,"trimal","r4s") %>% str_replace("notrim","r4s_muscle")

  r4s_ext="\\.eggnog_sptree\\.r4s_raw"
  files = map(seq(trimal_params),~find.file(directory = file.path(directory,r4s_dir[.x]),
                               fileext = paste0("\\.",trimal_params[.x],r4s_ext,"$"))) %>%
          set_names(nm = paste0("n_",trimal_params))
  n=map(files,n_distinct) %>% unique
  if( length(n) == 1){ return(files %>% as_tibble())  }
  return(files)
}

clade_regex="[0-9]+_[a-zA-Z]+_[0-9]+sp"
col_clade = paste0("clade_",c("id","name","ns"))

NODE_DIR = mz_out
clades = list.dirs(NODE_DIR,recursive = F) %>% basename %>% str_subset(clade_regex)
df_clade = tibble(clade_desc=clades) %>%
           separate(col='clade_desc', into=col_clade, sep = "_", remove=F) %>%
           rowwise() %>%
           mutate( resdir = file.path(NODE_DIR,clade_desc) )

muscle_files = pbmcapply::pbmcmapply(find.muscle,df_clade$resdir, mc.cores=14, SIMPLIFY=F)
trimal_files = pbmcapply::pbmcmapply(find.trimal,MoreArgs=list(params=c(10,20)), df_clade$resdir, mc.cores=14, SIMPLIFY=F )
r4s_files    = pbmcapply::pbmcmapply(find.r4s,MoreArgs=list(params=c(0,10,20)), df_clade$resdir, mc.cores=14, SIMPLIFY=F)

# purrr::flatten(muscle_files) %>%
#   str_detect(paste0("^",mz_out,"/.+/muscle/.+_1to1-orthologs\\.mu$")) %>%
#   mean

# trimal_files %>% purrr::flatten() %>% purrr::flatten() %>% unlist %>%
#   str_detect(paste0("^",mz_out,"/.+/trimal_gap[0-9]+/.+_1to1-orthologs\\.mu\\.trimal_gap[0-9]+$")) %>%
#   mean

test= r4s_files %>% purrr::flatten() %>% purrr::flatten() %>% unlist
test %>%
  str_detect(paste0("^",mz_out,"/.+/r4s_.+/.+_1to1-orthologs.+\\.eggnog_sptree\\.r4s_raw$")) %>%
  mean

df_clade$muscle_files = muscle_files
df_clade$trimal_files = trimal_files
df_clade$r4s_files = r4s_files

df_clade$n_muscle = map_int(df_clade$muscle_files,n_distinct)
df_clade$n_trimal = map_dfr(df_clade$trimal_files, ~map_df(.x,n_distinct))
df_clade$n_r4s = map_dfr(df_clade$r4s_files, ~map_df(.x,n_distinct))

cbind(df_clade$clade_desc, df_clade$n_muscle - df_clade$n_trimal)
cbind(df_clade$clade_desc, df_clade$n_muscle - df_clade$n_r4s)

errIO = pbmcapply::pbmclapply(unlist(r4s_files),mc.cores = 14,
                              function(x){ if( length(readLines(x,n=10))==0 ){ return(x) }}) %>%
        purrr::compact() %>% unlist()


if(NODE_DIR==fu_out){
  saveRDS(df_clade,file=here::here("data/eggnog","4751_Fungi-eggnog-results.rds"))
}else if(NODE_DIR==mz_out){
  saveRDS(df_clade,file=here::here("data/eggnog","33208_Metazoa-eggnog-results.rds"))
}
#save.image(here::here("data/eggnog","4751_Fungi-eggnog-results.rdata"))
#
# clade_desc = dirname(muscle_files) %>% str_extract(clade_regex)
# orthogroup = basename(muscle_files) %>% str_split_fixed(string = ., "_",n=4)
# orthoname = apply(orthogroup[,1:2],1,paste0,collapse="_")
# muscle_gaps = pbmcapply::pbmcmapply(count_gaps_fromfile, muscle_files,  mc.cores = 14)
# names(muscle_gaps) = orthoname
# muscle_gapstat = tibble(ortho=orthoname,
#                         clade=clade_desc,
#                         gaps=muscle_gaps) %>%
#                rowwise() %>%
#                mutate( agap = mean(gaps,na.rm=T),
#                        sgap = sd(gaps,na.rm=T),
#                        mgap = median(gaps,na.rm=T),
#                        gap5 = quantile(gaps,probs=0.05),
#                        gap20 = quantile(gaps,probs=0.2),
#                        gap80 = quantile(gaps,probs=0.8),
#                        gap95 = quantile(gaps,probs=0.95)) %>%
#               mutate(og=str_split_fixed(ortho,"_",2)[,2]) %>%
#               group_by(clade) %>% add_count(name='nclade')
#
# clade_desc = dirname(trimal_files) %>% str_extract(clade_regex)
# orthogroup = basename(trimal_files) %>% str_split_fixed(string = ., "_",n=4)
# orthoname = apply(orthogroup[,1:2],1,paste0,collapse="_")
# trimal_param = dirname(trimal_files) %>% basename()
#
# trim_gaps = pbmcapply::pbmcmapply(trimal_files, FUN=count_gaps_fromfile, mc.cores = 14)
# names(trim_gaps) = orthoname
# trim_gapstat = tibble(ortho=orthoname,
#                       clade=clade_desc,
#                       params=trimal_param,
#                       gaps=trim_gaps) %>%
#   rowwise() %>%
#   mutate( agap = mean(gaps,na.rm=T),
#           sgap = sd(gaps,na.rm=T),
#           mgap = median(gaps,na.rm=T),
#           gap5 = quantile(gaps,probs=0.05),
#           gap20 = quantile(gaps,probs=0.2),
#           gap80 = quantile(gaps,probs=0.8),
#           gap95 = quantile(gaps,probs=0.95)) %>%
#   mutate(og=str_split_fixed(ortho,"_",2)[,2]) %>%
#   group_by(clade,params) %>% add_count(name='nclade')
#
# pg=ggplot(muscle_gapstat) +
#   geom_density(aes(x=agap)) +
#   geom_text(aes(label=paste0('n=',nclade)),x=Inf,y=Inf,hjust='inward',vjust='inward',check_overlap = T,size=3) +
#   geom_density(data=trim_gapstat,aes(x=agap,color=clade_desc)) +
#   facet_wrap(~params,scales='free_y',nrow = 4) +
#   xlab('Average % gap in alignment') +
#   theme_clean()
# pg
#
# ggsave(path=here::here("plots"),
#        filename="4751_Fungi-gaps_frequency_eggnog.pdf",
#        plot=pg,
#        scale=1.5)

###
sgd_len = get.width(load.sgd.proteome()) %>% rename(s288c_len=len)

##### Wapinski Fungi lineage ------------------------------------------------------------
wapinski_dir = "/media/WEXAC/FUNGI/"
wapinski.rds =  here('output','sc-evorate-wapinski.rds')
wapinski.fasta = Rfast::read.directory(file.path(wapinski_dir,'fasta')) %>%
                 str_subset('\\.fasta$') %>% file.path(wapinski_dir,'fasta',.)
wapinski.msa = load_msa(wapinski.fasta,ref = 'Saccharomyces_cerevisiae')
wapinski.r4sfiles = Rfast::read.directory(file.path(wapinski_dir,'R4S')) %>%
               str_subset('raw\\.r4s$') %>% file.path(wapinski_dir,'R4S',.)
wapinski.r4s = load_r4s(wapinski.r4sfiles)

wapinski.evo = inner_join(wapinski.msa,wapinski.r4s,
                         by=c('id'='ID','msa_pos'='POS','ref_aa'='SEQ')) %>%
               group_by(id) %>% mutate( len_ref = max_(ref_pos), len_msa = max_(msa_pos)) %>%
               #dplyr::filter(!is.na(ref_pos)) %>%
               dplyr::rename(r4s_rate=SCORE) %>%
               dplyr::select(-c('QQ1','QQ2','STD','MSA')) %>%
               left_join(sgd_len, by=c('id'='orf'))


wapinski.evo[wapinski.evo$len_ref != wapinski.evo$s288c_len,]


YLR158C.fa=wapinski.fasta %>% str_subset("YLR158C")
Biostrings::readAAMultipleAlignment(YLR158C.fa)

YLR158C.msa = wapinski.msa %>% filter( id == "YLR158C")
YLR158C.r4s = wapinski.r4s %>% filter( ID == "YLR158C") %>% print(n=500)

YLR158C.fa$ref_aa
YLR158C.r4s$SEQ

##### Cerevisiae isolates ------------------------------------------------------
strains_dir = "/media/WEXAC_data/1011G/"
strains.rds =  here('output','sc-evorate-yk11.rds')

sc_evo$yk11 = preload(strains.rds,
                      load.evorate(alndir=file.path(strains_dir,'aln_s288c'), id_type = 'ORF', ext.seq = 'fasta', ref = 'S288C'),
                      'get evolutionary rates for 1011 isolated yeast...')

r4s_res=pbmcapply::pbmcmapply(r4s_files,FUN = read.R4S, mc.cores = 20, SIMPLIFY=F) # DO NOT SIMPLIFY TO KEEP AS A LIST
saveRDS(r4s_res,file.path(fudir,"rate4site-4751_fungi-clades.rds"))

#### LOAD EVOLUTIONARY DATA ####
load(here::here("data/eggnog","4751_Fungi-eggnog-results.rdata"))
r4s_res = readRDS(file.path(fudir,"rate4site-4751_fungi-clades.rds"))

#tens=seq(0,1,len=6)
brk = c(0,0.3,0.5,0.7,1)
deciles = brk %>% percent0 %>% paste.pair(s='-')

df_r4s = tibble(r4sfile = names(r4s_res),
                clade = str_extract(dirname(r4sfile),clade_regex),
                OG = str_split_fixed(basename(r4sfile),'_',n=3)[,2],
                r4s_res = r4s_res,
                param = str_extract(basename(r4sfile),"\\.[^\\.]*trim[^\\.]*\\.") ) %>%
         unnest(r4s_res) %>%
         mutate(fmsa_bin = cut(fmsa, brk, deciles) )

df_rate = df_r4s %>% group_by(clade,OG,ID,param) %>%
          mutate(ER=mean_(SCORE),
                 gap_first = mean_( !(SEQ %in% Biostrings::AA_STANDARD) )) %>%
          dplyr::select(OG,ID,clade,param,ER,gap_first) %>%
          distinct() %>%
          pivot_wider(id_cols=c(OG,ID,param), names_from=clade, values_from = ER, names_prefix = 'ER.')

df_rate_param= split( df_rate, df_rate$param)

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

df_er_ppm = left_join(df_rate, sc_abundance, by=c('ID'='protid')) %>% ungroup
dim(df_er_ppm)

cor_ppm = map_dfr(clades,~cor.sub.by(df_er_ppm, XX=paste0("ER.",.x), YY="ppm_int", BY="param") %>% mutate(xcol=paste0("ER.",.x),ycol="ppm_int"))
cc1=ggplot(cor_ppm, aes(x=param)) +
    geom_col(aes(y=r,fill=param)) +
    geom_text(aes(y=r,label=r %>% round(2)),size=3,vjust='outward') +
    facet_wrap(~xcol) + ylab("spearman (ER-PPM)") +
    theme_clean() +
    theme(axis.text.x = element_text(angle=30,hjust=1)) +
    scale_fill_metro_d()

cor_mpc = map_dfr(clades,~cor.sub.by(df_er_ppm, XX=paste0("ER.",.x), YY="ho2018_MPC", BY="param") %>% mutate(xcol=paste0("ER.",.x),ycol="ppm_int"))
cc2=ggplot(cor_mpc, aes(x=param)) +
    geom_col(aes(y=r,fill=param)) +
    geom_text(aes(y=r,label=r %>% round(2)),size=3,vjust='outward') +
    facet_wrap(~xcol)  + ylab("spearman (ER-MPC)") +
    theme_clean() +
    theme(axis.text.x = element_text(angle=30,hjust=1))+
    scale_fill_metro_d()

cor_mdpc = map_dfr(clades,~cor.sub.by(df_er_ppm, XX=paste0("ER.",.x), YY="ho2018_MDPC", BY="param") %>% mutate(xcol=paste0("ER.",.x),ycol="ppm_int"))
cc3=ggplot(cor_mdpc, aes(x=param)) +
    geom_col(aes(y=r,fill=param)) +
    geom_text(aes(y=r,label=r %>% round(2)),size=3,vjust='outward') +
    facet_wrap(~xcol) + ylab("spearman (ER-MDPC)") +
    theme_clean() +
    theme(axis.text.x = element_text(angle=30,hjust=1)) +
    scale_fill_fermenter()

ggsave(
  path=here::here("plots","paper_evo"),
  filename = "4751_Fungi-evorate_abundance-clades.pdf",
  plot = marrangeGrob(list(cc1,cc2,cc3), nrow=1, ncol=1),
  width = 21, height = 29.7, units = "cm"
)
