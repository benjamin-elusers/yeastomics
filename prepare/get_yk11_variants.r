#install.packages('rjson')
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("src","function_YK11.r"))

S288C = load.sgd.proteome(withORF=T,rm.stop=F) # Reference SGD protein sequences
WT = get.positions(S288C) %>%
  group_by(orf) %>%
  mutate(bin.pos=dplyr::ntile(wt_pos,100))

yk11_seq_file = here("data","proteome-1011-strains.rds")
if( file.exists(yk11_seq_file) ){
  YK11 = readRDS()  # 1011 strains proteomes sequences
}else{
  yk11_seq_dir = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"
  YK11 = load.1011.strains(seqdir = yk11_seq_dir)
}
P = tibble( orf=names(YK11),
            n_strains=lengths(YK11),
            len = sapply(widths(YK11),unique)) %>%
    left_join(get.width(S288C), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
    mutate( match_wt = len == len.s288c )

yk11_var_file = here("data","YK11_VAR.rds")
if( file.exists(yk11_var_file) ){
  VAR = readRDS(var.file)
}else{
  VAR = purrr::map_df(YK11, get.variants,verbose=F)
}

SNP = get.variants.to.SNP(VAR) %>%
  add_count(id,ref_pos,name='nvar') %>%
  group_by(id,ref_pos) %>%
  mutate( alt_cumfr=sum(alt_fr)) %>%
  left_join(P,by=c('id'='orf'))

PROT_SNP = left_join(SNP,WT, by=c('id'='orf','ref_pos'='wt_pos','len.s288c'='len')) %>%
  mutate(wt_low = alt_aa == wt_aa,
         wt_missing=is.na(wt_aa),
         early_stop = (alt_aa == "*" & ref_pos != len),
         dSTI.ref=get.score.mutation(ref_aa,alt_aa),
         dSTI.wt=get.score.mutation(wt_aa,alt_aa))

snp_count_per_orf = PROT_SNP %>%
  group_by(id,len.s288c) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(ref_pos))
snp_count_per_orf %>% arrange(n_snp)

write_rds(PROT_SNP,here("data",'YK11-SNP.rds'))
write_rds(snp_count_per_orf,here("data",'YK11-ORF-VAR.rds'))
