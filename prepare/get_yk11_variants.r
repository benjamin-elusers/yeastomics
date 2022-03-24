#install.packages('rjson')
source(here::here("src","__setup_yeastomics__.r"))
source(here::here("src","function_YK11.r"))

S288C = load.sgd.proteome(withORF=T,rm.stop=F) # Reference SGD protein sequences
WT = get.positions(S288C) %>%
  group_by(orf) %>%
  mutate(bin.pos=dplyr::ntile(wt_pos,100))

yk11_seq_file = here("data","proteome-1011-strains.rds")
if( file.exists(yk11_seq_file) ){
  YK11 = readRDS(yk11_seq_file)  # 1011 strains proteomes sequences
}else{
  yk11_seq_dir = "/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/"
  YK11 = load.1011.strains(seqdir = yk11_seq_dir)
}
P = tibble( orf=names(YK11),
            n_strains=lengths(YK11),
            len = widths(YK11)) %>%
    left_join(get.width(S288C), by=c('orf'='orf'),suffix=c('','.s288c')) %>%
    mutate( match_wt = len == len.s288c )

yk11_var_file = here("data","YK11_VAR.rds")
if( file.exists(yk11_var_file) ){
  VAR = readRDS(yk11_var_file)
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
  group_by(id,len,len.s288c,n_strains) %>%
  summarize(n_snp=sum(nvar),n_var=n_distinct(ref_pos))
snp_count_per_orf %>% arrange(n_snp)

write_rds(PROT_SNP,here("data",'YK11-SNP.rds'))
write_rds(snp_count_per_orf,here("data",'YK11-ORF-VAR.rds'))

PROT_SNP=read_rds(here("data",'YK11-SNP.rds'))
snp_count_per_orf=read_rds(here("data",'YK11-ORF-VAR.rds'))

### Must have run align_s288c_to_1011strains.r !!!
# Filter orf with best match to s288c from all strains
df_strains_s288c = read_rds(here("prepare","proteome_aligned_s288c_vs_1011strains.rds"))
#df_best = df_strains_s288c %>% group_by(ID1) %>% filter(SCORE.BLOSUM100 == max(SCORE.BLOSUM100))

proteomics_strains = c('AMH','CQC','BTT','CPI','BED','BAN','BPL','CMP')
best_strain = df_strains_s288c %>%
  group_by(strain) %>%
  summarize(aligned = sum(S),
            norf = n_distinct(ID1),
            proteome_len = sum(L),
            pid_proteome = 100*(aligned/proteome_len),
            nsnp = sum(NS),
            orf_with_snv = sum(NS != 0 & (I+D) == 0),
            orf_with_indel = sum((I+D) > 0)
            ) %>%
  mutate(is_proteomics = strain %in% proteomics_strains) %>%
  arrange(nsnp,desc(pid_proteome),desc(norf))%>%
  left_join(load.peter2018.data(1), by=c('strain'='standardized_name'))

write_rds(best_strain, here("prepare","best_strain_aligned_to_s288c_proteome.rds"))

# Show
median_snp = median(best_strain$nsnp)
ggplot(best_strain,aes(fill=is_proteomics)) +
  geom_col(aes(y=reorder(strain,pid_proteome),x=nsnp),orientation='y') +
  geom_vline(xintercept=median_snp,linetype='dashed') +
  facet_wrap(~ecological_origins, scales = 'free_y') + xlab('# of SNP') +
  scale_y_discrete(guide=guide_axis(n.dodge=2))+
  scale_fill_manual(values=c("TRUE"='red',"FALSE"='gray50')) +
  theme(legend.position = 'none', axis.text =  element_text(size = 5))

