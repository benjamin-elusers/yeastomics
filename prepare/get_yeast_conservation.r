source("https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/__setup_yeastomics__.r")

sc_identifiers = load.annotation()
evo_snp = preload( here('data','evorate-strains-snp.rds') ,load.evorate())
#resdir="/media/WEXAC_data/1011G/"
#ext.r4s = 'raw.r4s'
#ref='S288C'
#ID="ORF"
#ncores=parallelly::availableCores(which='max')-2

evo_snp_prot =  group_by(evo_snp,id) %>%
  summarize( len_ref = max(ref_pos), len_msa=max(msa_pos),
             n_strains = max(matched), f_strains = max(matched)/max(total),
             n_matched = sum(matched==n_strains), n_mismatched = sum(mismatched!=0), n_indel = sum(indel!=0 ),
             f_matched = mean(matched==n_strains), f_mismatched = mean(mismatched!=0), f_indel = mean(indel!=0 ),
             r4s=mean(r4s_rate),
             iq=mean(iq_rate),iq_ml=mean(iq_mlrate), iq_cat=mean(iq_cat), iq_hi = mean(iq_rate_hicat),
             leisr_mle = mean(leisr_mle[leisr_mle!=0]),
             leisr_low=mean(leisr_low[leisr_low!=0]),
             leisr_up=mean(leisr_up[leisr_up!=0]), leisr_global=mean(leisr_global), leisr_local=mean(leisr_local) )

write_delim(evo_snp,here::here('output','evolution-snp-residue.tsv'),delim = '\t')
write_delim(evo_snp_prot,here::here('output','evolution-snp-protein.tsv'),delim = '\t')

evo_full = preload(here('data','evorate-fungi-orthologs.rds'),
                    load.evorate(resdir="/media/WEXAC_data/FUNGI/",ref='Saccharomyces_cerevisiae',ID="ORF"))

evo_full_prot = group_by(evo_full, id) %>%
  summarize( len_ref = max(ref_pos), len_msa=max(msa_pos),
             n_species = max(matched), f_species = max(matched)/max(total),
             n_matched = sum(matched==n_species), n_mismatched = sum(mismatched!=0),  n_indel = sum(indel!=0),
             f_matched = mean(matched==n_species), f_mismatched = mean(mismatched!=0), f_indel = mean(indel!=0),
             r4s=mean(r4s_rate),
             iq=mean(iq_rate),iq_ml=mean(iq_mlrate), iq_cat=mean(iq_cat), iq_hi = mean(iq_rate_hicat),
             leisr_mle = mean(leisr_mle), leisr_low=mean(leisr_low), leisr_up=mean(leisr_up),
             leisr_global=mean(leisr_global), leisr_local=mean(leisr_local) )


write_delim(evo_full,here::here('output','evolution-fungi-residue.tsv'),delim = '\t')
write_delim(evo_full_prot,here::here('output','evolution-fungi-protein.tsv'),delim = '\t')


evo_yeast = left_join(evo_snp_prot,evo_full_prot, by=c('id','len_ref'),suffix=c('.yk11','.fungi')) %>%
  mutate(HAS_ORTHOLOG = !is.na(len_msa.fungi) ) %>%
  left_join(sc_identifiers,by=c('id'='ORF')) %>%
  dplyr::mutate( f_snp = n_mismatched/len_msa.yk11, pid.fungi=1-f_mismatched) %>%
  dplyr::rename(orf=id,n_snp = n_mismatched) %>%
  dplyr::select(-f_mismatched) %>%
  relocate(orf,UNIPROT,GENENAME,SGD,HAS_ORTHOLOG, len_ref,
           len_msa.yk11, n_snp,f_snp, len_msa.fungi,pid.fungi) %>%
  dplyr::select(-paste0(er_fungi_worst,'.fungi'),-paste0(er_strains_worst,'.yk11'))

