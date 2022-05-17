source("https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/__setup_yeastomics__.r")

evo_snp = load.evorate()
evo_snp_prot = 
  group_by(evo_snp, id) %>% 
  summarize( len_ref = max(ref_pos), len_msa=max(msa_pos),
             n_strains = max(matched), f_strains = max(matched)/1012, n_mismatched = sum(mismatched!=0), 
             r4s=mean(r4s_rate), 
             iq=mean(iq_rate),iq_ml=mean(iq_mlrate), iq_cat=mean(iq_cat), iq_hi = mean(iq_rate_hicat), 
             leisr_mle = mean(leisr_mle[leisr_mle!=0]), 
             leisr_low=mean(leisr_low[leisr_low!=0]), 
             leisr_up=mean(leisr_up[leisr_up!=0]), leisr_global=mean(leisr_global), leisr_local=mean(leisr_local) )

write_delim(evo_snp,here::here('output','evolution-snp-residue.tsv'),delim = '\t')
write_delim(evo_snp_prot,here::here('output','evolution-snp-protein.tsv'),delim = '\t')

evo_full = load.evorate(resdir="/media/WEXAC_data/FUNGI/",ref='Saccharomyces_cerevisiae',ID="ORF")

evo_full_prot = 
  group_by(evo_full, id) %>% 
  summarize( len_ref = max(ref_pos), len_msa=max(msa_pos),
             f_mismatched = mean(mismatched!=0), f_indel = mean(indel!=0),
             r4s=mean(r4s_rate), 
             iq=mean(iq_rate),iq_ml=mean(iq_mlrate), iq_cat=mean(iq_cat), iq_hi = mean(iq_rate_hicat), 
             leisr_mle = mean(leisr_mle), leisr_low=mean(leisr_low), leisr_up=mean(leisr_up),
             leisr_global=mean(leisr_global), leisr_local=mean(leisr_local) )


write_delim(evo_full,here::here('output','evolution-fungi-residue.tsv'),delim = '\t')
write_delim(evo_full_prot,here::here('output','evolution-fungi-protein.tsv'),delim = '\t')

