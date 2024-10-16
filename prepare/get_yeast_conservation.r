source("https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/__setup_yeastomics__.r")
sgd_len = get.width(load.sgd.proteome()) %>% rename(s288c_len=len)
evodir = "/data/benjamin/Evolution"

##### S. cerevisiae isolates ---------------------------------------------------
yk11_dir = file.path(evodir,"YK11")
#yk11.rds =  here('output','evorate-yk11-msa-r4s.rds')
yk11_rds =  file.path(yk11_dir,'evorate-yk11-msa-r4s.rds')

if(!file.exists(yk11_rds)){

  yk11.fasta = find.fasta(file.path(yk11_dir,'aln_s288c'),"\\.fasta$")

  yk11.seq = load_seq(yk11.fasta,ref ='S288C', id_type = 'ORF', ncores = 14)
  yk11.msa = load_msa(yk11.seq,ref = 'S288C')

  yk11.r4sfiles = Rfast::read.directory(file.path(yk11_dir,'R4S')) %>%
                  str_subset('raw\\.r4s$') %>% file.path(yk11_dir,'R4S',.)
  yk11.r4s = load_r4s(yk11.r4sfiles)

  # yk11.leisrfiles = Rfast::read.directory(file.path(yk11_dir,'LEISR')) %>%
  #   str_subset('LEISR\\.json$') %>% file.path(yk11_dir,'LEISR',.)
  #yk11.leisr = load_leisr(yk11.leisrfiles)
  # yk11.iqtreefiles = Rfast::read.directory(file.path(yk11_dir,'IQTREE')) %>%
  #   str_subset('\\.mltree$') %>% file.path(yk11_dir,'IQTREE',.)
  #yk11.iqtree = load_iqtree(yk11.iqtreefiles)

  yk11.evo = inner_join(yk11.msa,yk11.r4s,
                            by=c('id'='ID','msa_pos'='POS')) %>%
    group_by(id) %>% mutate( len_ref = max_(ref_pos), len_msa = max_(msa_pos),
                             fid = matched/total, fmis=mismatched/total) %>%
    dplyr::rename(r4s_rate=SCORE, r4s_aa=SEQ, r4s_nmsa=nmsa, r4s_fmsa=fmsa, r4s_fgap=fgap, r4s_nseq=nseq) %>%
    dplyr::select(-c('QQ1','QQ2','STD','MSA')) %>%
    left_join(sgd_len, by=c('id'='orf')) %>%
    relocate(id,msa_pos,ref_pos,ref_gap,ref_aa,r4s_aa,
             matched,fid, mismatched,fmis, indel,ins,del, total,
             r4s_nmsa, r4s_fmsa, r4s_fgap, r4s_nseq)

  yk11.data = list()
  yk11.data[["seq"]] = yk11.seq
  yk11.data[["msa"]] = yk11.msa
  yk11.data[["r4s"]] = yk11.r4s
  yk11.data[["evo"]] = yk11.evo

  saveRDS(yk11.data, yk11.rds)
  saveRDS(yk11.data, yk11_rds)
}else{
  cat('reading precomputed yk11 data...\n')
  yk11.data = readRDS(yk11.rds)
}

##### Wapinski Fungi lineage ---------------------------------------------------
wapinski.rds =  here('output','evorate-wapinski-msa-r4s.rds')

if(!file.exists(wapinski.rds)){

  wapinski_dir = file.path(evodir,"WAPINSKI")
  wapinski.fasta = Rfast::read.directory(file.path(wapinski_dir,'fasta')) %>%
                   str_subset('\\.fasta$') %>% file.path(wapinski_dir,'fasta',.)

  wapinski.seq = load_seq(wapinski.fasta,ref ='Saccharomyces_cerevisiae', id_type = 'ORF')
  wapinski.msa = load_msa(wapinski.seq,ref = 'Saccharomyces_cerevisiae')
  wapinski.r4sfiles = Rfast::read.directory(file.path(wapinski_dir,'R4S')) %>%
                   str_subset('raw\\.r4s$') %>% file.path(wapinski_dir,'R4S',.)
  wapinski.r4s = load_r4s(wapinski.r4sfiles)

  wapinski.evo = inner_join(wapinski.msa,wapinski.r4s,
                            by=c('id'='ID','msa_pos'='POS')) %>%
                 group_by(id) %>% mutate( len_ref = max_(ref_pos), len_msa = max_(msa_pos),
                                          fid = matched/total, fmis=mismatched/total) %>%
                 dplyr::rename(r4s_rate=SCORE, r4s_aa=SEQ, r4s_nmsa=nmsa,
                               r4s_fmsa=fmsa, r4s_fgap=fgap, r4s_nseq=nseq) %>%
                 dplyr::select(-c('QQ1','QQ2','STD','MSA')) %>%
                 left_join(sgd_len, by=c('id'='orf')) %>%
                 relocate(id,msa_pos,ref_pos,ref_gap,ref_aa,r4s_aa,
                         matched,fid, mismatched,fmis, indel,ins,del, total,
                         r4s_nmsa, r4s_fmsa, r4s_fgap, r4s_nseq)

  wapinski.data = list()
  wapinski.data[["seq"]] = wapinski.seq
  wapinski.data[["msa"]] = wapinski.msa
  wapinski.data[["r4s"]] = wapinski.r4s
  wapinski.data[["evo"]] = wapinski.evo

  saveRDS(wapinski.data, wapinski.rds)
  wapinski_rds =  file.path(wapinski_dir,wapinski.rds)
  saveRDS(wapinski.data, wapinski_rds)
}else{
  cat('reading precomputed wapinski data...\n')
  wapinski.data = readRDS(wapinski.rds)
}
##### Eggnog Fungi lineage -----------------------------------------------------
fungi.rds =  here('output','evorate-eggnog_fungi-msa-r4s.rds')

fu_dir = file.path(evodir,"EGGNOG","4751_Fungi")
clades_dir = list.dirs(fu_dir,recursive = F,full.names = F) %>% str_subset("sp$")

fungi.data=list()
i=1

for( clade in clades_dir){
  cat(i,")",clade,"\n")
  clade_rds =  file.path(fu_dir,sprintf('evorate-eggnogV5_fungi-%s-msa-r4s.rds',clade))
  if(!file.exists(clade_rds)){

    fasta_dir = file.path(fu_dir,clade,'fasta')
    fasta_ext = '\\.fa$'
    clade.fasta = list.files(fasta_dir, fasta_ext, full.names=T, recursive=F)
    ncores=14
    message(sprintf("Retrieving orthologs identifiers by order of reference using %s cpus",ncores))
    clade.ids = pbmcapply::pbmclapply(FUN=get.fasta.names, clade.fasta, mc.cores=ncores) %>%
                purrr::map(.,pluck,1) %>%
                set_names(nm=basename(clade.fasta) %>%str_replace('.fa$','.mu'))

    muscle_dir=file.path(fu_dir,clade,'muscle')
    muscle_ext='\\.mu$'
    clade.ali = list.files(muscle_dir,muscle_ext,full.names=T, recursive = F)

    seqnames = intersect(basename(clade.ali),names(clade.ids))
    clade.seq = load_seq(clade.ali,ref = clade.ids, id_type = 'filename')
    clade.msa = load_msa(clade.seq[seqnames],ref=clade.ids[seqnames], id_type = 'filename')

    r4s_dir = file.path(fu_dir,clade,'r4s_muscle')
    r4s_ext = '\\.r4s_raw$'
    clade.r4sfiles = list.files(r4s_dir,r4s_ext,full.names=T, recursive=F)
    clade.r4s = load_r4s(clade.r4sfiles) %>%
                 mutate( IDFILE = str_replace(IDFILE,"(?<=\\.mu).+$",""))

    clade.evo = inner_join(clade.msa,clade.r4s,
                           by=c('idfile'='IDFILE','msa_pos'='POS')) %>%
                separate(sep = "_", into = c('clade_name','OG','protid'), remove=F, col = 'idfile') %>%
                group_by(clade_name,OG,id) %>%
                mutate( len_ref = max_(ref_pos), len_msa = max_(msa_pos),
                        fid = matched/total, fmis=mismatched/total) %>%
                dplyr::rename(r4s_rate=SCORE,r4s_aa=SEQ, r4s_nmsa=nmsa,
                              r4s_fmsa=fmsa, r4s_fgap=fgap, r4s_nseq=nseq) %>%
                dplyr::select(-c('QQ1','QQ2','STD','MSA')) %>%
                relocate(idfile,clade_name,OG,
                         protid,id,ID,
                         msa_pos,ref_pos,ref_gap,ref_aa,r4s_aa,
                         matched,fid, mismatched,fmis, indel,ins,del, total,
                         r4s_nmsa, r4s_fmsa, r4s_fgap, r4s_nseq)

    clade.data = list()
    clade.data[[clade]][["seq"]] = clade.seq
    clade.data[[clade]][["msa"]] = clade.msa
    clade.data[[clade]][["r4s"]] = clade.r4s
    clade.data[[clade]][["evo"]] = clade.evo

    saveRDS(clade.data[[clade]],clade_rds)
  }else{
    cat('--> Reading precomputed data...\n\n')
    clade.data = readRDS(clade_rds)
    #fungi.data
  }
  i=i+1
}

##### Eggnog Metazoa lineage ---------------------------------------------------
metazoa.rds =  here('output','evorate-eggnog_metazoa-msa-r4s.rds')

mz_dir = "/data/benjamin/Evolution/EGGNOG/33208_Metazoa"
clades_dir = list.dirs(mz_dir,full.names = F, recursive = F) %>% str_subset("sp$")

CLADES_DIR = clades_dir[clades_dir %pin% c('Metazoa','Mammalia','Euarchontoglires','Laurasiatheria')]

#i=str_which(clades_dir,'Mammalia') #

metazoa.data=list()
i=1
#clade = clades_dir[i]
for( clade in CLADES_DIR){
  cat(i,")",clade,"\n")
  clade_rds =  file.path(mz_dir,sprintf('evorate-eggnogV5_metazoa-%s-msa-r4s.rds',clade))

  if(!file.exists(clade_rds)){

    fasta_dir = file.path(mz_dir,clade,'fasta')
    fasta_ext = '\\.fa$'
    clade.fasta = list.files(fasta_dir, fasta_ext, full.names=T, recursive=F)
    ncores=14
    message(sprintf("Retrieving orthologs identifiers by order of reference using %s cpus",ncores))
    clade.ids = pbmcapply::pbmcmapply(FUN=get.fasta.names, clade.fasta, mc.cores=14,SIMPLIFY=F) %>%
      purrr::map(.,pluck,1) %>%
      set_names(nm=basename(clade.fasta) %>%str_replace('.fa$','.mu'))
    saveRDS(clade.ids,file.path(mz_dir,sprintf('%s-ids.rds',clade)))

    muscle_dir=file.path(mz_dir,clade,'muscle')
    muscle_ext='\\.mu$'
    clade.ali = list.files(muscle_dir,muscle_ext,full.names=T, recursive = F)

    seqnames = intersect(basename(clade.ali),names(clade.ids))
    clade.seq = load_seq(clade.ali,ref = clade.ids, id_type = 'filename')
    saveRDS(clade.seq,file.path(mz_dir,sprintf('%s-seq.rds',clade)))

    clade.msa = load_msa(clade.seq[seqnames],ref=clade.ids[seqnames], id_type = 'filename')
    saveRDS(clade.msa,file.path(mz_dir,sprintf('%s-msa.rds',clade)))
    rm(clade.seq)
    rm(clade.msa)
    gc()

    r4s_dir = file.path(mz_dir,clade,'r4s_muscle')
    r4s_ext = '\\.r4s_raw$'
    clade.r4sfiles = list.files(r4s_dir,r4s_ext,full.names=T, recursive=F)
    errIO = pbmcapply::pbmclapply(clade.r4sfiles, mc.cores=14,
                                  function(x){ if( length(readLines(x,n=10))==0 ){ return(x) }}) %>%
            purrr::compact() %>% unlist() %>% unlink()
    clade.r4sfiles = list.files(r4s_dir,r4s_ext,full.names=T, recursive=F)

    clade.r4s = load_r4s(clade.r4sfiles,ncores = 10) %>% # 9426
      mutate( IDFILE = str_replace(IDFILE,"(?<=\\.mu).+$",""))
      #group_by(IDFILE) %>% nest(col=-IDFILE)
    saveRDS(clade.r4s,file.path(mz_dir,sprintf('%s-r4s_muscle.rds',clade)))

    clade.msa=readRDS(file.path(mz_dir,sprintf('%s-msa.rds',clade)))
    clade.seq=readRDS(file.path(mz_dir,sprintf('%s-seq.rds',clade)))

    clade.evo = inner_join(clade.msa,clade.r4s,
                           by=c('idfile'='IDFILE','msa_pos'='POS')) %>%
      separate(sep = "_", into = c('clade_name','OG','protid'), remove=F, col = 'idfile') %>%
      group_by(clade_name,OG,id) %>%
      mutate( len_ref = max_(ref_pos), len_msa = max_(msa_pos),
              fid = matched/total, fmis=mismatched/total) %>%
      dplyr::rename(r4s_rate=SCORE,r4s_aa=SEQ, r4s_nmsa=nmsa,
                    r4s_fmsa=fmsa, r4s_fgap=fgap, r4s_nseq=nseq) %>%
      dplyr::select(-c('QQ1','QQ2','STD','MSA')) %>%
      relocate(idfile,clade_name,OG,
               protid,id,ID,
               msa_pos,ref_pos,ref_gap,ref_aa,r4s_aa,
               matched,fid, mismatched,fmis, indel,ins,del, total,
               r4s_nmsa, r4s_fmsa, r4s_fgap, r4s_nseq)
    clade.data = list()
    clade.data[[clade]][["seq"]] = clade.seq
    clade.data[[clade]][["msa"]] = clade.msa
    clade.data[[clade]][["r4s"]] = clade.r4s
    clade.data[[clade]][["evo"]] = clade.evo

    saveRDS(clade.data[[clade]],clade_rds)
    rm(clade.data,clade.seq,clade.msa,clade.r4s,clade.evo)
    gc()
  }else{
    cat('--> Reading precomputed data...\n\n')
    #clade.data = readRDS(clade_rds)
    #metazoa.data
  }
  i=i+1
}

## ANNOTATION ------------------------------------------------------------------
sc_annotation = load.annotation()
sc_identifiers = sc_annotation %>% dplyr::select(UNIPROT,ORF,GENENAME,SGD,OG) %>%
                  dplyr::filter(!duplicated(ORF) & !duplicated(UNIPROT) & !duplicated(SGD) & !duplicated(GENENAME))
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
saveRDS(evo_snp_prot,here::here('output','evolution-snp-protein.rds'))

write_delim(evo_snp,gzfile(here::here('output','evolution-snp-residue.tsv.gz')),delim = '\t')
write_delim(evo_snp_prot,gzfile(here::here('output','evolution-snp-protein.tsv.gz')),delim = '\t')

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
saveRDS(evo_full_prot,here::here('output','evolution-fungi-protein.rds'))
write_delim(evo_full,gzfile(here::here('output','evolution-fungi-residue.tsv.gz')),delim = '\t')
write_delim(evo_full_prot,gzfile(here::here('output','evolution-fungi-protein.tsv.gz')),delim = '\t')


## CORRELATION EVORATE
fungi_rate = evo_full_prot %>% dplyr::select(r4s:leisr_local) %>% as.matrix
strains_rate = evo_snp_prot %>% dplyr::select(r4s:leisr_local) %>% as.matrix
er_fungi_cor = cor(fungi_rate,use='pairwise.complete',met='spearman')
er_strains_cor = cor(strains_rate,use='pairwise.complete',met='spearman')

library(ggcorrplot)
p_fungi = ggcorrplot(er_fungi_cor,type='upper',method='circle',lab = T,lab_size = 3, title='fungi evo. rate',ggtheme = theme_classic())
p_snp = ggcorrplot(er_strains_cor,type='upper',method='circle',lab = T,lab_size = 3, title='strains evo. rate',ggtheme = theme_classic())
ggsave(p_fungi, path = here('plots'),filename='cor-evolution-fungi.png',scale=1.2)
ggsave(p_snp, path = here('plots'),filename='cor-evolution-snp.png',scale=1.2)
er_fungi_worst   = colnames(er_fungi_cor)[ abs(er_fungi_cor[1,]) < 0.7 ]
er_strains_worst = colnames(er_strains_cor)[ abs(er_strains_cor[1,]) < 0.7 ]

evo_yeast = left_join(evo_full_prot,evo_snp_prot, by=c('id','len_ref'),suffix=c('.fungi','.yk11')) %>%
  mutate(HAS_ORTHOLOG = !is.na(len_msa.fungi) ) %>%
  left_join(sc_identifiers,by=c('id'='ORF')) %>%
  dplyr::mutate( f_snp = n_mismatched.yk11/len_msa.yk11, pid.fungi=1-f_mismatched.fungi) %>%
  dplyr::rename(orf=id,n_snp = n_mismatched.yk11) %>%
  dplyr::select(-f_mismatched.fungi) %>%
  relocate(orf,UNIPROT,GENENAME,SGD, OG, HAS_ORTHOLOG, len_ref,
           len_msa.yk11, n_snp,f_snp, len_msa.fungi,pid.fungi) %>%
  dplyr::select(-paste0(er_fungi_worst,'.fungi'),-paste0(er_strains_worst,'.yk11'))

write_delim(evo_yeast,gzfile(here::here('output','evolution-yeast-protein.tsv.gz')),delim = '\t')
saveRDS(evo_yeast,here::here('data','evolution-yeast-protein.rds'))

