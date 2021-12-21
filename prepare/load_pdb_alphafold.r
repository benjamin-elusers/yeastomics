source("https://raw.githubusercontent.com/benjamin-elusers/yeastomics/main/src/__setup_yeastomics__.r")
library(bio3d)

##### 1. UNIPROT PROTEOME #####
scprot=load.dubreuil2019.data(4)
scres=load.dubreuil2019.data(3)
scevo=load.dubreuil2021.data(3)

saved.alphafold = here::here('prepare','alphafold-uniprot-yeast.rds')
id_uniprot = scprot$UNIPROT
if(!file.exists(saved.alphafold )){
  #stickiness = get.stickiness()[aa_af]
  #group_by(uni,plddt_bin) %>%
  #mutate(sti.afbin = mean(stickiness))
  afprot = get.alphafold.proteome(scprot$UNIPROT)

  # Average stickiness per region to wide format
  # afprot %>% pivot_wider(id_cols=c(uni,chain,resn,resi,plddt,aa_af), names_from=plddt_bin, names_glue = "sti.af_{plddt_bin}",values_from = sti.af )
  # CHECKPOINT
  write_rds(afprot,here::here('prepare','alphafold-uniprot-yeast.rds'))
}
afprot=read_rds(saved.alphafold )


get_alphafold_content = function(afprot, plddt_cutoff=30){
  #afprot=read_rds( here::here('prepare','alphafold-yeast-uniprotKB_dubreuil2019.rds') )
  # By proteins
  AF.content = afprot %>% mutate( plddt30 = plddt<=30 ) %>%
       group_by(uni) %>%
       summarise(fct_count(plddt_bin,prop=T)) %>%
       dplyr::rename(bin=f,f=p,L=n) %>%
       pivot_wider(id_cols=c(uni,bin),names_from=bin,names_glue="{.value}.af_{bin}",values_from=c(f,L))
  return(AF.content)
}

test = get_alphafold_content(afprot)

get_alphafold_stickiness = function(){
  afprot=read_rds( here::here('prepare','alphafold-yeast-uniprotKB_dubreuil2019.rds') )
  AF.sti = group_by(afprot, uni,plddt_bin,sti.afbin) %>%
           pivot_wider(id_cols=c(uni,plddt_bin),names_from=plddt_bin,names_glue="stickiness.af_{plddt_bin}",values_from = sti.afbin,values_fn = {unique})
}

get_d2p2_content = function(MIN_D2P2=7){

  # Get D2P2 predictions of disorder
  sc_d2p2 = read_rds(here::here('data','d2p2-yeast-uniprotKB.rds')) %>%
            get.d2p2.diso(.,as.df = T)
    #        mutate(d2p2.seg = find.consecutive(d2p2.diso>=MIN_D2P2, TRUE, min=3),
    #               d2p2.gap = find.consecutive(d2p2.diso>=MIN_D2P2, FALSE, min=1)) %>%
    #group_by(d2p2.seg) %>% mutate( d2p2.seglen = sum_(d2p2.seg!=0)) %>%
    #group_by(d2p2.gap) %>% mutate( d2p2.gaplen = sum_(d2p2.gap!=0))

  D2P2.content = sc_d2p2 %>% dplyr::filter(has.d2p2) %>%
    dplyr::select(-c(has.d2p2,d2p2.size)) %>%
    group_by(d2p2.id) %>%
    summarise(L.d2p2 = sum_(d2p2.diso>=MIN_D2P2),
              f.d2p2 = mean_(d2p2.diso>=MIN_D2P2))
              #d2p2.nseg = n_distinct(d2p2.seg) - 1*(sum_(d2p2.seg==0)>0),
              #d2p2.Lseg = sort(unique(d2p2.seglen[d2p2.seglen>0])),
              #d2p2.iseg = row_number(d2p2.Lseg),
              #d2p2.Lsegmax = max(d2p2.seglen))
  return(D2P2.content)
}

prot = scprot %>% dplyr::select(UNIPROT,STRING,GENENAME,
                         standard,medium,high,rejected,isMB,
                         mean.mpc,median.mpc,bins.abundance,ppm,
                         prot.size,starts_with(c('L.','f.','stickiness.')))

d2p2.content = get_d2p2_content(MIN_D2P2=7)
af.content = get_alphafold_content()

df= prot %>%
  left_join(d2p2.content,by=c('UNIPROT'='d2p2.id')) %>%
  left_join(af.content,by=c('UNIPROT'='uni'))

ggplot(df, aes(y=f.d2p2,x=f.af_D)) +
  geom_point() + facet_wrap(~bins.abundance)
cor.sub.by(df,'f.d2p2','f.af_D',BY = 'bins.abundance')
spearman.toplot(df$f.d2p2,df$f.af_D)


AF = left_join(get_alphafold_content(),get_alphafold_stickiness()) %>%
  # Reorder columns by type of bins of alphafold confidence scores
  dplyr::select(uni, str_order(str_sub(colnames(.),start=-4)) )
head(AF)




tmp1 = sc_d2p2 %>%
        group_by(d2p2.id,d2p2.size,d2p2.diso,is_d2p2=d2p2.diso!=0) %>%
        summarise(,fct_count(as.character(d2p2.diso),prop=T))

tmp1 %>% pivot_wider(id_cols=c(d2p2.id,d2p2.size,is_d2p2), names_from='d2p2.diso',values_from=n)


tmp2 = tmp1 %>% group_by(d2p2.id,d2p2.size,nod2p2=d2p2.diso==0) %>%
            summarise(cumsum(n))
            #pivot_wider(id_cols=c(d2p2.id,d2p2.size),names_from=f,names_glue="d2p2_{f}",values_from=n)
head(tmp)


# Get uniprot sequence
sc_uniseq = scres %>% select(id,prot.size,resnum,resname)
D2P2.sti = inner_join(sc_d2p2,sc_uniseq,c('d2p2.id'='id','d2p2.resi'='resnum')) %>%
            dplyr::rename(aa_d2p2=resname) %>%
            mutate(stickiness = get.stickiness()[aa_d2p2]) %>%
            group_by(d2p2.id,d2p2_cons=(d2p2.diso>=7)) %>%
            summarise(L.d2p2=sum_(d2p2.diso>=7), fct_count(as.character(d2p2.diso)),
                      stickiness.d2p2_raw = ifelse(d2p2_cons,mean_(stickiness),NA),
                      stickiness.d2p2 = ifelse(L.d2p2>=10,mean_(stickiness),NA))
head(D2P2.sti)


SC = left_join(scprot,AF,by=c('UNIPROT'='uni')) %>% select(-starts_with('Wsti')) %>%
     left_join(D2P2,c('UNIPROT'='d2p2.id'))
SC
write_rds(SC,file = here::here('prepare','sc_alphafold.rds'))
##
# CHECKPOINT
SC=read_rds(here::here('prepare','sc_alphafold.rds'))

library(ggplot2)

SC.f1 = SC %>% dplyr::filter(!is.na(bins.abundance) & !is.na(f)) %>%
        group_by(bins.abundance,f) %>% add_count(name='nprot') %>%
        group_by(f) %>% add_count(name='n') %>%
        dplyr::rename(alphafold_type = f, alphafold_frac = p) %>%
        left_join(SC.f1,sc_d2p2, c('UNIPROT'='d2p2.id')


levels(SC.f1$alphafold_type) = c('Very low (pLDDT < 50)','Low (70 > pLDDT > 50)','Confident (90 > pLDDT > 70)','Very high (pLDDT > 90)')

F1 = ggplot(SC.f1,aes(x=as.factor(bins.abundance)))
F1 = F1 + geom_boxplot(aes(fill=alphafold_type,y=alphafold_frac),varwidth=T,notch=T,outlier.shape=NA,outlier.size = 0, coef = 0)
F1 = F1 + facet_wrap(~alphafold_type,drop = T,scales = 'free_y',nrow = 4)
F1 = F1 + geom_text(mapping=aes(y=Inf,label=nprot),size=3,hjust=1,angle=90,vjust=1,check_overlap = T)
F1 = F1 + geom_text(mapping=aes(y=-Inf,x=Inf,label=n),check_overlap = T,vjust='inward',hjust='inward')
F1 = F1 + theme_minimal() + theme(legend.position = 'bottom', strip.text = element_blank()) + ggpubr::grids()
F1

F2 = ggplot(SC.f1,aes(x=bins.abundance))
F2 = F2 + geom_boxplot(aes(fill=as.factor(bins.abundance),x=as.factor(bins.abundance),y=f.IUP20),varwidth=T,notch=T,outlier.shape=NA,outlier.size = 0, coef = 0)
F2 = F2 + theme_minimal() + theme(legend.position = 'bottom', strip.text = element_blank()) + ggpubr::grids()
F2 = F2 + scale_fill_grey(start=0.7,end = 0.3) + ylim(0,50)
F2
library(cowplot)
save_plot(plot =F1,filename = "~/Desktop/alphafold-vs-abundance.png",base_aspect_ratio = 0.8,base_height = 10)
save_plot(plot =F2,filename = "~/Desktop/iup20-vs-abundance.png",scale=1)
F3= ggplot(SC.f1 %>% filter(alphafold_type=='Very low (plDDT < 50)'),aes(x=as.factor(bins.abundance)))
F3 = F3 + geom_boxplot(aes(fill=alphafold_type,y=alphafold_frac),varwidth=T,notch=T,outlier.shape=NA,outlier.size = 0, coef = 0)
F3 = F3 + geom_text(mapping=aes(y=Inf,label=nprot),size=3,hjust=1,angle=90,vjust=1,check_overlap = T)
F3 = F3 + geom_text(mapping=aes(y=-Inf,x=Inf,label=n),check_overlap = T,vjust='inward',hjust='inward')
F3 = F3 + theme_minimal() + theme(legend.position = 'bottom', strip.text = element_blank()) + ggpubr::grids()
F3

save_plot(plot =F3,filename = "~/Desktop/alphafold-lowconf-vs-abundance.png",base_aspect_ratio = 0.8,base_height = 10)



spearman(100*SC.f1$alphafold_frac[SC.f1$alphafold_type=='D'],SC.f1$f.domain[SC.f1$alphafold_type=='D'])

cor.sub.by(SC.f1, 'sti.avg', 'ppm', c('isMB','alphafold_type'))
#cor.sub.by(SC.f1, 'sti.avg', 'ppm', c('f','high'))
