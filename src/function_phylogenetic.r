source("src/utils.r",local = T)
# Phylogenetic data ------------------------------------------------------------
# load.phylogenetic = function() {
#   load("/media/elusers/users/benjamin/A-PROJECTS/01_PhD/06-phd-final-report/data/EVO/SC-phylo-ali-prot_data.Rdata")
#   remove(mySC, stat.msa, mySC.phy)
#   return(SC)
#   #scinfo = readRDS('yeast-proteome-sgd_infos.rds')
# }

get.ygob.pair = function(ygob = load.ygob.ohnologs() ){
  return(ygob[, c('orf', 'dup.orf')])
}


load.ygob.ohnologs = function() {
  message("REF: Byrne and Wolfe, 2005, Genome Research")
  message("The Yeast Gene Order Browser: Combining curated homology and syntenic context reveals gene fate in polyploid species")
  # http://ygob.ucd.ie/
  ygob.ohno = read.csv("http://ygob.ucd.ie/browser/ohnologs.csv", # s. cerevisiae ohnologs
                       stringsAsFactors = F)[, -1]
  colnames(ygob.ohno) = c('WGD_anc', 'orf1', 'gname1', 'orf2', 'gname2', 'pid', 'rlen')
  regexORF = "^[Y][A-P][LR][0-9]{3}[WC](?:-[A-Z])?$"
  isORF1 = grepl(regexORF, ygob.ohno$orf1)
  isORF2 = grepl(regexORF, ygob.ohno$orf2)
  regexYGOB = "^(Scer_YGOB_)(Y\\w+|[^Y]\\w+)$"
  # Scer_YGOB_SDC25 = YLL016W (SDC25)
  # Scer_YGOB_YDR134C= YDR134C (CCW22)
  fake = data.frame(
    old.orf = c("Scer_YGOB_SDC25", "Scer_YGOB_YDR134C"),
    orf = c("YLL016W", "YDR134C"),
    gname = c("SDC25", "CCW22"),
    stringsAsFactors = F
  )
  ygob  = ygob.ohno %>%
    mutate(
      fake.orf = orf1 %in% fake$old.orf,
      gname1 = replace(gname1, fake.orf, fake$gname),
      orf1 = replace(orf1, fake.orf, fake$orf),
      pid = as.numeric(gsub("^(\\d+)\\%", "\\1", pid)),
      ref = 1
    ) %>%
    dplyr::rename(orf = "orf1", gname = "gname1", dup.orf = "orf2", dup.gname = "gname2" ) %>%
    dplyr::select(WGD_anc, orf, gname, ref, dup.orf, dup.gname, pid, rlen)
  return(ygob)
}

load.pombe.orthologs = function() {
  url.orthologs = "ftp://ftp.pombase.org/pombe/orthologs/cerevisiae-orthologs.txt"
  sp.sc = read.delim(url.orthologs, comment.char = "#", stringsAsFactors = F,
                     header = F, sep = '\t', row.names = NULL,
                     col.names = c('PombaseID','ORFS'))
  orthologs = sp.sc %>%
    separate_rows(ORFS,sep = "\\|") %>%
    mutate( PombaseID=str_trim(PombaseID),
            ORF =  str_trim(str_remove(ORFS,"\\((FUSION-)?(N|C)\\)")),
            splitted = str_extract(ORFS,"(?<=\\()(N|C)(?=\\))")) %>%
    dplyr::select(-ORFS) %>%
    # Remove pombe genes with no orthologs in cerevisiae (MEL genes don't exist in yeast)
    filter(ORF != "NONE" & !grepl("MEL[1256]",ORF) ) %>%
    group_by(PombaseID) %>% mutate(pombe.1 = n()==1 ) %>%
    group_by(ORF) %>% mutate(cerevisiae.1 = n()==1 ) %>%
    rowwise %>% mutate(ortho.1to1 = pombe.1 & cerevisiae.1)

  return(orthologs)
}

get.ortho.pair = function(ortho = load.pombe.orthologs() ){
  return(ortho[,c('PombaseID','ORF')])
}

get.ygob.pair = function(ygob = load.ygob.ohnologs() ){
  return(ygob[, c('orf', 'dup.orf')])
}


compare.to.ancestors = function(ancestor, current){
  ancestral.genome = read.delim(file = ancestor, header = T,sep = '\t',stringsAsFactors = F)
  # Extract yeast ancestor genes
  ancestral.genome$proxy_yeast = gsub(".+(YEAST.+)","\\1",ancestral.genome$proxy_genes)
  ancestral.genome$proxy_yeast[ !grepl("YEAST",ancestral.genome$proxy_genes) ] = NA
  ancestral.genome$proxy_yeast_sgd = str_replace(ancestral.genome$proxy_yeast, ".+SGD\\=(S[0-9]+).*",replacement = "\\1")
  ancestral.genome$proxy_yeast_uni = str_replace(ancestral.genome$proxy_yeast, ".+UniProtKB\\=([A-Z0-9]+).*",replacement = "\\1")

  tmp.prot = AAStringSet(gsub("\\-","",ancestral.genome$protein_sequence))
  names(tmp.prot) = ancestral.genome$proxy_yeast_sgd
  anc.prot = tmp.prot[!is.na(names(tmp.prot))]
  aafreq.anc = alphabetFrequency(tmp.prot,as.prob = T)
  aacount.anc = alphabetFrequency(tmp.prot,as.prob = F)

  ancestor.ali = align.pair.prot(p1=current[names(anc.prot)],p2=anc.prot)
  aafreq.anc.ali = alphabetFrequency(unaligned(subject(ancestor.ali)),as.prob=T)
  aacount.anc.ali = alphabetFrequency(unaligned(subject(ancestor.ali)),as.prob=F)

  anc.stat = data.frame(sgdid=names(anc.prot),
                        pid1.anc=pid(ancestor.ali,'PID1'),
                        pid2.anc=pid(ancestor.ali,'PID2'),
                        pid3.anc=pid(ancestor.ali,'PID3'),
                        pid4.anc=pid(ancestor.ali,'PID4'),
                        nX.anc = aacount.anc.ali[,'X'],
                        pX.anc = 100*aafreq.anc.ali[,'X'],
                        pG.anc = 100*rowSums(aafreq.anc.ali[,c('-','+')]),
                        L.anc = width(unaligned(subject(ancestor.ali))),
                        S.anc = nmatch( ancestor.ali ),
                        N.anc = nmismatch( ancestor.ali )
  )

  return(anc.stat)
}

get.enog.4891 = function(taxid_info='http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv',
                         node_info='./data/EVO/eggNOG-data/by_level_taxa/4891.proteomes'){
  sp.4891 = data.frame(
    taxid=c('4952','5476', # outgroup
            '4956','4950','28985','45285','33169','381046', # prewgd
            '4932','5478','1071379','113608','588726','27288','27289','36033','432096', # postwgd
            '1041607','42374','273371','36914','340170','45596','4909','4920','4922','4924','4929','4959','5479','5480','5482','1005962'  # nonwgd
    ),
    abb = c('Y. lipolytica','C. albicans', # outgroup
            'Z. rouxii','T. delbrueckii','K .lactis','E. cymbalariae','E. gossypii','L. thermotolerans', # prewgd
            'S. cerevisiae','C. glabrata','T. blattae','T. phaffii','K. naganishii','N. castellii','N. dairenensis','V. polyspora','K. africana', # postwgd
            "W. ciferrii", "C. dubliniensis", "C. orthopsilosis", "L. elongisporus", "S. passalidarum", "C. tenuis", # nonwgd
            "P. kudriavzevii", "M. farinosa", "K. pastoris", "S. stipitis", "M. guilliermondii", # nonwgd
            "D. hansenii", "C. maltosa", "C. parapsilosis", "C. tropicalis", "O. parapolymorpha" # nonwgd
    ),
    clade=c('out','out',
            'prewgd','prewgd','prewgd','prewgd','prewgd','prewgd',
            'postwgd','postwgd','postwgd','postwgd','postwgd','postwgd','postwgd','postwgd','postwgd',
            'nonwgd','nonwgd','nonwgd','nonwgd','nonwgd','nonwgd','nonwgd','nonwgd','nonwgd','nonwgd',
            'nonwgd','nonwgd','nonwgd','nonwgd','nonwgd','nonwgd'
    )
  )
  library(rotl)
  taxid = read.delim(taxid_info,header=T,stringsAsFactors = F)
  node = read.delim(node_info,header=F,stringsAsFactors = F, col.names = c('taxid','sciname','rank','ancestors','lineage'))
  LCA = Reduce(intersect,strsplit(node$ancestors,","))
  LCA.taxid = Reduce(intersect,strsplit(node$lineage,","))
  node$LCA = tail(LCA,1)
  node$LCA.taxid = tail(LCA.taxid,1)
  MRCA = find.common.ancestor(node$ancestors)
  node$MRCA = sapply(MRCA,b4.last)
  node$ancestors = toString(LCA)
  node$lineage = toString(LCA.taxid)
  node.info = merge(node,sp.4891, by='taxid')
  return(node.info)
}

find.common.ancestor= function(lineage){
  L = strsplit(lineage,',')
  LCA = Reduce(intersect, L)
  MRCA = sapply(L,function(x){ setdiff(unlist(but.last(x)), LCA) })
  return(MRCA)
}
