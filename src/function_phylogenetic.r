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

read.R4S = function(r4s, id=NULL,verbose=T){
  require(tidyverse)
  if(is.null(id)){ id = basename(r4s) }
  if(verbose){ message(sprintf('reading r4s results %s\n',basename(r4s))) }
  #Rates were calculated using the expectation of the posterior rate distribution
  #Prior distribution is Gamma with 16 discrete categories

  #SEQ: the amino acid in the reference sequence in one letter code.
  #SCORE: The conservation scores. lower value = higher conservation.
  #QQ-INTERVAL: the confidence interval for the rate estimates. The default interval is 25-75 percentiles
  #STD: the standard deviation of the posterior rate distribution.
  #MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.

  #POS SEQ  SCORE    QQ-INTERVAL     STD      MSA DATA
  #The alpha parameter 0.05
  #The likelihood of the data given alpha and the tree is:
  #LL=-9966.53
  #1     M 0.0002242   [9.65e-19,3.989e-06] 0.001014 1011/1011
  #2     V  0.0168   [0.002948,0.01779] 0.003824 1011/1011
  #3     L 0.0004615   [9.65e-19,4.613e-05] 0.002097 1011/1011
  #4     T 0.006544   [0.0004101,0.01779] 0.006853 1011/1011
  #5     I 0.0002452   [9.65e-19,3.989e-06] 0.001131 1011/1011
  r4s.col = c('POS','SEQ','SCORE','QQ_INTERVAL','STD','MSA')

  # Make sure QQ-INTERVAL does not have space after the comma
  clean_r4s = readLines(r4s) %>% gsub("(,\\s)",",",x = .)

  df.r4s = readr::read_table2(file = clean_r4s, comment = '#', col_names = r4s.col) %>%
           mutate(
             ID = id,
             QQ = str_remove_all(QQ_INTERVAL, pattern = "\\[|\\]"),
             QQ1 = as.double(str_split_fixed(QQ,',',n=2)[,1]),
             QQ2 = as.double(str_split_fixed(QQ,',',n=2)[,2])
           ) %>% select(ID,POS,SEQ,SCORE,QQ1,QQ2,STD,MSA)

  return(df.r4s)
}

read.R4S.param = function(r4s, as.df=F){
  require(tidyverse)
  r4s.param = readLines(r4s) %>% grep(pattern="^#", x = ., value = T)
  if(as.df){
    df.param = tibble(
      id = basename(r4s),
      prior =  r4s.param %>% grep(x=.,"Prior",v=T) %>% gsub(pat="#Prior distribution is (\\w)",repl="\\1",x=.),
      alpha = r4s.param %>% grep(x=.,"alpha parameter",v=T) %>% gsub(pat="#The alpha parameter ",repl=""),
      LL =  r4s.param %>% grep(x=.,"LL",v=T) %>% gsub(pat="#LL=",repl=""),
      r4s.avg =  r4s.param %>% grep(x=.,"Average",v=T) %>% gsub(pat="#Average = ",repl=""),
      r4s.std=  r4s.param %>% grep(x=.,"Deviation",v=T) %>% gsub(pat="#Standard Deviation = ",repl="")
    )
    return(df.param)
  }
  return(r4s.param)
}
