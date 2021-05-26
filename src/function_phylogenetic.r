#source("src/utils.r",local = T)
#source("src/function_alignment.r",local = T)

# Phylogenetic data ------------------------------------------------------------
# load.phylogenetic = function() {
#   load("/media/elusers/users/benjamin/A-PROJECTS/01_PhD/06-phd-final-report/data/EVO/SC-phylo-ali-prot_data.Rdata")
#   remove(mySC, stat.msa, mySC.phy)
#   return(SC)
#   #scinfo = readRDS('yeast-proteome-sgd_infos.rds')
# }

get.pair.prot = function(prot, pair){
  seq.1 = prot[pair[[1]]]
  seq.2 = prot[pair[[2]]]
  return(list(s1=seq.1,s2=seq.2))
}

pairwise.alignment.to.df  = function(ali){
  library(dplyr)
  N=length(ali)
  ali.list=list()
  doing="Global alignnment of pairwise sequences..."
  message(doing)
  tic(doing)
  for( i in 1:N ){
    prg = sprintf("%.1f %% [%5s/%5s]        \r", 100*(i/N), i,N)
    cat(prg)
    A = ali[i]
    P=alignedPattern(A)
    S=alignedSubject(A)
    p=unaligned(P)
    s=unaligned(S)

    S1=unlist(strsplit(x = toString(P),''))
    S2=unlist(strsplit(x = toString(S),''))
    resi = seq_along(S2)
    m.ali = tibble( pos=resi, aligned = (S1==S2), nid = nmatch(A),
                    s1 = names(P), s2 = names(S),
                    aa.s1 = S1, gap.s1 = (aa.s1=='-'),
                    aa.s2= S2, gap.s2 = (aa.s2=='-'),
                    l.s1 = width(p), l.s2 = width(s), l=width(P),
                    ol.12 = l.s1/l.s2, ol.21 = l.s2/l.s1,
                    pid.aligned = pid(A,type = "PID2" ),
                    pid.short = pid(A,type = "PID3" ),
                    pid.long  = 100* (nid / max(l.s1,l.s2)) ) %>%
      group_by(gap.s2) %>% mutate( pos.s2 = ifelse(gap.s2, NA, no=row_number() ) ) %>%
      group_by(gap.s1) %>% mutate( pos.s1 = ifelse(gap.s1, NA, no=row_number() ) ) %>%
      ungroup() %>% distinct()
    ali.list[[i]] = m.ali
  }
  toc()
  return( bind_rows(ali.list) )
}

align.pair.prot  = function(p1,p2,mat = "BLOSUM62",tomatrix=F, opening=10, extend=4 ){
  if(length(p1) != length(p2))
    stop("p1 and p2 should have the same length!")
  ali = pairwiseAlignment(p1, p2, substitutionMatrix = mat, gapOpening=opening, gapExtension=extend)
  if(tomatrix){
    return( pairwise.alignment.to.df(ali) )
  }else{
    return(ali)
  }
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

get.ygob.pair = function(ygob = load.ygob.ohnologs() ){
  return(ygob[, c('orf', 'dup.orf')])
}

get.sc.ohno = function(myseq) {
  ygob = load.byrne2005.data()
  bogy = ygob %>%
    dplyr::rename( orf = dup.orf, gname = dup.gname, dup.orf = orf, dup.gname = gname ) %>%
    mutate(ref = 2) %>%
    dplyr::select(WGD_anc, orf, gname, ref, dup.orf, dup.gname, pid, rlen)
  ohno = ygob %>% bind_rows(bogy)

  if(missing(myseq)){
    warn("No input proteome sequences! Cannot compare ohnolog sequences...")
    return(ohno)
  }
  ohno.seq = get.pair.prot(prot=myseq, pair = get.ygob.pair(ohno))
  ohno.ali = align.pair.prot(p1=ohno.seq$s1, p2=ohno.seq$s2, mat='BLOSUM62')
  aafreq = alphabetFrequency(alignedSubject(ohno.ali))
  gaps = rowSums(aafreq[,c("-","+")])

  library(org.Sc.sgd.db)
  pair = get.ygob.pair(ohno) %>%
    mutate( L1 = width(ohno.seq$s1),
            L2 = width(ohno.seq$s2),
            RLEN = round(pmin(L1,L2) / pmax(L1,L2), 2),
            PID1 = round(pid( ohno.ali, "PID1" ),1),
            PID2 = round(pid( ohno.ali, "PID2" ),1),
            PID3 = round(pid( ohno.ali, "PID3" ),1),
            PID4 = round(pid( ohno.ali, "PID4" ),1),
            SCORE.B100 = score( ohno.ali ),
            S = nmatch( ohno.ali ),
            N = nmismatch( ohno.ali ),
            G = gaps,
            SGDID = AnnotationDbi::select(org.Sc.sgd.db,keys = ohno$orf,columns ="SGD",keytype = 'ORF')[,2],
            SGDID.dup = AnnotationDbi::select(org.Sc.sgd.db,keys = ohno$dup.orf,columns ="SGD",keytype = 'ORF')[,2]
    ) %>% right_join(ohno)

  return(pair)
}

get.ortho.pair = function(ortho = load.pombe.orthologs() ){
  return(ortho[,c('PombaseID','ORF')])
}

compare.to.ancestors = function(ancestor, current){
  # ancestor is the data table containing ancestral sequences (AncientGenomes.org)
  # current is a Biostrings object of proteome sequences of a particular species (S. cerevisiae for example)
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

load.eggnog.node = function(node=NULL,to.matrix=F,show.nodes=F){
  eggnog = lst(baseurl="http://eggnog.embl.de/download/",
                v4.5=paste0(baseurl,"eggnog_4.5/"),
                latest=paste0(baseurl,"latest/"),
                tax_level = paste0(v4.5,"eggnog4.taxonomic_levels.txt"),
                tax_info =paste0(latest,"e5.taxid_info.tsv"),
                level_info =  paste0(latest,"e5.level_info.tar.gz")
               )

  eggnog_nodes=readr::read_delim(eggnog$tax_level,delim="\t",col_types = cols(.default="c")) %>%
    janitor::clean_names(replace = c("#"="")) %>%
    mutate(node_full = sprintf(" %-s (%s)", paste(tax_id,tolower(level_name),sep="."), nog_prefix)) %>%
    dplyr::select(-ends_with("_count")) %>% arrange(tax_id)

  eggnog_taxons=readr::read_delim(eggnog$tax_info,delim="\t", col_types = "ccccc") %>%
                janitor::clean_names(replace = c("#"=""))

  # Testing
  # node.test = list(null=NULL,empty=c(),none="",num=4751,name='Fungi', nog='fuNOG')
  # node = node.test$nog
  nodes=dplyr::select(eggnog_nodes,tax_id:level_name)
  if(show.nodes){ return(eggnog_nodes) }
  node_exists= !purrr::is_empty(node)
  is_eggnog_node=F

  if(node_exists){
    node_pos  = which( node == nodes, arr.ind=T) # check if the node exists
    if( nrow(node_pos) == 1 ){
      inod=node_pos[,'row']
      node_type = names(nodes)[ node_pos[,'col'] ] # find what column the node is found (tax_id, nog_prefix or level_name)
      message("(",node,") is a valid eggnog node of type [",node_type,"]")
      is_eggnog_node=T
    }else if( nrow(node_pos)>1){
      warning("(",node,") is not unique! Please select a valid unique node",immediate.=T)
      is_eggnog_node=F
    }
  }
  nodes_list = eggnog_nodes$node_full %>%  gtools::mixedsort()

  if( !is_eggnog_node ){
    warning(node," is not a valid eggnog taxonomic level!",immediate.=T)
    inod = menu( choices=nodes_list )
  }
  node_info = eggnog_nodes[inod,] %>% dplyr::select(-node_full)
  print(node_info)

  node_members = sprintf("%s/per_tax_level/%s/%s_members.tsv.gz",eggnog$latest,node_info$tax_id,node_info$tax_id)
  df.ortho = readr::read_delim(node_members,"\t",progress=T,
                    col_names = c('node','OG','ns','np','ortho','taxons'),
                    col_types=cols(.default="c")) %>%
    mutate(has_yeast = str_detect(string = taxons,pattern = "4932")) %>%
    separate_rows(ortho,sep=",") %>%
    extract(ortho,into=c('taxid','protid'),regex='^(^[0-9]+)\\.(.+)') %>%
    group_by(OG,taxid) %>% mutate(is_1to1=n()==1) %>%
    left_join(node_info, by=c('node'='tax_id'))

  if( to.matrix ){
    ortho.m= df.ortho %>%
              filter(is_1to1) %>%
              pivot_wider(id_cols=c(OG,np),names_from=taxid, values_from=protid)
    return(ortho.m)
  }
  return(df.ortho)
}

#NODES = get.eggnog.node(node = 4751,to.matrix = F)
#colnames(NODES)

get.enog.4891 = function(){
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
  library(tidyverse)
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

  # Make sure QQ-INTERVAL does not have any space in between brackets
  # LOOKBEHIND AND LOOKAHEAD library perl regex engine
  # gsub(x = "   34     I 0.005863   [0.0001698, 0.004] 0.00562 1011/1011",  pattern = "(?<=\\[)([^\\]]*)( +)","\\1",perl = T)
  # gsub(x = "  131     L    2.02   [0.4701,  2.02]       0 1011/1011",  pattern = "(?<=\\[)([^\\]]*)(\\s\\s?)([^\\s\\]]*)(?=\\])","\\1\\3",perl = T)
  # test="  131     L    2.02   [0.4701,  2.02]       0 1011/1011"
  # test="   34     I 0.005863   [0.0001698, 0.004] 0.00562 1011/1011"

  clean_r4s =readLines(r4s) %>% str_replace_all(string = ., pattern = "\\s+(?=[^\\[\\]]*\\])", replacement="")
  #gsub(x = .,  pattern = "(?<=\\[)([^\\]]*)( +)([^\\]]*)(?=\\])","\\1\\3",perl = T)
    #gsub("(?<=\\[)(\\s+)","",x = .,perl = T) # Remove spaces after bracket
    #gsub("(\\s+)(?=\\])","",x = .,perl = T) # Remove spaces before bracket
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
  library(tidyverse)
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

