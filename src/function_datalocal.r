source("src/utils.r",local = T)
source("src/function_sequence.r",local = T)
source("src/function_phylogenetic.r",local = T)

read.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con)
  return(txt)
}

subname=function(name,sep="\\.",lc=F){ # extracts substring until first separator
  b4sep = sprintf("([^%s]+).+",sep)
  part1 = sub(b4sep, "\\1", x=name)
  if(lc){ tolower(part1) }
  return(part1)
}

strfind = function(strings, patterns){ # search multiple patterns in character vectors
  sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = T) })
}

# Local proteome data ----------------------------------------------------------
load.emmanuel.data = function(toolbox="/data/elevy/70_R_Data/bin/RToolBox_yeast_general.R"){
  source(toolbox,local = T)
  library(AnnotationDbi)
  library(org.Sc.sgd.db)
  library(GO.db)
  library(tidyverse)

  SC = get.proteome.table() %>%
    mutate(
      has_len = !is.na(len),
      rel_diso1 = diso1/len, rel_diso2 = diso2/len, diso05 = diso05/len,
      has_diso05 = !is.na(diso05), has_diso1 = !is.na(diso1), has_diso2 = !is.na(diso2),
      has_tox = !is.na(over.tox),
      has_loc_ymd = !is.na(loc.ymd), has_loc_dtt = !is.na(loc.dtt),
      has_loc_h2o2 = !is.na(loc.h2o2), has_loc_starv = !is.na(loc.starv),
      has_viability = !is.na(viable)
    )
  return(SC)
}

load.1011.strains= function(seqdir="/media/elusers/users/benjamin/A-PROJECTS/01_PhD/02-abundance-evolution/strains1011/data/sequences/Proteome_1011/",
                             .recursive=F){
  if( !dir.exists(seqdir) ) stop("Directory of proteome sequences not found!")
  fastas =list.files(path=seqdir, pattern = "fasta", full.names = T, ignore.case = T, include.dirs = F,recursive = .recursive)
  return(read.proteomes(fastas,strip.fname=T))
}

# Evolution Sequence/Structure -------------------------------------------------
# Precomputed data #
#==================#
# PROTEOME=readRDS('./data/PROTEIN-EVORATE.rds')
# ALIGNED.DATA = readRDS('./data/RESIDUE-EVORATE.rds')
#==================#

load.wapinsky2007.data = function(path.data="./data/"){
  library(tictoc)
  doing = "Get rate4site data for yeast based on wapinsky 2007 fungi lineage"
  tic(doing)
  #rate4site-yeast-fungi_lineage.tsv.gz
  r4s.dataset=sprintf("%s/rate4site-yeast-fungi_lineage.tsv.gz", path.data)
  r4s.orf = read.csv(file = r4s.dataset, sep='\t', skip=0, header=T, fill=T,
                     strip.white=T, stringsAsFactors = F,
                     col.names = c('r4s_orf','r4s_size','r4s_pos','r4s_aa','r4s_res',
                                   'r4s_std','r4s_nmsa','species','r4s_qq1','r4s_qq2',
                                   'R4S_num','R4S_aa','R4S_res','R4S_qq1','R4S_qq2','R4S_std','R4S_nmsa'))
  toc()
  return(r4s.orf)
}

path.r4s = "/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/rate4site-3.2.0"
load.rate4site_1011.data = function(path.res = paste0(path.r4s,"/src/rate4site/RUN-1011G"),
                                    only_results=T){
  library(tictoc)
  library(purrr)
  library(progress)
  library(tidyverse)

  doing = "Get rate4site data for 1011 isolated yeast strains"
  message(doing)
  tic(doing)

  resdir = list.dirs(path.res,full.names = T,recursive = F)
  orfs = get.orf.filename(resdir)
  names(resdir)=orfs

  #seqin = read.proteomes(paste0(res, "/seqin"))
  #names(seqin)=orfs
  #ws = sapply(widths(seqin),unique)
  #r4s = tibble( orf = orfs, len=ws, strains=lengths(seqin) )

  S288C  = load.sgd.proteome(withORF = T,rm.stop = F)
  wr= setNames(widths(S288C),names(S288C))
  sgd = tibble( orf= names(S288C), len.s288c=wr )

  r4s.seq = sgd %>% #left_join(r4s,sgd, by=c('orf'='orf'), suffix=c('.1k11','.s288c')) %>%
            mutate( r4s.file = sprintf("%s/%s_raw.r4s",resdir[orf],orf),
                    found=file.exists(r4s.file),
                    size=file.size(r4s.file))

  get.r4s.raw = function(file,name,.verb=F,.pb=NULL){
    if(!.pb$finished){ .pb$tick() }
    r4s.raw = NA
    if(file.exists(file)){ r4s.raw = read.R4S(r4s=file, id=name,verbose = .verb) }
    return(r4s.raw)
  }

  pb =  progress::progress_bar$new(total = nrow(r4s.seq), width = 70,
                         format = " (:spin) reading r4s [:bar] :percent (elapsed: :elapsed # eta: :eta)")

  r4s.data = r4s.seq %>%
    rowwise() %>%
    mutate( r4s.res = pmap( list(file=r4s.file, name=orf),
                            get.r4s.raw, .pb = pb, .verb=F)
    )

  r4s.df = r4s.data %>% unnest(cols = r4s.res) %>% select(-r4s.res)
  toc(log=T)
  if(only_results){ res = r4s.df %>% select(-c(r4s.file,found,size) ) }
  return(res)
}

load.aligned.data = function(data.path='../data/'){
  message("
          Merged dataset with residue-level informations from:
          (1) SGD S288C proteome sequence
          (2) Uniprot reference yeast proteome (aligned to SGD)
          (3) Dubreuil et al. (2019) for Abundance, Disorder (IUP+D2P2) & stickiness (+aaindex)
          (4) Rate4Site based on fungi lineage Wapinsky et al. (2007) of 14 yeast species
          (5) 3Dcomplex for yeast quaternary structures (based on X-Ray from PDB)"
  )
  dataset = file.path(data.path, "sgduni-yeast-aligned-datasets.rds")
  if( !file.exists(dataset) ){ stop(sprintf("Cannot find the dataset at : %s",dataset)) }
  return(readRDS(dataset))
}

# Evolutionary rate per regions on quaternary structures
get.evo3d =  function(path.data='../data/'){
  ALIGNED.RES = get.aligned.data(path.data) %>%
    mutate(iup=IUP20!=0, noiup = (IUP20==0),
           d2p2=d2p2_diso>=7, nod2p2= d2p2_diso<3,
           dom = domain!=0 & noiup, nodom = antidomain!=0) %>%
    dplyr::select( ends_with('.id'),
                   starts_with('l_'),starts_with('pos_'),
                   starts_with('aa_'), starts_with('gap_'),
                   starts_with('has_'),
                   starts_with('R4S_',ignore.case = F),
                   ends_with('.3d'),
                   c(iup,noiup,d2p2,nod2p2,dom,nodom,l_evo, l_pdb, gap_pdb),
                   -c(code.id, chain_name.id, org_ref.id,accession.id,uniprotAcc.id,
                      R4S_nmsa, R4S_orf, R4S_orf_sd, R4S_uni, R4S_uni_sd, R4S_std,
                      aa_sgd, gap_sgd,
                      aa_uniprot, gap_uniprot,
                      pos_dubreuil, aa_dubreuil, l_dubreuil, gap_dubreuil,
                      noseq.3d, npdb.3d,
                      homo.3d, best_BU.3d,hydro_kyte.3d, stickiness_ec.3d ) ) %>%
    filter(has_R4S)
  EVO3D= ALIGNED.RES %>% filter(has_R4S & has_pdb)
  return(EVO3D)
}

# Evolutionary rate per regions on quaternary structure with ASA cutoff
get.evo3d.byasa = function(asa,aligned.data=NULL){
  message(sprintf('Surface is ASA above %s%%.\n',asa))
  if(is.null(aligned.data)){  aligned.data = get.evo3d() }
  byasa = aligned.data %>%
    group_by(UNIPROTKB.id,SGD.id,ORF.id,code.chain.id) %>%
    mutate( CUTOFF_LEN=20, CUTOFF_ASA = asa,
            L_evo = sum_(has_R4S), L_pdb= sum_(has_pdb),
            ASA_under = asa_rel_in_BU.3d <= 25 & asa_rel_alone.3d <=25,  ASA_over = asa_rel_in_BU.3d > CUTOFF_ASA,
            L_buried=sum_(ASA_under), L_surface=sum_(ASA_over),
            L_iup = sum_(iup),  L_d2p2 = sum_(d2p2),
            L_noiup = sum_(noiup),  L_nod2p2 = sum_(nod2p2),
            L_dom = sum_(dom), L_nodom = sum_(nodom),
            R_full = ifelse(L_evo>=CUTOFF_LEN, mean_(R4S_norm), NA),
            R_pdb = ifelse(L_pdb>=CUTOFF_LEN, mean_(R4S_norm[!gap_pdb]), NA),
            R_surface   = ifelse(L_surface>=CUTOFF_LEN,   mean_(R4S_norm[ASA_over]),   NA),
            R_buried    = ifelse(L_buried>=CUTOFF_LEN,    mean_(R4S_norm[ASA_under]),    NA),
            R_iup   = ifelse(L_iup>=CUTOFF_LEN,   mean_(R4S_norm[iup]),   NA),
            R_d2p2   = ifelse(L_d2p2>=CUTOFF_LEN,   mean_(R4S_norm[d2p2]),   NA),
            R_noiup   = ifelse(L_noiup>=CUTOFF_LEN,   mean_(R4S_norm[noiup]),   NA),
            R_nod2p2   = ifelse(L_nod2p2>=CUTOFF_LEN,   mean_(R4S_norm[nod2p2]),   NA),
            R_dom   = ifelse(L_dom>=CUTOFF_LEN,   mean_(R4S_norm[dom]),   NA),
            R_nodom   = ifelse(L_nodom>=CUTOFF_LEN,   mean_(R4S_norm[nodom]),   NA)
    )    %>%
    dplyr::select( UNIPROTKB.id, SGD.id, code.chain.id, ORF.id,
                   has_iup, has_pdb, has_R4S,
                   starts_with(c('R_','L_','CUTOFF_'))
    ) %>%
    filter(!is.na(UNIPROTKB.id)) %>%
    ungroup() %>% distinct() %>%
    rename_with(function(x){ evo.cols[x]},.cols=names(evo.cols))
}

# Disorder predictions ---------------------------------------------------------
fetch.d2p2 = function(id,quiet=F){ # Get the d2p2 predictions for single id
  library(rjson)
  # id='YDR134C'
  d2p2.url = sprintf('http://d2p2.pro/api/seqid/["%s"]',id)
  res = rjson::fromJSON(readLines(d2p2.url, warn=FALSE))
  nodata = length(res[[id]]) == 0
  if( !nodata ){
    # response$YCL051W[[1]][[3]]$disorder$consensus
    d = res[[id]]
    sp = d[[1]][[2]]
    if(!quiet){ message(sprintf("= %s from %s =\n",id,sp)) }
    consensus = d[[1]][[3]]$disorder$consensus
    record = list(id=id,from=sp,diso=d[[1]][[3]]$disorder, structure=d[[1]][[3]]$structure)
    return(record)
  }else{
    message(sprintf("id='%s' not recognized!\n",id))
    return(list(id=id,from=NULL,diso=NULL,structure=NULL))
  }
}

load.d2p2 = function(ids,saved){ # Get the d2p2 predictions for multiple ids
  if( file.exists(saved) ){
    if( file_ext(saved) != 'rds' ){ warning("File type not recognized ! (should be RDS object)") }
    return( readRDS(saved) )
  }
  if(!is.vector(ids)){ stop('Input must be a vector of identifiers as character...') }

  library(tools)
  library(tictoc)

  N=length(ids)
  doing =sprintf("Loading d2p2 predictions for %s protein identifiers\n",N)
  tic(doing)
  message(doing)
  d2p2 = list()
  for( i in 1:N ){
    ID = ids[i]
    prg = sprintf("%.1f %% [%5s/%5s]        \r", 100*(i/N), i,N)
    d2p2[[ID]] = fetch.d2p2(ID,quiet=T)
    cat(prg)
  }
  toc()
  message(sprintf("Saving D2P2 predictions to : %s\n",saved))
  saveRDS(d2p2, saved)
  return(d2p2)
}

get.d2p2.pred = function(d2p2){ return(d2p2$diso$consensus) }
get.d2p2.id = function(d2p2){ return(d2p2$id) }
get.d2p2.diso = function(d2p2,as.df=F){
  diso.L = sapply(d2p2, function(X){ setNames( object=list(get.d2p2.pred(X)) , nm=get.d2p2.id(X) ) })
  if( as.df ){
    resi.L = lapply(diso.L,function(X){ seq_along(X) })
    id.col = rep(names(diso.L),times=lengths(diso.L))
    len.col = rep(lengths(diso.L),times=lengths(diso.L))
    names(len.col)=NULL
    diso.df = tibble( d2p2_id=id.col,
                      d2p2_resi=unlist(resi.L,use.names = F),
                      d2p2_diso=unlist(diso.L,use.names = F),
                      d2p2_size=len.col,
                      has_d2p2=T)
    return(diso.df)
  }
  return(diso.L)
}

# Gene expression --------------------------------------------------------------
get.codons4tai = function(){
  codon.ord = c('TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG',
                'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG',
                'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG',
                'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',
                'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG',
                'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG',
                'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG',
                'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'
  )
  return(codon.ord)
}

load.trna.table = function(trna.tab="/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/codonR/tRNA/",
                            species="4932"){
    library(tidyverse)
    if( !dir.exists(trna.tab) ){ stop("Can't find the directory containing tables of trna counts!") }
    trna.filenames = list.files(path = trna.tab, pattern = '\\.trna')
    trna.files = file.path(trna.tab,trna.filenames)
    available.trna.tables = c(
      "1064592" = "Naumovozyma_castellii_CBS_4309",
      "1071378" = "Naumovozyma_dairenensis_CBS_421",
      "283643"  = "Cryptococcus_neoformans_var_neoformans_B-3501A",
      "284590"  = "Kluyveromyces_lactis_NRRL_Y-1140",
      "284593"  = "Candida_glabrata_CBS_138",
      "322104"  = "Scheffersomyces_stipitis_CBS_6054",
      "4896"    = "Schizosaccharomyces_pombe",
      "4932"    = "Saccharomyces_cerevisiae",
      "4950"    = "Torulaspora_delbrueckii",
      "4956"    = "Zygosaccharomyces_rouxii",
      "573826"  = "Candida_dubliniensis_CD36",
      "9606"    = "Homo_sapiens"
    )

    sp = match.arg(species, choices = names(available.trna.tables),several.ok = F)
    message(sprintf("reading tRNA counts for: %s\n",available.trna.tables[sp]))
    # # Name of the file = Taxon ID + initials of genus/species + ".trna"
    taxid = subname(trna.files,"_")
    # # First row starts with # followed by the name of the species
    spnames = sapply(trna.files,read.delim,sep='\t',header=F,nrows=1) %>%
      unlist %>% as.character %>%
      str_sub(start=2) %>% str_replace_all(" ","_")

    trna.counts = lapply(setNames(trna.files,taxid),
                         read.delim,sep='\t',skip=1,header=F,blank.lines.skip=T )
    mytrna = trna.counts %>% pluck(sp)

    return(mytrna)
}


load.trna.adaptation = function(inputseq,
                            trnas="data/GtRNA-counts.rds",
                            sp='4932',
                            exclude.counts=T
){
  message("REF: M. Dos Reis et al., 2004, Nucleic Acids Research")
  message("Solving the riddle of codon usage preferences: a test for translational selection")
  # https://doi.org/10.1093/nar/gkh834
  url.codonR="https://raw.githubusercontent.com/mariodosreis/tai/master/R/tAI.R"
  ord.codonR = get.codons4tai()

  if( !require(tAI) ){ devtools::install_github("mariodosreis/tai") } # Install this package first
  library(tAI)
  tai.R = sprintf("%s/tAI.R",url.codonR)
  # Can use local address
  #if(file.exists(codonR)){ tai.R = sprintf("%s/tAI.R",codonR) }
  if( ! class(inputseq) == "DNAStringSet" ){
    if(is.list(inputseq)){ dna = list2str(inputseq) }else{
      stop("Input sequences must be a list of character vectors or a Biostrings 'DNAStringSet'!")
    }
  }
  dna=inputseq
  codons = Biostrings::trinucleotideFrequency(dna,step = 3)[,ord.codonR]
  # IMPORTANT! Remove Methionine and STOP codons
  codons.out = which(colnames(codons) %in% c('ATG','TCA','TAA','TAG'))
  codons.60 = codons[,-codons.out]

  # Genomic tRNA database (copy/paste the species-specific tRNA gene counts to a text file)
  trna.count = readRDS(trnas)[[sp]]
  colnames(trna.count)=c('AA3','tot','cod','count','acod')

  trna.ord = match( ord.codonR, trna.count$acod )
  trna = trna.count[trna.ord,]

  ws =  tAI::get.ws(tRNA = trna$count, sking = 0)
  tai = tAI::get.tai(codons.60, ws)
  res = data.frame(prot=names(dna),tAI=tai)
  if(!exclude.counts){ res = cbind(res,codons) }
  return(res)
}

load.codon.usage= function(cds){
  if( !require(coRdon) ){ BiocManager::install("coRdon") } # Install this package first
  library(coRdon)
  cT=codonTable(cds)
  orf2ko = sub("ko:","",keggLink("ko",'sce'))
  names(orf2ko) =  sub("^.+:","",names(orf2ko))

  # Add KO identifiers to select ribosome
  cT=setKO(cT,ann=orf2ko[cT@ID])

  #cT.m = Biostrings::trinucleotideFrequency(CDS,step = 3)
  CU  = bind_cols(ID=cT@ID,
                  CU_milc=MILC(cT,ribo=T,self = F)[,1],
                  CU_enc=ENC(cT),
                  CU_b=B(cT,ribo=T,self = F)[,1],
                  CU_mcb=MCB(cT,ribo=T,self = F)[,1],
                  CU_encprime=ENCprime(cT,ribo=T,self = F)[,1],
                  CU_scuo=SCUO(cT),
                  CU_melp=MELP(cT,ribo=T)[,1],
                  CU_e=E(cT,ribo=T)[,1],
                  CU_cai=CAI(cT,ribo=T)[,1],
                  CU_gcb=GCB(cT,ribo=T),
                  CU_fop=Fop(cT,ribo=T)[,1],
                  CU_count= cT@counts %>% as_tibble(),
  )
  return(CU)
}


# Amino Acid features ----------------------------------------------------------
get.aascales=function(){
  data.frame(
    AA=get.AA1(),
    aggrescan=get.aggrescan(),
    camsol=get.camsol(),
    foldamyloid=get.foldamyloid(),
    kytedoolittle=get.kytedoolittle(),
    pawar=get.pawar_ph7(),
    roseman=get.roseman(),
    stickiness=get.stickiness(),
    voronoi_sickiness=get.voronoi_stickiness(),
    wimleywhite=get.wimleywhite()
  )
}


# Protein-Protein interactions (intact dataset from Hugo) ------------------------------
load.intact.yeast = function(only.direct=T,
                             rm.useless=T,
                             orga="cerevisiae",
                             intact.data=sprintf("/media/elusers/users/hugo/07_3DComplex_scripts/PPIs_analysis/INTACT/PPIs_%s.txt",orga),
                             from.MACOS=F){

  #/media/elusers/users/hugo/07_3DComplex_scripts/scripts_PPIs_networks/stack_studies_PPIs_nored_PMID_28122020_INTACT.R
  # intact.url="ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt"

  # search.keyword <- function(x, pos, keyword='cerevisiae|:4932|:559292|yeast') {
  #   colA = grepl(keyword,x[['Taxid.interactor.A']])
  #   colB = grepl(keyword,x[['Taxid.interactor.B']])
  #   colH = grepl(keyword,x[['Host.organism.s.']])
  #   x[colA & colB & colH,]
  #
  # }
  # test =readr::read_delim_chunked(file=intact.url, delim="\t",
  #                           callback = readr::DataFrameCallback$new(search.keyword),
  #                           chunk_size = 1000)
  #                             )
  #sc.intact = vroom(pipe("curl ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt | grep -w cerevisiae"),)
  if(from.MACOS){ intact.data=gsub("/media","/Volumes",intact.data) }

  IntAct = read.csv(intact.data,sep = "\t", quote = "", stringsAsFactors = F)
  colnames(IntAct) = c("protA","protB","altA","altB","aliasA","aliasB",
                       "method","pub.author1","pub.id",
                       "taxA","taxB","typeAB",
                       "src.db","int.id","conf",
                       "expansion.method","bio.roleA","bio.roleB",
                       "exp.roleA","exp.roleB",
                       "typeA","typeB",
                       "xrefA","xrefB","xrefAB",
                       "annotA","annotB","annotAB",
                       "host","parAB",
                       "created","updated",
                       "checksumA","checksumB","checksumAB",
                       "negative","featA","featB",
                       "stoichioA","stoichioB",
                       "id.met.protA","id.met.protB")
  keep.cols = colnames(IntAct)
  if(rm.useless){
    cat("-> rm columns that are useless...\n")
    cols.to.remove = c("created","updated","checksumA","checksumB","checksumAB",
                       "aliasA","aliasB","taxA","taxB","host",
                       "featA","featB","annotA","annotB",
                       "xrefA","xrefB","xrefAB","expansion.method")
    keep.cols = setdiff(colnames(IntAct),cols.to.remove)
  }
  INTACT.0 = IntAct[,keep.cols]
  cat(sprintf("(0) TOTAL INTERACTIONS : %s \n",nrow(INTACT.0)))
  # Highligh PPIs obtained with one of the methods given in input
  #INTACT$reliable = INTACT$method %in% used.methods


  # 1. remove small molecules and non-protein
  cat("-> rm EBI and small molecules \n")
  INTACT.1 =   INTACT.0[!(grepl("EBI",   INTACT.0$protA) | grepl("EBI",   INTACT.0$protB)),]
  cat(sprintf("(1) TOTAL INTERACTIONS : %s \n",nrow(INTACT.1)))

  sort_interactions = function(ppi,pair=c('protA','protB')){
    # Remove duplicated pairs and make undirected graph
    sorted.pair = t(apply(ppi[,pair], 1, sort)) # sort within pair
    ppi[,pair] = sorted.pair
    cat("-> sort pairs...\n")
    ppi$is.dup = duplicated(sorted.pair)
    cat("-> check duplicates...\n")
    ord.pair = order( paste0(sorted.pair[,1],"_",sorted.pair[,2]) )
    ppi = ppi[ord.pair,]
    cat("-> sort lexically...\n")
    # Make the is.dup the first column
    col_is.dup <- which(colnames(ppi)=="is.dup")
    ppi.sorted <- ppi[,c(col_is.dup,1:(ncol(ppi)-1))]
    cat("-> change columns order...\n")
    return(ppi.sorted)
  }

  # 2. sort the PPIs, so that B-A becomes A-B
  INTACT.2 = sort_interactions(ppi=INTACT.1)
  cat(sprintf("(2) TOTAL INTERACTIONS : %s \n",nrow(INTACT.2)))

  # 3. keep only uniprot/orf (redundant with step 1.)
  keep.uniprot.orf = function(ppi,pair=c('protA','protB')){
    uni_regex="([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
    orf_regex="([Y][A-P][LR][0-9]{3}[WC](?:-[A-Z])?)|(Q[0-9]{4})|(R[0-9]{4}[WC])"
    prot_regex = paste0(uni_regex,"|",orf_regex)

    p1 = grepl(prot_regex, x=ppi[,pair[1]] )
    p2 = grepl(prot_regex, x=ppi[,pair[2]] )

    ppi.prot = ppi[p1 & p2,]
    ppi.prot[,pair[1]] = gsub("^[a-z]+:","",ppi.prot[,pair[1]])
    ppi.prot[,pair[2]] = gsub("^[a-z]+:","",ppi.prot[,pair[2]])

    return(ppi.prot)
  }
  cat("-> Filter UniProt/ORF proteins...\n")
  INTACT.3 = keep.uniprot.orf(INTACT.2)
  cat(sprintf("(3) TOTAL INTERACTIONS : %s \n",nrow(INTACT.3)))

  # 4. remove duplicated interactions (with different Pubmed ID)
  INTACT.4 = INTACT.3[!duplicated(INTACT.3[,c("pub.id","protA","protB")]),]
  cat("-> Filter duplicated interactions within PMID...\n")
  cat(sprintf("(4) TOTAL INTERACTIONS : %s \n",nrow(INTACT.4)))

  if(only.direct){
    get.biophysical.methods = function(){
      if(!require(rols)){ BiocManager::install("rols") }
      library(rols)
      MI <- Ontology("MI")
      # MI:0013 Biophysical methods
      biophysical.methods = descendants(term(MI,"MI:0013"))
      return( as(biophysical.methods, "data.frame") )
    }


    # 5. filter for direct method
    biophy_meth= as.character(get.biophysical.methods()['id'])
    biophy_meth_regex= paste0(biophy_meth,'"',collapse="|")
    biophy_meth =

    #c("MI:0114","MI:0276","MI:0071","MI:0028","MI:0808","MI:0020","MI:0826","MI:0016","MI:0038","MI:0397")
    cat("-> Filter direct physical method...\n")
    INTACT.5 = INTACT.4[grep(biophy_meth_regex,x=INTACT.4$method), ]
    cat(sprintf("(5) TOTAL INTERACTIONS : %s \n",nrow(INTACT.5)))
  }else{
    INTACT.5 = INTACT.4
  }
  return(INTACT.5)
}

# Quaternary structures (3d-complex) -------------------------------------------
load.3dcomplex.yeast = function(limit = F, n = 1000) {
  library(tictoc)
  library(RMySQL)
  doing="Get 3d complex yeast protein structures by residues..."
  message(doing)
  ## Establish a connection with mySQL database
  V6 = dbConnect( MySQL(), user = "elevy", password = "Mysql1!", dbname = "3dcomplexV6", host = "mata")
  ## Prepare mySQL query
  # A) SELECT = combine vector of fields for each table
  chain = sprintf( "CH.%s",
                   c("org_ref", "uniprotAcc", "gene_name", "seqid", "overlap", "ident", "length_atom", "length_full")
  )
  residue_new = sprintf( "R.%s",
                         c("code", "chain", "resnum", "letter", "resname", "resnum_fasta", "asa_rel_alone", "asa_rel_in_BU", "`r4s.score`",
                           "`ali.nb.r4s`", "`ali.nb.seq`"
                         )
  )
  propensities = sprintf("P.%s", c("`hydro_kyte`", "`stickiness_ec`"))
  complex = sprintf("C.%s",
                    c("resol", "accession", "nsub2", "sym", "best_BU", "QSBIO_err_prob", "homo")
  )
  select.fields = toString(c(residue_new, propensities, chain, complex))
  # B) FROM = string with table names (including join statement)
  from.tables = " residue_new AS R
  JOIN chain as CH ON R.code = CH.code AND R.chain = CH.chain_name
  JOIN complex as C ON R.code = C.code
  JOIN propensities as P ON R.letter = P.letter"
  # C) WHERE = string with filtering criteria on fields
  cond = "CH.org_ref='sc'"
  # D) Build final query as a string
  query = sprintf("SELECT\n%s\nFROM\n%s\nWHERE\n%s\n", select.fields, from.tables, cond)

  ## Run mySQL query if test on 1000 rows is successfully passed
  Qtest = dbGetQuery(V6, paste0(query, " LIMIT ", n))
  if (limit) { treturn(Qtest) }

  if (exists("Qtest") && nrow(Qtest) > 0) {
    tic(doing)
    Q = dbGetQuery(V6, query)
    toc()
  } else{
    stop(sprintf("Something went wrong with your query!\n\n%s\n", query))
  }
  return(Q)
  # MISSING RESIDUES
  # V5 = dbConnect(MySQL(), user="elevy", password="Mysql1!", dbname="3dcomplexV5", host="mata")
  # chain = sprintf("CH.%s",c("chain_name", "seqid", "overlap", "ident","length_atom", "length_full", "gene_name", "org_ref", "abundance","abundanceMax", "abundanceID"))
  # missing = sprintf("M.%s",c("code", "resnum", "resname","fullAA", "atomAA", "type"))
  # residue_new = sprintf("R.%s",c("resnum", "resnum_fasta", "letter", "asa_rel_in_BU", "`r4s.score`" ))
  # select.fields = toString(c(chain,missing,residue_new))
  # cond = "M.code = CH.code AND M.chain = CH.chain_name and CH.org_ref='sc'"
  # Qm = dbGetQuery(V5,query)
}

get.mapping.3dcomplex.yeast = function(limit = F, n = 1000) {
  library(tictoc)
  library(RMySQL)
  doing="Get mapping between uniprot and pdb id for yeast protein complexes..."
  message(doing)
  #'select C.code, CH.chain_name, resol, seqid, ident, overlap from complex C, chain CH where C.code = CH.code and org_ref = "sc" and ident > 90 and CH.overlap > 20 and resol < 2.5 order by CH.overlap DESC'

  ## Establish a connection with mySQL database
  V6 = dbConnect( MySQL(), user = "elevy", password = "Mysql1!", dbname = "3dcomplexV6", host = "mata")
  ## Prepare mySQL query
  # A) SELECT = combine vector of fields for each table

  chain = sprintf( "CH.%s",
                   c("org_ref", "code", "chain_name", "seqid", "ident", "overlap", "length_atom", "length_full")
  )
  complex = sprintf("C.%s",
                    c("resol", "accession")#"nsub2", "sym", "best_BU", "QSBIO_err_prob", "homo")
  )
  select.fields = toString(c(chain, complex))
  # B) FROM = string with table names (including join statement)
  from.tables = "complex as C, chain  as CH"
  #JOIN chain as CH ON CH.code = C.code
  #"
  # C) WHERE = string with filtering criteria on fields
  cond = "CH.code = C.code AND CH.org_ref='sc' AND ident>90 AND CH.overlap > 20 AND resol < 3
  ORDER BY CH.overlap DESC
  "
  # D) Build final query as a string
  query = sprintf("SELECT\n%s\nFROM\n%s\nWHERE\n%s\n", select.fields, from.tables, cond)
  ## Run mySQL query if test on 1000 rows is successfully passed
  Qtest = dbGetQuery(V6, paste0(query, " LIMIT ", n))
  if (limit) { return(Qtest) }

  if (exists("Qtest") && nrow(Qtest) > 0) {
    tic(doing)
    Q = dbGetQuery(V6, query)
    toc()
  } else{
    stop(sprintf("Something went wrong with your query!\n\n%s\n", query))
  }
  return(Q)
}






