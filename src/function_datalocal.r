source("src/utils.r",local = T)
source("src/function_sequence.r",local = T)

# Local proteome data ----------------------------------------------------------
load.emmanuel.data = function(toolbox="/data/elevy/70_R_Data/bin/RToolBox_yeast_general.R"){
  source(toolbox,local = T)
  require(AnnotationDbi)
  require(org.Sc.sgd.db)
  require(GO.db)
  require(tidyverse)

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

load.1011.strains = function(seqdir){
  if( !dir.exists(seqdir) ) stop("Directory of proteome sequences not found!")
  fastas =list.files(path=seqdir, pattern = "fasta", full.names = T, ignore.case = T, include.dirs = F)
  return(read.proteomes(fastas))
}

# Evolution Sequence/Structure -------------------------------------------------
# Precomputed data #
#==================#
# PROTEOME=readRDS('./data/PROTEIN-EVORATE.rds')
# ALIGNED.DATA = readRDS('./data/RESIDUE-EVORATE.rds')
#==================#
load.wapinsky2007.data = function(path.data="./data/"){
  require(tictoc)
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
  message(sprintf('Surface is ASA above %s%%.',asa))
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
  require(rjson)
  # id='YDR134C'
  d2p2.url = sprintf('http://d2p2.pro/api/seqid/["%s"]',id)
  res = rjson::fromJSON(readLines(d2p2.url, warn=FALSE))
  nodata = length(res[[id]]) == 0
  if( !nodata ){
    # response$YCL051W[[1]][[3]]$disorder$consensus
    d = res[[id]]
    sp = d[[1]][[2]]
    if(!quiet){ message(sprintf("= %s from %s =",id,sp)) }
    consensus = d[[1]][[3]]$disorder$consensus
    record = list(id=id,from=sp,diso=d[[1]][[3]]$disorder, structure=d[[1]][[3]]$structure)
    return(record)
  }else{
    message(sprintf("id='%s' not recognized!",id))
    return(list(id=id,from=NULL,diso=NULL,structure=NULL))
  }
}

load.d2p2 = function(ids){ # Get the d2p2 predictions for multiple ids
  if(!is.vector(ids)){
    stop('Input must be a vector of identifiers as character...')
  }
  N=length(ids)
  d2p2 = list()
  for( i in 1:N ){
    ID = ids[i]
    prg = sprintf("%.1f %% [%5s/%5s]        \r", 100*(i/N), i,N)
    d2p2[[ID]] = fetch.d2p2(ID,quiet=T)
    cat(prg)
  }
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

load.codon.usage = function(inputseq,
                            codonR="/data/benjamin/NonSpecific_Interaction/Data/Evolution/eggNOG/codonR/",
                            exclude.counts=T
){
  message("REF: M. Dos Reis et al., 2004, Nucleic Acids Research")
  message("Solving the riddle of codon usage preferences: a test for translational selection")
  # https://doi.org/10.1093/nar/gkh834
  url.codonR="https://raw.githubusercontent.com/mariodosreis/tai/master/R/tAI.R"
  ord.codonR = get.codons4tai()

  if( !require(tAI) ){ devtools::install_github("mariodosreis/tai") } # Install this package first
  require(tAI)
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
  trna.genes=sprintf("%s/tRNA/4932_sc.trna",codonR)
  trna.raw = read.delim(trna.genes, col.names=c('AA3','tot','cod','count','acod'),
                        header=F,comment.char = "#", sep = '\t', stringsAsFactors = F)
  trna.ord = match( ord.codonR, trna.raw$acod )
  trna = trna.raw[trna.ord,]

  ws =  tAI::get.ws(tRNA = trna$count, sking = 0)
  tai = tAI::get.tai(codons.60, ws)
  res = data.frame(prot=names(dna),tAI=tai)
  if(!exclude.counts){ res = data.frame(res,codons) }
  return(res)
}

# Amino Acid features ----------------------------------------------------------
# path$scales = "/media/elusers/users/emmanuel/PAPERS/013_disorder_stickiness/data/Scales/"
get.scales = function(SCALE_PATH = path$scales, byAA = T) {
  scale.files = list.files(SCALE_PATH, pattern = 'scale.csv')
  names(scale.files) = subname(scale.files, sep = "\\.", lc = T)
  L.scales = lapply(scale.files, function(x) { read.csv(paste0(SCALE_PATH, x), stringsAsFactors = F) })
  df.scales = Reduce(f = merge, x = L.scales)
  if (!byAA) {
    numcol = sapply(df.scales, is.numeric)
    aa.val = t(df.scales[, numcol])
    colnames(aa.val) = df.scales$AA
    df.scales = data.frame(aa.val)
  }
  return(df.scales)
}

# Quaternary structures (3d-complex) -------------------------------------------
load.3dcomplex.yeast = function(limit = F, n = 1000) {
  require(tictoc)
  require(RMySQL)
  message("Get 3d complex yeast protein structures by residues...")
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
  if (limit) { return(Qtest) }

  if (exists("Qtest") && nrow(Qtest) > 0) {
    Q = dbGetQuery(V6, query)
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
  require(tictoc)
  require(RMySQL)
  message("Get mapping between uniprot and pdb id for yeast protein complexes...")
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
    Q = dbGetQuery(V6, query)
  } else{
    stop(sprintf("Something went wrong with your query!\n\n%s\n", query))
  }
  return(Q)
}

