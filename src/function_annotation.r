#source('src/utils.r',local = T)

load.sgd.features = function(by.chr=T){ # Gene/Protein features from SGD
  library(stringr)
  sgd_feat.url = "http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
  sgd.feat = read.delim2(sgd_feat.url, sep='\t', quote = "",
                         header=F, fill=T, strip.white=T,stringsAsFactors = F,
                         col.names = c('sgdid','type','qual','name',
                                       'gname','alias','parent','sgdid2',
                                       'chr','start','end','strand','gpos',
                                       'coordv','seqv','desc')  )
  ord = 1:nrow(sgd.feat)
  if(by.chr){ ord = gtools::mixedorder(sgd.feat$chr) }
  return(sgd.feat[ord,])
}

load.sgd.orf = function(sgd.feat){
  if(missing(sgd.feat)){ sgd.feat = load.sgd.features() }
  ORF             = sgd.feat$type == 'ORF'
  #Verified        = sgd.feat$qual == 'Verified'
  #Uncharacterized = sgd.feat$qual == 'Uncharacterized'
  #notDubious      = c(Verified,Uncharacterized)
  return( sgd.feat$name[ORF] )
}

# Molecular Interaction controlled vocabulary
get.MI.annotation= function(id="MI:0013",
                            relation=c('descendants','children','ancestors','parents','siblings'),
                            include=T){
  library(rols)
  ol = rols::Ontologies()
  MI <- ol[['mi']]
  ID = term(MI,id)
  related = match.arg(relation, choices =c('descendants','children','ancestors','parents'), several.ok = F )
  message("looking for ",relation," of ",id," [",termLabel(ID),"] ...")
  if(related == 'children'   ){ MI.annot = children(ID)    }
  if(related == 'descendants'){ MI.annot = descendants(ID) }
  if(related == 'ancestors'  ){ MI.annot = ancestors(ID)   }
  if(related == 'parents'    ){ MI.annot = parents(ID)     }
  if(related == 'siblings'    ){ MI.annot = children(parents(ID)[[1]])     }
  if(include){ return( rbind( as(ID,"data.frame"), as(MI.annot, "data.frame") ) ) }
  return( as(MI.annot, "data.frame") )
}

# Biological annotations (mapped to Uniprot) -----------------------------------

load.uniprot.features = function(tax=559292,refdb='UNIPROTKB'){ # Gene/Protein features from Uniprot
  library(UniProt.ws)
  # 559292 S. cerevisiae 288C (maintained by SGD)

  uniprot  = UniProt.ws::UniProt.ws(taxId=tax)
  refkey=match.arg(arg=refdb,choices=UniProt.ws::keytypes(uniprot),several.ok = F)

  doing=sprintf('Retrieving uniprot features for %s ID...',refdb)
  message(doing)
  # TAKES ~50sec for about 7000 ids
  tic(doing)

  ids = keys(uniprot, refdb)
  FEATURES = c("SGD","UNIPROTKB","REVIEWED","EXISTENCE","SCORE","LENGTH",
               "FAMILIES","FEATURES","PATHWAY","DOMAIN","DOMAINS","INTERACTOR","GO-ID","GO","GENES",
               "PROTEIN-NAMES","SUBCELLULAR-LOCATIONS","COMMENTS","KEYWORDS")

  up.feat = UniProt.ws::select(x=uniprot, keys=ids,  keytype = refkey,
                               columns = FEATURES) %>%
    dplyr::filter(!is.na(UNIPROTKB)) %>%
    dplyr::rename(SUBLOC =`SUBCELLULAR-LOCATIONS`, PNAME = `PROTEIN-NAMES`) %>%
    mutate(L = as.numeric(LENGTH)) %>% dplyr::select(-LENGTH) %>%
    distinct()
  toc(log=T)
  return(up.feat)
}

get.uniprot.pmid = function(uniprot) {
  library(AnnotationDbi)
  library(org.Sc.sgd.db)
  library(dplyr)
  pmid = AnnotationDbi::select(x = org.Sc.sgd.db,
                               columns = c('PMID'),
                               keys = uniprot, keytype = 'UNIPROT') %>%
    group_by(UNIPROT) %>%
    mutate(n_pub = n_distinct(PMID), PMIDS = toString(PMID)) %>%
    dplyr::select(-PMID) %>%
    distinct()
  return(pmid)
}

get.uniprot.go = function(uniprot) {
  library(AnnotationDbi)
  library(org.Sc.sgd.db)
  library(dplyr)
  library(tidyr)
  library(GO.db)
  go = AnnotationDbi::select(x = org.Sc.sgd.db,
                             columns = c('GO', 'ONTOLOGY', 'EVIDENCE'),
                             keys = uniprot, keytype = 'UNIPROT'
  ) %>%
    mutate(obsolete = GO %in% keys(GOOBSOLETE)) %>%
    filter(!obsolete & !is.na(GO)) %>%
    mutate(goterm = Term(GO), onto=ONTOLOGY) %>% arrange(UNIPROT,GO) %>%
    group_by(UNIPROT) %>% mutate( ALL= n_distinct(GO) ) %>%
    group_by(UNIPROT,ONTOLOGY) %>% add_count(name='ngo') %>%
    tidyr::pivot_wider(names_from=onto, values_from=ngo, values_fn=list(ngo=sum), values_fill=list(ALL=0)) %>%
    group_by(UNIPROT) %>% fill(MF,CC,BP,.direction = 'downup') %>%
    group_by(GO) %>% mutate(shared=n_distinct(UNIPROT)) %>%
    group_by(UNIPROT,GO) %>% mutate( EVI_ALL = str_c(EVIDENCE,collapse = '|'))

  goevidence = go %>% mutate(checked = TRUE) %>% ungroup() %>%
    pivot_wider( id_cols = c(UNIPROT, GO), names_from = EVIDENCE,
                 values_from = checked, values_fn = list(checked = sum)
    ) %>% mutate( n_evi = rowSums(dplyr::select(.,IDA:RCA),na.rm=T) )

  goinfo = dplyr::select(go,-EVIDENCE) %>%
    left_join(goevidence %>% dplyr::select(-c(IDA:RCA))) %>%
    mutate(has_biodata = replace_na(ND<1,replace = TRUE)) %>%
    dplyr::select(-ND)%>%
    distinct() %>% ungroup()
  return(goinfo)
}

get.uniprot.sgd = function(uniprot) {
  library(AnnotationDbi)
  library(org.Sc.sgd.db)
  library(dplyr)
  library(GO.db)

  get.sgd_annot = function(desc, annotation = "", sep = ";") {
    library(stringr)
    found = str_extract_all(tolower(desc), pattern = sprintf("[^%s]*%s[^%s]+", sep, annotation, sep))
    sapply(found, paste, collapse = ";")
  }

  parse.sgd_annot = function(desc, exclude = c()) {
    annotations = unlist(strsplit(x = tolower(desc), split = ";"))
    excluded = unlist(strsplit(x = paste(exclude), split = ";"))
    paste(annotations[!annotations %in% excluded], collapse = ";")
  }

  sgdinfo = AnnotationDbi::select( x=org.Sc.sgd.db,
                                   columns = c('ORF', 'GENENAME', 'DESCRIPTION'),
                                   keys = uniprot, keytype = 'UNIPROT') %>%
    distinct() %>%
    mutate(FUNCTION = tolower(str_extract(DESCRIPTION, "[^;]+"))) %>%
    mutate(ROLE = get.sgd_annot(DESCRIPTION, annot = "(role| act|library)", sep = ";,")) %>%
    mutate(LOC = get.sgd_annot(DESCRIPTION, annot = "(localize[sd]|localization)"))  %>%
    mutate(ORTHO = get.sgd_annot(DESCRIPTION, annot = "(ortholog|paralog|homolog)")) %>%
    mutate(COMPLEX = get.sgd_annot(DESCRIPTION, annot = "(complex|subunit|oligomer|mer)"))

  sgdinfo$OTHER = sapply(1:nrow(sgdinfo), function(x) {
    parse.sgd_annot(sgdinfo$DESCRIPTION[x],
                    exclude = sgdinfo[x, c("FUNCTION", "ROLE", "LOC", "ORTHO", "COMPLEX")])
  })

  return(sgdinfo %>% dplyr::select(-DESCRIPTION))
}


parse.uniprot.subcellular_locations = function(subloc){
  library(stringr)
  library(tidyverse)
  # subloc should be a named vector
  # Names = identifiers
  # Values = subcellular locations uniprot annotation

  # UNIPROT SUBCELLULAR LOCATIONS contains extra information such as:
  # Annotation always starts with "SUBCELLULAR LOCATION:"
  # Annotation of isoforms/subunits are preceded by square brackets "[]"
  # For isoform "SUBCELLULAR LOCATION:" is repeated within the annotation
  # Evidence Code and PMID wrapped around curly brackets "{ECO...|PMID...}"
  # The sentence followed by "Note=" may contain additional information about localization
  loctag  =  "(SUBCELLULAR LOCATION: )"
  evidence = "(\\{[^\\{\\}]+\\})"
  isoform = "(\\[[^\\[\\]]+\\]:)"
  note = "(Note\\=.+\\.)"
  end_note = sprintf("(%s$)",note)
  punc='(insoluble|aggregate|foci|punctate|granules)'

  LOC = subloc  %>%
    # ERASE NON-ESSENTIAL PARTS OF ANNOTATION (keeping terms corresponding to subcellular compartments)
    str_replace_all(pattern=loctag,replacement="") %>%
    str_replace_all(pattern=evidence,replacement="") %>%
    str_replace_all(pattern=end_note,replacement="") %>%
    str_replace_all(pattern=punc,replacement="") %>%
    str_replace_all(pattern=isoform,replacement="") %>%
    str_replace_all(pattern=paste0(note,"(?=[;,])"),replacement="") %>%
    str_to_lower() %>%
    # REMOVE SPACE IN BETWEEN SEPARATORS OF ANNOTATIONS
    str_replace_all(pattern="(?<=[;.,])( +)|( +)(?=[;.,])",replacement="") %>%
    # SPLIT BY TERMS (hopefully subcellular compartments only)
    stringi::stri_split_regex(pattern="[;.,]",omit_empty = T) %>%
    # GIVE BACK THE PROTEINID TO EACH SET OF TERMS
    set_names(nm=names(subloc))

  uni.loc = enframe(LOC,name='id',value='loc') %>% bind_cols(SUBLOC = subloc) %>%
    unnest(cols=c(loc)) %>%
    mutate( loc = str_replace_all( str_trim(loc,side="both"), c("-" = ".", "[ ]+" = "_")) ) %>%
    dplyr::filter(!is.na(loc)) %>%
    distinct() %>%
    mutate( HAS_FOCI = str_detect(string = SUBLOC, pattern = punc) ) %>%
    mutate( HAS_ISOFORM = str_detect(string = SUBLOC, pattern = isoform) ) %>%
    dplyr::select(-SUBLOC)

  return(uni.loc)
}

parse.uniprot.pathways = function(pathways){
  library(stringr)
  library(tidyverse)
  # pathways should be a named vector
  # Names = identifiers
  # Values = pathways from uniprot annotation

  # separate_rows(UNI.ALL[,c('UNIPROTKB','PATHWAY')], PATHWAY,sep="PATHWAY:  ")

  # UNIPROT PATHWAYS contains extra information such as:
  # Annotation always starts with "PATHWAY:"
  # Evidence Code and PMID wrapped around curly brackets "{ECO...|PMID...}"
  pathtag  =  "(PATHWAY: )"
  evidence = "(\\{[^\\{\\}]+\\})"
  #pathways = split(UNI.ALL$PATHWAY[!duplicated(UNI.ALL$UNIPROTKB)], UNI.ALL$UNIPROTKB[!duplicated(UNI.ALL$UNIPROTKB)])

  PATH = pathways  %>%
    # ERASE NON-ESSENTIAL PARTS OF ANNOTATION (keeping terms corresponding to subcellular compartments)
    str_replace_all(pattern=pathtag,replacement="") %>%
    str_replace_all(pattern=evidence,replacement="") %>%
    str_to_lower() %>%
    # REMOVE SPACE IN BETWEEN SEPARATORS OF ANNOTATIONS
    #str_replace_all(pattern="(?<=[;.,])( +)|( +)(?=[;.,])",replacement="") %>%
    # SPLIT BY TERMS (hopefully subcellular compartments only)
    #stringi::stri_split_regex(pattern="[;.,]",omit_empty = T) %>%
    # GIVE BACK THE PROTEINID TO EACH SET OF TERMS
    set_names(nm=names(pathways)) %>%
    purrr::discard(is.na)

  uni.path = enframe(PATH,name='id',value='path') %>% bind_cols(PATH = path) %>%
    unnest(cols=c(path)) %>%
    mutate( path = str_replace_all( str_trim(path,side="both"), c("-" = ".", "[ ]+" = "_")) ) %>%
    dplyr::filter(!is.na(path)) %>%
    distinct() %>%
    dplyr::select(-SUBLOC)

  return(uni.loc)
}

get.uniprot.localization = function(annot,loc_to_columns=T){
  if(missing(annot)){
    stop("Use 'load.uniprot.features()' function to load annotations from uniprot for yeast.")
  }
  library(tidyverse)
  library(stringr)
  library(hablar)
  uni.annot = annot %>%
    # dplyr::rename_if( SUBLOC = `SUBCELLULAR-LOCATIONS` ) %>%
    filter(one2one) %>%
    dplyr::select(SGD,UNIPROTKB,EXISTENCE,SCORE,FAMILIES,SUBLOC) %>%
    hablar::convert( hablar::chr(SCORE) ) %>%
    mutate( SCORE = fct_relabel(SCORE, function(x){ gsub(" out of ","/",x) }),
            has_LOC = !is.na(SUBLOC))

  # GET UNIPROT ANNOTATION OF SUBCELLULAR LOCATIONS AS NAMED VECTOR
  subloc = uni.annot  %>% dplyr::filter(has_LOC) %>% dplyr::select(UNIPROTKB,SUBLOC) %>% deframe()
  uni.loc = parse.uniprot.subcellular_locations(subloc)

  if( loc_to_columns ){
    uni.isloc  = uni.loc %>%  mutate(seen=1) %>%
      pivot_wider( id_cols = c('id','HAS_FOCI','HAS_ISOFORM'),
                   names_from = 'loc',
                   values_from = 'seen', values_fn=list(seen = sum),values_fill = list(seen=0)) %>%
      group_by(id,HAS_FOCI,HAS_ISOFORM)
    return(uni.isloc)
  }
  return(uni.loc)
}


## _1_ protein families -------------------------------------------------------
url.similar="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/similar.txt"

## _2_ subcellular locations --------------------------------------------------
url.subcell="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/subcell.txt"
#@uni2loc = get.uniprot.localization(sgd2uni,loc_to_columns = F)
#uni_isloc = get.uniprot.localization(sgd2uni,loc_to_columns = T)
#LOC = unique(uni2loc$loc)
#uni_isloc %>% ungroup() %>% summarise_at(LOC,sum) %>% pivot_longer(everything()) %>% arrange(value) %>% filter(value > 30)

## _3_ pathways ---------------------------------------------------------------
get.uniprot.pathways = function(ORG="YEAST"){
  library(tictoc)
  library(tidyverse)
  tic("Reading uniprot pathways...")
  url.pathway="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pathway.txt"
  pathways.txt=readLines(url.pathway,)
  pathways.division=grep("\\.$",pathways.txt,v=T)
  pathways.members=grep("\\.$",pathways.txt,v=T)
  toc()

  if(!missing(ORG) & ORG != "" ){ message(sprintf("FILTER FOR: %s\n",ORG))}

  pw.members=c()
  pw.lvl=NULL
  begin=F
  tic("Parse Uniprot pathways...")
  PATHWAY = tribble(~name,~ids)
  for( pw in pathways.txt){
    is.division = grepl("\\.$",pw)
    regex.entry = "([A-Z0-9]{1,5}_[A-Z0-9]{1,5})"
    regex.AC = UNIPROT.nomenclature()
    is.break = grepl("^$",pw)
    is.members = grepl(paste0(regex.entry,"\\s+\\(",regex.AC,"\\)"),pw)#,"\\s+",regex.AC,",?\\s+"),pw)
    if( is.division ){
      begin=T
      pw.lvl = pw
      #cat(sprintf("START OF [%s]\n",paste0(pw.lvl,collapse=" ->")))
    }else if(is.members & begin){
      AC = str_extract_all(pw,regex.AC)
      ENTRY = str_extract_all(pw,regex.entry)
      is.org = grep(ORG,unlist(ENTRY))
      pw.members = c(pw.members,unlist(AC)[is.org])
    }else if(is.break & begin){
      #cat(sprintf("(%s)\n",paste0(pw.members,collapse=" ")))
      #cat(sprintf("END OF [%s]\n",paste0(pw.lvl,collapse=" ->")))
      PATHWAY=add_row(PATHWAY, name=pw.lvl, ids = ifelse(length(pw.members)>0,paste0(pw.members,collapse=","),NA))
      pw.lvl=NULL
      pw.members=c()
      begin=F
    }else{
      #cat("\n")
      #cat( sprintf("unrecognized format (no data?): %s \n",pw) )
    }
  }
  toc()
  return(PATHWAY)
}

get.KEGG = function(sp='sce',type=c('pathway','module'),as.df=F){

  library(KEGGREST)
  library(tidyverse)

  genes.grp = keggLink(type,sp)
  df1 =enframe(genes.grp, name="sp.id", value="grp") %>%
        mutate( id = sub(paste0(sp,":"),"",sp.id) )
  if( type == 'module'){
    df1$grp = sub(paste0(sp,"_"),"",df1$grp)
    sp = ""
  }

  grp.desc = keggList(type,sp)
  df2 =enframe(grp.desc, name = 'grp', value='desc')
  if(type == 'pathway'){ df2$desc  = sub(" -[^-]+$","",df2$desc) }

  ## returns gene-pathway
  df = left_join(df2,df1,by='grp') %>%
       filter(!is.na(id))
  if(as.df){ return(df) }
  return( split(df$id,df$desc,drop = T) )
}

#PP = get.uniprot.pathways(ORG='YEAST')
#test = separate(PP,name,sep = ";", into = c('div1','div2','div3')) %>% filter(!is.na(ids)) %>% group_by(div1,div2) %>% summarise( allids= paste(ids,collapse=","), n = str_count(string = allids,",")+1)
#  pathways.list = map( strsplit(pathways.division,";"), function(x){ x %>% str_replace("\\.","") %>% str_trim() })

#first = unique(unlist(map(pathways.list,pluck,1)))
#second = unique(unlist(map(pathways.list,pluck,2)))
#names(pathways.list) = unlist(  )

## _4_ functional categories --------------------------------------------------
load.yeast.biofunctions = function(){
  # load gene classification of biological functions (costanzo 2010, Van Leeuwen 2016)
  A = load.costanzo2010.data()
  B = load.vanleeuwen2016.data()
}
