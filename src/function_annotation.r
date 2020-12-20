source('utils.r')
# Biological annotations (mapped to Uniprot) -----------------------------------
get.uniprot.pmid = function(uniprot) {
  require(AnnotationDbi)
  require(org.Sc.sgd.db)
  require(dplyr)
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
  require(AnnotationDbi)
  require(org.Sc.sgd.db)
  require(dplyr)
  require(tidyr)
  require(GO.db)
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
  require(AnnotationDbi)
  require(org.Sc.sgd.db)
  require(dplyr)
  require(GO.db)

  get.sgd_annot = function(desc, annotation = "", sep = ";") {
    require(stringr)
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
    mutate(ROLE = get.sgd_annot(DESCRIPTION, annot = "(role| act|require)", sep = ";,")) %>%
    mutate(LOC = get.sgd_annot(DESCRIPTION, annot = "(localize[sd]|localization)"))  %>%
    mutate(ORTHO = get.sgd_annot(DESCRIPTION, annot = "(ortholog|paralog|homolog)")) %>%
    mutate(COMPLEX = get.sgd_annot(DESCRIPTION, annot = "(complex|subunit|oligomer|mer)"))

  sgdinfo$OTHER = sapply(1:nrow(sgdinfo), function(x) {
    parse.sgd_annot(sgdinfo$DESCRIPTION[x],
                    exclude = sgdinfo[x, c("FUNCTION", "ROLE", "LOC", "ORTHO", "COMPLEX")])
  })

  return(sgdinfo %>% dplyr::select(-DESCRIPTION))
}

get.mapping.sgd.to.uniprot = function(tax=559292,input=id_sgd){
  # Get mapping for uniprot identifiers not included in UniRef proteome
  # Check that duplicated sequences are not associated to different uniprot identifiers
  require(tictoc)
  require(UniProt.ws)
  uniprot  = UniProt.ws::UniProt.ws(taxId=tax) # 559292 S. cerevisiae 288C (maintained by SGD)
  tic('Retrieving mapping between uniprot accession and SGD ID...')
  # TAKES ~50sec for about 7000 ids
  mapped = UniProt.ws::select(x=uniprot,
                              keys=unique(input),  keytype = "SGD",
                              columns = c("SGD","UNIPROTKB")) %>%
    filter(!is.na(UNIPROTKB)) %>% distinct()
  toc(log=T)
  return(mapped)
}

load.sgd.to.uniprot = function(tax=559292,input_id=id_sgd){
  if(missing(input_id)){ id_sgd = unique(load.sgd.features()[['sgdid']]) }
  # Download uniprot information (sequence,annotations...) for mapped SGD identifiers
  sgd_mapped = get.mapping.sgd.to.uniprot(input=input_id)
  require(UniProt.ws)
  uniprot  = UniProt.ws::UniProt.ws(taxId=tax) # 559292 = S. cerevisiae 288C (maintained by SGD)
  tic('Retrieving infos (sequence,annotations...) for SGD mapped uniprot entry...')
  # TAKES ~220sec for about 7000 ids
  sgd2uni_info = UniProt.ws::select(
    x=uniprot, keys=sgd_mapped$UNIPROTKB, keytype = "UNIPROTKB",
    columns = c("SGD","UNIPROTKB","REVIEWED","EXISTENCE","SCORE","LENGTH","FAMILIES",
                "PROTEIN-NAMES","SUBCELLULAR-LOCATIONS","COMMENTS","KEYWORDS",
                "SEQUENCE")
  )
  sgd2uni_info = sgd2uni_info %>%
    dplyr::rename(SUBLOC =`SUBCELLULAR-LOCATIONS`, PNAMES = `PROTEIN-NAMES`) %>%
    transmute(L = as.numeric(LENGTH))
  toc(log=T)
  SGD2UNI = inner_join(sgd_mapped,sgd2uni_info) %>% distinct()
  return(SGD2UNI)
}

parse.uniprot.subcellular_locations = function(subloc){
  require(stringr)
  require(tidyverse)
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
    filter(!is.na(loc)) %>%
    distinct() %>%
    mutate( HAS_FOCI = str_detect(string = SUBLOC, pattern = punc) ) %>%
    mutate( HAS_ISOFORM = str_detect(string = SUBLOC, pattern = isoform) ) %>%
    dplyr::select(-SUBLOC)

  return(uni.loc)
}

get.uniprot.localization = function(annot,loc_to_columns=T){
  if(missing(annot)){
    message(sprintf("Requires annotation table with 2 columns:\n(1) UNIPROTKB/SGD identifier\n(2) Uniprot subcellular location annotations"))
    stop("Use 'load.sgd.to.uniprot()' function to load annotations from uniprot for yeast.")
  }
  require(tidyverse)
  require(stringr)
  require(hablar)
  uni.annot = annot %>%
    # dplyr::rename_if( SUBLOC = `SUBCELLULAR-LOCATIONS` ) %>%
    filter(one2one) %>%
    dplyr::select(SGD,UNIPROTKB,EXISTENCE,SCORE,FAMILIES,SUBLOC) %>%
    hablar::convert( hablar::chr(SCORE) ) %>%
    mutate( SCORE = fct_relabel(SCORE, function(x){ gsub(" out of ","/",x) }),
            has_LOC = !is.na(SUBLOC))

  # GET UNIPROT ANNOTATION OF SUBCELLULAR LOCATIONS AS NAMED VECTOR
  subloc = uni.annot  %>% filter(has_LOC) %>% dplyr::select(UNIPROTKB,SUBLOC) %>% deframe()
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


