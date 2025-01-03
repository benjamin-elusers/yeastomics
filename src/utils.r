# author: Benjamin Dubreuil
# created: 20/12/2020
# comments: general purpose functions
# ```{r}Sys.time()```
#
# I/O --------------------------------------------------------------------------

url.exists <- function(url) {
  # Test if url exists

  x <- httr::status_code(httr::HEAD(url))
  if (identical(x, 200L)) return(TRUE)
  cnt <- 1
  while (identical(x, 502L)) {
    ## server not response
    ## try again
    x <- httr::status_code(httr::HEAD(url))
    if (identical(x, 200L)) return(TRUE)
    cnt <- cnt + 1
    if (cnt > 10) break
  }
  return(FALSE)
}

safe_download = function(url, path, ext='.txt', debug=F){

  file = path
  if( dir.exists(path) ){
    filename = basename(url)
    file = file.path(path,paste0(filename,ext))
  }

  if( file.exists(file) && debug ){
    message(sprintf("already downloaded (%s)",file))
  }

  if(url.exists(url) && !file_exists(file)){
    tryCatch(
      download.file(url , file, mode = "w", quiet=T),
      error = function(e){
        message(sprintf("URL does not exist (%s)",url))
      })
  }
}

check.url <- function(url) {
  # Test if url exists and accessible

  idx <- vapply(url, url.exists, logical(1))
  url[!idx] <- NA
  return(url)
}

# Retrieve text content from URL
open.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con,skipNul=T)
  #closeAllConnections()
  close.connection(con)
  return(textConnection(txt))
}

read.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con,skipNul=T)
  #closeAllConnections()
  close.connection(con)
  return(txt)
}

# Function to determine file type and load appropriately
load_datafile <- function(datafile) {
  ext <- tools::file_ext(tolower(datafile))
  if (ext == "rds") {
    return(readRDS(datafile))
  } else if (ext %in% c("rda", "rdata", "data")) {
    load(datafile, envir = .GlobalEnv) # Loads into the global environment
    return(invisible(NULL)) # Return NULL since the data is in the environment
  } else {
    stop("Unsupported file format: ", ext)
  }
}

# Execute loading call and save data unless saved data already exists
preload <- function(saved.file, loading.call, doing = "Creating data...") {
  library(tictoc)
  cat(doing, "\n")
  
  # If the file doesn't exist, create and save it
  if (!file.exists(saved.file)) {
    tic(doing)
    res <- eval(substitute(loading.call))
    
    # Ensure the directory exists
    save_dir <- dirname(saved.file)
    if (!dir.exists(save_dir)) {
      message("Directory does not exist. Creating path: ", save_dir)
      dir.create(save_dir, recursive = TRUE)
    }
    
    # Save the result in RDS format
    saveRDS(res, saved.file)
    toc()
    return(res)
  } else {
    # Load the existing file
    return(load_datafile(saved.file))
  }
}

get.last.file = function(path,pattern){
  library(tidyverse)
  list.files(path, pattern, full.names=T) |>
    file.info() |>
    dplyr::slice(which.max(mtime)) |>
    rownames()
}

fallback = function(url,archived){
  # if current download is dead/blocked use the latest archive
  if( httr::http_error(url) ){
    warning(sprintf("falling back to last release archived (%s)...",archived),immediate. = T)
    url = archived
  }else{
    message("URL...OK")
  }
  return(url)
}

# Testing/Subsetting -----------------------------------------------------------
is.whole  <- function(x){ all(floor(x) == x) }      # checks if a value has decimal part (not necessarily integer e.g. 1.0 is whole)
is.string <- function(x){  length(x) == 1 & is.character(x) } # checks if a value is a single string of letters
#is.binary <- function(x){ all( (1*x) %in% 0:1) }    # checks if a value is 0/1 (binary/logical)

is_num_bin <- function(x){ is.numeric(x) & length(unique(na.omit(x))) == 2 } # checks if a numeric vector has only 2 values
is_fac_bin <- function(x){ is.factor(x) & nlevels(na.omit(x)) == 2 }         # checks if a factor vector has only 2 levels
# checks if a value is binary (2 unique outcomes -- by default numeric values must be either 0 and 1)
is_binary  <- function(x,xvals=c(0,1)){ is_num_bin(x) & all(range(na.omit(x)) %in% xvals[1:2]) | is_fac_bin(x) | is.logical(x) }
is_number <- function(x){ grepl("^[0-9]+$",x) } # checks if it contains only number
is_frequency <- function(x){ all(is.numeric(x) & x>=0 & x<=100) } # checks if values are between 0 and 100

is.even   <- function(x){ as.integer(x) %% 2 == 0 } # checks if a value is even
is.odd    <- function(x){ as.integer(x) %% 2 != 0 } # checks if a value is odd
is.dup    <- function(x){ x %in% x[duplicated(x)] } # detects duplicates

is.in = function(el,set,case.sensitive=T,withNames=T){
  found = el %in% set
  if(!case.sensitive){ found = toupper(el) %in% toupper(set) }
  if(withNames) return(setNames(found,el))
  return(found)
}

is.outlier <- function(x,thr=0.1,coef=1.5) { # return outliers outliers from distribution
  whiskers = coef * IQR(x,na.rm=T)
  return(x < quantile(x, thr,na.rm=T) - whiskers | x > quantile(x, 1-thr,na.rm=T) + whiskers)
}

find_na_rows = function(df,as.indices=F){
  # return dataframe rows with NA values (or just indices of NA row)
  NA_codon = rowsNA(df) > 0
  if(as.indices){ return(which(NA_codon)) }
  return(df[which(NA_codon),,drop = F]) # drop=F to keep the same dimensions (no coercion)
}

find.consecutive = function(x,val,minilen=2){  # return stretch of identical consecutive values in vector
  y=rle(x)
  V=y$values
  L=y$lengths

  seg = V == val & L >= minilen
  y$values[seg] = seq_along(V[seg])
  y$values[!seg] = 0

  return(inverse.rle(y))
}

# find.consecutive = function(x){
#   return (cumsum(c(T, diff(x) == 0L) & !x))
# }

mid       <- function(x){ ifelse(is.even(length(x)), length(x)/2, (length(x)+1)/2) }

evens     <- function(x){ x[seq(1,length(x),by=2)] } # returns values at even indexes
odds      <- function(x){ x[seq(2,length(x),by=2)] }  # returns values at odd indexes

is.last   <- function(x){ which(x) == length(x) }
last      <- function(x){ tail(x, n = 1) }             # returns last element of a vector
b4.last   <- function(x){ head( tail(x,n=2) , n = 1) } # returns element before last
but.last  <- function(x){ head(x, n = -1) }            # returns vector w/o last element
from.last <- function(x,n=2){ head(x, n = -abs(n)) }   # returns vector w/o last element

unfactor  <- function(x){ as.numeric(as.character(x)) } # forces conversion of factors to numbers
compact   <- function(x){ Filter(Negate(is.null), x) }  # removes null element in list
invert    <- function(x){ setNames(names(x),make.unique(x)) }

# Counting ---------------------------------------------------------------------
ulen    <- function(x){ return( length(unique(x)) ) } # gets length of unique values
normin  <- function(x){ return( scale(x,center=min(x,na.rm=T),scale=1)[,1] ) } # adds minimum to values (keeps original range)
norm01  <- function(x){ return( scale(x,center=min(x,na.rm=T),scale=diff(range(x,na.rm=T)))[,1] ) } # scales values between 0 and 1
maxrep  <- function(x,n){ sort(table(x), decreasing = TRUE)[1:n] } # gets the Nth most repeated values
sum.na  <- function (x,notNA=FALSE){ sum(is.na(x) == !notNA) } # returns sum of NA or not-NA values
geomean <- function(x) {  exp(mean(log(x[x != 0 & !is.na(x)]))) } # returns geometrical mean
geosd   <- function(x) {  exp(sd(log(x[x != 0 & !is.na(x)]))) } # returns geometrical standard deviation
min_above <- function(x,above=0,...){ min(x[x>above],...) }
max_below <- function(x,below=0,...){ max(x[x<below],...) }

rowsNA = function(m){ # counts how many rows have NAs
  # input must be 2D array (matrix/dataframe)
  nas = apply(as.matrix(m),1,sum.na)
  if(sum(nas)>0){ message(sprintf('%s rows contain NA values.',sum(nas>0))) }
  return(nas)
}

sum_mean <- function (x,na.rm=T,as.list=T){ # returns sum and mean with/out NA value
  res = list( 'n'= length(x), 'sum'=sum(x,na.rm=na.rm), 'mean' = mean(x,na.rm=na.rm), 'na' = sum(is.na(x)))
  if( !as.list ){  return(unlist(res)) }
  return(res)
}

getmode <- function(v) { # returns the mode (most frequent value)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# Strings ----------------------------------------------------------------------
str2chr  <- function(x){ return( unlist(strsplit(x,split='')) ) } # split string by character
concat   <- function(...){ return( paste(...,collapse='') ) }     # concatenate characters to string

towords  <- function(x){ return( unlist(strsplit(x,split="\\s")) ) } # split string by any white space
xxS      <- function(x,sx,s='.'){ paste0(x,s,sx) } # Add suffix to a string
Pxx      <- function(x,px,s='.'){ paste0(px,s,x) } # Add prefix to a string
starting <- function(str,pre){ str[ startsWith(str,pre) ] }
ending   <- function(str,suf){ str[ endsWith(str,suf) ] }

trim.lead  <- function (x){ sub("^\\s+", "", x) } # returns string w/o leading whitespace
trim.trail <- function (x){ sub("\\s+$", "", x) }# returns string w/o trailing whitespace
trim       <- function (x){ gsub("^\\s+|\\s+$", "", x) } # returns string w/o leading or trailing whitespace
trim.NaN   <- function(d){ d[apply(d,2,function(x) any(is.nan(x))),] } # returns a vector without NaN

paste.even <- function(x,s='-'){ paste(evens(x),odds(x),sep=s) }
paste.odd  <- function(x,s='-'){ c(paste.even(but.last(x),s), last(x)) }

rm.dup.str = function(str,sep='.'){ # removes duplicated word in string
  words = strsplit(str,paste0("\\",sep))
  words.unique = lapply(words,unique)
  new.str = unlist(lapply(words.unique,paste,collapse=sep))
  return(new.str)
}

repchar <- function(char,times){ # replicates a character multiple times as a string
  return(paste( rep(char, times), collapse='' ))
}

subname=function(name,sep="\\.",lc=F){ # extracts substring until first separator
  b4sep = sprintf("^([^%s]+)%s",sep,sep)
  part1 = sub(b4sep, "\\1", x=name)
  if(lc){ tolower(part1) }
  return(part1)
}

get.longest = function(S, s='\\.'){
  # get the longest string in a list of splitted string
  library(stringr)
  L = str_split(string = S, pattern = s)
  long=sapply(L,function(x){ nc=nchar(x); which.max(nc)})
  sapply(1:length(L),function(i){ L[[i]][long[i]] })
}

strfind = function(strings, patterns, index=F){ # search multiple patterns in character vectors
  sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = !index) })
}

file_ext <- function (fn) {
  # remove a path
  splitted    <- strsplit(x=fn, split=.Platform$file.sep)[[1]]
  # or use .Platform$file.sep in stead of '/'
  fn          <- splitted [length(splitted)]
  ext         <- ''
  splitted    <- strsplit(x=fn, split='\\.')[[1]]
  l           <-length (splitted)
  if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l]
  # the extention must be the suffix of a non-empty name
  ext
}

catn = function(x, ...){ cat(x,"\n",...) }


match_strings = function(SP1, SP2, max_strings=20, use_soundex=T, manual = F, verbose=T){
  # Matching strings based on similarity
  library(tidystringdist)
  library(stringdist)
  library(sjmisc)
  sp1 =  tolower(SP1)
  sp2 =  tolower(SP2)
  paired_name = expand_grid(s1=sp1,s2=sp2) |>
                # Keep identically matched names
                rowwise() |> mutate(is_identical = identical(s1,s2)) |>
                group_by(s1) |> filter( if_any(is_identical) | sum(is_identical)==0 )

  matched_name = paired_name |>
                 # Keep matched names that are nested (i.e. one string is a substring of the other)
                 rowwise() |> mutate(is_substring = str_contains(s2,s1) - str_contains(s1,s2)) |>
                 group_by(s1) |> filter( if_any(is_substring, ~ . != 0) | all(is_substring==0)) |>
                 mutate( is_matched = (n() == 1))

  twins = matched_name |> filter(is_matched) |> mutate(verified = 'by_substring')
  unmatched  = matched_name |> filter(!is_matched)

  similarities = tidy_stringdist(df=unmatched,v1=s1,v2=s2) |>
                 group_by(s1) |>
                 # remove names that have been previously matched (identical or substring)
                 mutate( already_matched = s2 %in% twins$s2, verified='by_similarity') |>
                 filter( !already_matched ) |>
                 filter( ifelse(use_soundex, soundex == min(soundex), T) ) |>
                 add_count(name='n1') |> mutate( is_matched = (n1 == 1) )  |>
                 slice_min(order_by = lcs, n = max_strings)

  if(verbose){
    cat(sprintf('%s identical/substring matches...\n',n_distinct(twins$s1)))
    cat(sprintf('%s unmatched names...\n',n_distinct(unmatched$s1)))
  }

  ndup = sum(similarities$n1 != 1)
  results = bind_rows(twins,similarities)
  if(manual && ndup>1){
    cat('matching remaining names manually...\n\n')
    DUP = similarities |> filter(n1 != 1)
    dup_name = DUP$s1 |> unique()
    name_chosen=c()
    for(name in dup_name){
      name_options = DUP |> filter(s1 == name) |> pull(s2) |> setdiff(name_chosen)
      n_options = n_distinct(name_options)
      x=1
      if(n_options > 1){
        x = menu(name_options,graphics = F, title = paste0("'",name,"' corresponds to:"))
      }
      name_chosen = append(name_chosen, name_options[x])
    }
    MANUAL = tibble(s1=dup_name,s2=name_chosen) |>
                left_join(DUP, by=c('s1','s2')) |>
                mutate(verified = 'manually/taxid/synonym')
    results = bind_rows(twins,MANUAL) |>
              add_count(name='n') |>
              mutate( is_matched = (n1 == 1) ) |>
              group_by(s1) |>
              add_count(name='n1')
  }

  return(results)
}

# Convert to opposite case (uppercase to lowercase and vice versa)
toggle_case = function(s){
  AZ=concat(LETTERS)
  az=concat(letters)
  svec = str2chr(as.character(s))
  if( all(grepl("[[:upper:]]",svec)) ){
    S = tolower(s) # Faster if all characters are uppercase
  }else if( all(grepl("[[:lower:]]",svec)) ){
    S = toupper(s)  # Faster if all characters are lowercase
  }else{
    # If characters used mixed casefold
    S = sapply(svec,function(chr){
        ifelse( grepl("[[:upper:]]",chr), chartr(AZ,az,chr), chartr(az,AZ,chr) )
    })
  }
  return(concat(S))
}

# Format values ----------------------------------------------------------------

RoundUpToNearest <- function(nb, roundto=1){ # Round up numbers to unit or decimal specified
  if (roundto==1){ return(nb) }
  else{ return( ceiling(nb / roundto) * roundto ) }
}

RoundDownToNearest <- function(nb, roundto=1){ # Round down numbers to unit or decimal specified
  if (roundto==1){ return (nb) }
  else{ return( floor(nb / roundto) * roundto ) }
}


percent <- function(x, digits = 2, format = "f", ...) {
  if(all(x>=0) & all(x<=1)){
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }else if(all(x>=0) & all(x<=100)){
    paste0(formatC(x, format = format, digits = digits, ...), "%")
  }
}

percent0 <- function(x, format = "f", ...) {
  percent(x,digits=0,format=format)
}

toUnit <- function(x, digits = 0, format = "f", units,...) {
  if( missing(units) ){ print("You need to specify a string for the units") }
  paste0(formatC(x, format = format, digits = digits, ...), as.character(units) )
}

rank0 <- function(x) {
  sapply( x, function(d){
    d.chr = as.character(d); ndigits = nchar(d.chr);
    last.digit = as.numeric(substr(d.chr,ndigits,ndigits))
    twolast.digit = as.numeric(substr(d.chr,ndigits-1,ndigits))
    if( !(last.digit %in% c(1,2,3)) | twolast.digit %in% c(11,12,13) ){ return(paste0(d,"th")) }
    if( last.digit == 1 & twolast.digit != 11){ return(paste0(d,"st")) }
    if( last.digit == 2 & twolast.digit != 12){ return(paste0(d,"nd")) }
    if( last.digit == 3 & twolast.digit != 13){ return(paste0(d,"rd")) }
  }
  )
}

toUnits <- function(x){
  if( x >= 1e+24 ){ return(c(pre="Y",x=format(x/1e+24,digit=1)))       }
  else if( x >= 1e+21  &  x < 1e+24 ){ return(c(pre="Z",x=format(x/1e21,digit=1)))   }
  else if( x >= 1e+18  &  x < 1e+21 ){ return(c(pre="E",x=format(x/1e18,digit=1)))   }
  else if( x >= 1e+15  &  x < 1e+18 ){ return(c(pre="P",x=format(x/1e15,digit=1)))   }
  else if( x >= 1e+12  &  x < 1e+15 ){ return(c(pre="T",x=format(x/1e12,digit=1)))   }
  else if( x >=  1e+9  &  x < 1e+12 ){ return(c(pre="G",x=format(x/1e9,digit=1)))    }
  else if( x >=  1e+6  &  x < 1e+9  ){ return(c(pre="M",x=format(x/1e6,digit=1)))    }
  else if( x >=  1e+3  &  x < 1e+6  ){ return(c(pre="k",x=format(x/1e3,digit=1)))    }
  else if( x >=  1e+2  &  x < 1e+3  ){ return(c(pre="h",x=format(x/1e2,digit=1)))    }
  else if( x >=  1e+1  &  x < 1e+2  ){ return(c(pre="da",x=format(x/1e1,digit=1)))   }
  else if( x <   1e+1  &  x > 1e-1  ){ return(c(pre=" ",x=format(x,digit=1)))        }
  else if( x <=  1e-1  &  x > 1e-2  ){ return(c(pre="d",x=format(x/1e-1,digit=1)))   }
  else if( x <=  1e-2  &  x > 1e-3  ){ return(c(pre="c",x=format(x/1e-2,digit=1)))   }
  else if( x <=  1e-3  &  x > 1e-6  ){ return(c(pre="m",x=format(x/1e-3,digit=1)))   }
  else if( x <=  1e-6  &  x > 1e-9  ){ return(c(pre="mu",x=format(x/1e-6,digit=1)))  }
  else if( x <=  1e-9  &  x > 1e-12 ){ return(c(pre="n",x=format(x/1e-9,digit=1)))  }
  else if( x <= 1e-12  &  x > 1e-15 ){ return(c(pre="p",x=format(x/1e-12,digit=1)))  }
  else if( x <= 1e-15  &  x > 1e-18 ){ return(c(pre="f",x=format(x/1e-15,digit=1)))  }
  else if( x <= 1e-18  &  x > 1e-21 ){ return(c(pre="a",x=format(x/1e-18,digit=1)))  }
  else if( x <= 1e-21  &  x > 1e-24 ){ return(c(pre="z",x=format(x/1e-21,digit=1)))  }
  else if( x <= 1e-24 ){ return(c(pre="y",x=format(x/1e-24,digit=1)))  }
}

decompose_int = function(x,desc=F){
  if(is.na(x)){ return(NA) }
  if(!is.numeric(x)){ stop("Input is not a valid number!") }
  x_int = abs(as.integer(x)) # convert to positive integer
  nd = nchar(x_int) # number of digits
  P = seq(1,nd) # position in number
  by10 = rev(as.integer(10^(P-1))) # power of 10
  decomposed = sapply(P, function(p){  as.integer( substring(x_int,p,p) ) * by10[p] })

  if( sum(decomposed) == x_int ){
    return(sort(decomposed,decreasing = desc))
  }else{
    print(sort(decomposed,decreasing = desc))
    stop('Decomposition is wrong!')
  }
}

to_iupac_multiplier = function(x){
  if(!is.numeric(x)){ stop("must give a positive integer!") }
  if(is.na(x)){ return(NA) }

  if( x==0 || x > 9999 ){
    poly = set_names(x,sprintf('(%s)–',x))
    return(poly)
  }
  X = decompose_int(x)
  N = length(X)
  #Multiplier # Number
  # Exceptions for 11,100,1000
  # Exceptions for 12,20,200,2000
  iupac_multi=c(
  "mono-" = 1, "di-" = 2, "tri-" = 3, "tetra-" = 4, "penta-" = 5, "hexa-" = 6,
  "hepta-" = 7, "octa-" = 8, "nona-" = 9, "deca-" = 10,
  "undeca-" = 11, "dodeca-" = 12, "trideca-" = 13, "tetradeca-" = 14,
  "pentadeca-" = 15, "hexadeca-" = 16, "heptadeca-" = 17, "octadeca-" = 18,
  "nonadeca-" = 19, "icosa-" = 20, "henicosa-" = 21, "docosa-" = 22, "tricosa-" = 23,
  "triaconta-" = 30, "hentriaconta-" = 31, "dotriaconta-" = 32,
  "tetraconta-" = 40, "hentetraconta-" = 40, "dotetraconta-" = 40,
  "pentaconta-" = 50, "henpentaconta-" = 50, "dopentaconta-" = 50,
  "hexaconta-" = 60, "henhexaconta-" = 60, "dohexaconta-" = 60,
  "heptaconta-" = 70, "henheptaconta-" = 70, "doheptaconta-" = 70,
  "octaconta-" = 80, "henoctaconta-" = 80, "dooctaconta-" = 80,
  "nonaconta-" = 90, "hennonaconta-" = 90, "dononaconta-" = 90,
  "hecta-" = 100, "dicta-" = 200, "tricta-" = 300, "tetracta-" = 400,
  "pentacta-" = 500, "hexacta-" = 600, "heptacta-" = 700, "octacta-" = 800,
  "nonacta-" = 900,
  "kilia-" = 1000, "dilia-" = 2000, "trilia-" = 3000, "tetralia-" = 4000,
  "pentalia-" = 5000, "hexalia-" = 6000, "heptalia-" = 7000, "octalia-" = 8000,
  "nonalia-" = 9000)

  multipliers = iupac_multi[iupac_multi %in% X]


  if( X[1] %in% c(1,2) & length(X) > 1 ){
    first = iupac_multi[ iupac_multi == sum(X[1:2]) ]
    multipliers = c(first,multipliers[-c(1:2)])
  }
  return(multipliers)
}

to_oligomer = function(nsub){
  if(is.na(nsub)){ return(NA) }
  multis =  to_iupac_multiplier(nsub)
  affixes = paste0(names(multis),collapse='')
  affix = gsub('-','',affixes)
  oligomer = paste0(affix,'mer')
  return(oligomer)
}

# Sequences --------------------------------------------------------------------
load.proteome = function(url,nostop=T) {
  library(Biostrings)
  p = Biostrings::readAAStringSet(filepath = url)
  if(nostop){ # Remove the trailing star from amino acid sequence (stop codon)
    star = Biostrings::subseq(p,start=-1) == '*'
    p = Biostrings::subseq(p,start=1,  end=Biostrings::width(p)-star)
  }
  return(p)
}

load.genome = function(url) {
  library(Biostrings)
  g = Biostrings::readDNAStringSet(filepath = url)
  return(g)
}
SGD.nomenclature = function(){ return("([Y][A-P][LR][0-9]{3}[WC](?:-[A-Z])?)|(Q[0-9]{4})|(R[0-9]{4}[WC])") }
get.orf = function(x){ stringr::str_extract(x,SGD.nomenclature()) }
# Amino acids ------------------------------------------------------------------
# SEQ must be a character vector
is.pos  <- function(SEQ){ SEQ == "K" | SEQ == "R" }
is.neg  <- function(SEQ){ SEQ == "D" | SEQ == "E" }
pos     <- function(SEQ){ sum(is.pos(SEQ)) }
neg     <- function(SEQ){ sum(is.neg(SEQ)) }
fpos    <- function(SEQ){ mean(is.pos(SEQ)) }
fneg    <- function(SEQ){ mean(is.neg(SEQ)) }
charged <- function(SEQ){ return( pos(SEQ) + neg(SEQ) ) }
fcr     <- function(SEQ){ return( fpos(SEQ) + fneg(SEQ)  ) }
npcr    <- function(SEQ){ return( abs( fpos(SEQ) - fneg(SEQ) )  ) }
charge.asym <- function(p,n){ #p=positives n=negatives (sums)
  if( (p+n) == 0){ return(0) }
  return( (p-n)**2 / (p+n) )
}

# Fasta format
seq2fasta <- function(seq, maxline=70, as.string=F){
  single = unlist(strsplit(as.character(seq),split=''))
  string = paste(single,collapse='')
  L = nchar(string)
  fasta.nrows = (L %/% maxline)

  fasta.seq = c()
  for(i in 0:(fasta.nrows)){
    start=(i*maxline) +1
    fasta.seq= c( fasta.seq, substr(string,start,start+70) )
  }
  if(as.string){ return( paste0(fasta.seq,collapse='') ) }
  return( paste0(fasta.seq,collapse='\n') )
}

# Graphics/Colors --------------------------------------------------------------
set.contrast.text <- function(color){
  ifelse( mean(col2rgb(color)) > 127, "black", "white")
}

# perceived brightness
#brightness = function(r, g, b){  return (r * 299 + g * 587 + b * 114) / 1000 }
#brighterThan = function(r, g, b, x) { return (brightness(r,g,b) > x); }

# Desaturate colors by specified proportion
desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

nice_par = function(mar = c(3, 3, 2, 1), mgp = c(2, 0.4, 0), tck = -0.01,
                    cex.axis = 0.9, las = 1, mfrow = c(1, 1), ...) {
  par(mar = mar, mgp = mgp, tck = tck, cex.axis = cex.axis, las = las,
      mfrow = mfrow, ...)
}

# Environment ------------------------------------------------------------------
script.to.env = function(src,nm){ # load script object into separate environment
  if( !(nm %in% search()) ){
    env = attach(what = NULL,name = nm);
  }else{
    env = as.environment(nm)
  }
  sys.source(file=src,env)
}

detachAllPackages <- function() {
  basic.packages <- paste0("package:",c("stats","graphics","grDevices","utils","datasets","methods","base"))
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

showMemoryUse <- function(sort="size", decreasing=FALSE, limit) {
  objectList <- ls(parent.frame())
  oneKB <- 1024
  oneMB <- 1048576
  oneGB <- 1073741824
  memoryUse <- sapply(objectList, function(x) as.numeric(object.size(eval(parse(text=x)))))
  memListing <- sapply(memoryUse, function(size) {
    if (size >= oneGB) return(paste(round(size/oneGB,2), "GB"))
    else if (size >= oneMB) return(paste(round(size/oneMB,2), "MB"))
    else if (size >= oneKB) return(paste(round(size/oneKB,2), "kB"))
    else return(paste(size, "bytes"))
  })
  memListing <- data.frame(objectName=names(memListing),memorySize=memListing,row.names=NULL)
  if (sort=="alphabetical") memListing <- memListing[order(memListing$objectName,decreasing=decreasing),]
  else memListing <- memListing[order(memoryUse,decreasing=decreasing),] #will run if sort not specified or "size"
  if(!missing(limit)) memListing <- memListing[1:limit,]
  print(memListing, row.names=FALSE)
  return(invisible(memListing))
}

ls.objects <- function (pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

quiet <- function(x) {
  # Make a function quiet
  sink("/dev/null")
  on.exit(sink())
  invisible(force(x))
}

get_os <- function(){
  # Check OS type (distinguish MacOS and Linux)
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

get_source_path <- function(){
  if(.Platform$GUI == "RStudio"){
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
}


# shorthands -------------------------------------------------------------------
replace_name <- function(.x, ..., .strict = TRUE) {
  pos <- tidyselect::eval_rename(quote(c(...)), .x, strict = .strict)
  names(.x)[pos] <- names(pos)
  .x
}

lf <- function(path = ".", maxdepth = 0L, pattern = NULL, all.files = FALSE, include.dirs = FALSE) {
  # https://stackoverflow.com/questions/70957718/how-to-list-files-until-n-level-deep-subdirectory-in-r
  # example: lf(".", maxdepth = n, pattern = "[.]xlsx$")
    fn <- list.files(path, pattern = pattern, all.files = all.files, full.names = TRUE, no.. = TRUE)
    dn <- list.dirs(path, full.names = TRUE, recursive = FALSE)
    fn <- fn[match(fn, dn, 0L) == 0L]
    if (!all.files) {
      dn <- dn[grepl("^[^.]", basename(dn))]
    }
    if (length(dn) == 0L) {
      return(fn)
    }
    if (maxdepth < 1L) {
      if (include.dirs) {
        return(c(fn, dn))
      } else {
        return(fn)
      }
    }
    l <- lapply(dn, lf, maxdepth = maxdepth - 1L, pattern = pattern, all.files = all.files, include.dirs = include.dirs)
    if (include.dirs) {
      l <- Map(c, dn, l, USE.NAMES = FALSE)
    }
    c(fn, unlist(l, FALSE, FALSE))
}

find.files = function(directory, pattern, full=T){
# shorthand for non-recursive list.files() (depth=0)
  list.files(directory, pattern=pattern, full.names=full, recursive=F)
}

find.fasta = function(directory,full=T){
  find.files(directory=directory,pattern="\\.(fa|fas|fasta)$",full=full)
}


get_rows_by_keyword = function(word,df){
  library(tidyverse)
  rows = df |> dplyr::filter(dplyr::if_any(everything(),stringr::str_detect, paste0("(?i)",as.character(word))))
  return(rows)
}

find_keywords = function(df,keywords,strict=T){
  library(tidyverse)

  count_keywords = list()
  for( kw in  unique(keywords) ){
    count_kw = sapply(df,
                function(x){ str_count( string=tolower(x), pattern=paste0("(?i)",as.character(kw))) }) |> rowSums()
    count_keywords[[kw]]=count_kw
  }

  df_count_keywords = df |> bind_cols(count_keywords) |>
    rowwise() |>
    mutate(keyword_matched=sum(c_across(all_of(keywords)))) |>
    filter(keyword_matched>0) |>
    arrange(desc(keyword_matched)) |> ungroup()

  if(strict){
   MAX_K = max(df_count_keywords$keyword_matched)
   return( df_count_keywords |> dplyr::filter( keyword_matched == MAX_K) )
 }
  return( df_count_keywords )
}

get_match = function(x,y){ x[match(x,y)] }

load.package <- function(name) { # Quietly Load Libraries
  if(length(name)>1){
    for(i in 1:length(name)){ load.package(name[i]) }
  }else if(length(name)==1){
    suppressMessages(suppressWarnings(library(name, quietly = T, warn.conflicts = F, character.only = T)))
  }else{
    warning('no package to load...')
  }
}

lsos <- function(..., n=10) { # checks top10 memory consumption from R objects
  ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

# initialize empty dataframe with columns names
make_dataframe = function(col_names,n_rows=0,def_val=NA){
  nc = length(col_names)
  df = setNames( data.frame(matrix(data = def_val,nrow = n_rows, ncol=nc)), nm = col_names)
  return(df)
}

ht = function(d, n = 6){ rbind(head(d, n), tail(d, n)) }  # ht == headtail
hh = function(d){ d[1:5, 1:5] } # show the first 5 rows & first 5 columns of a data frame

paste.paired <- function(x,s='-'){ # join by pair (last value remains alone if odd number of values)
  if( is.even(length(x)) ){ return(paste.even(x,s)) }
  return(paste.odd(x,s))
}

paste.pair = function(x,s='_'){ # join by pair (upper limit of pair i is the lower_limit of pair i+1)
  pair=c()
  for (i in seq_along(x[-1]) ){
    pair[i] = paste(x[i],x[i+1],sep=s)
  }
  return(pair)
}

between_ <- function(x,l,r){
  if(any(is.na(c(l,r)))){ return(NA) }
  x>min(l,r) & x<max(l,r)
}

get.even  <- function(x){ x[is.even(x)] } # returns even values
get.odd   <- function(x){ x[is.odd(x)]  } # returns odd values
get.dup   <- function(x){ x[is.dup(x)]  } # returns duplicates
get.mid   <- function(x){ x[mid(x)] }     # returns the middle value
get.half1 <- function(x){ x[seq(1,mid(x))] }  # returns first half of vector
get.half2 <- function(x){ x[seq(mid(x)+1,length(x))] }   # returns last half of vector

round2near <- function(nb,roundto){ # returns the nearest number up to the specified decimal
  a = RoundDownToNearest(nb,roundto);
  b = RoundUpToNearest(nb,roundto);
  if( abs(nb-a) < abs(nb-b) ){ return(a) }
  else if( abs(nb-a) > abs(nb-b) ){ return(b) }
  else{ return(nb) }
}
Round2Nearest <- function(...){ round2near(...) }

table_ = function(x,...){ addmargins(table(x),...) }
mad_      <- function(...){ mad(na.rm=T,...) } # Median absoluted deviations without NA
quantile_ <- function(x,...){  quantile(x,...,na.rm=T) } # quantile with no error message for missing values
range_ <- function(x,...){  range(x,...,na.rm=T) } # range with no error message for missing values


slope <- function(x, y){
  mean_x <- mean(x,na.rm=T)
  mean_y <- mean(y,na.rm=T)
  nom <- sum((x - mean_x)*(y-mean_y),na.rm=T)
  denom <- sum((x - mean_x)^2,na.rm=T)
  m <- nom / denom
  return(m)
}

# the slope formula is just
# covariance(x, y) / variance(x)
slope2 <- function(x, y){
  return(cov(x, y,use = 'pairwise', 'pearson')/var(x,na.rm = T))
}

intercept <- function(x, y, m){
  b <- mean(y,na.rm=T) - (m * mean(x,na.rm=T))
  return(b)
}
# Correlation with spearman method (ranks) and pairwise value
scor <- function(x,y,met='spearman',use='pairwise.complete.obs'){
  res = cor.test(x,y,method = met, use=use, exact=F) |>
        replace_name("r"="estimate","p"="p.value")
  res$pv = res$p.value
  res$p.value = ifelse(res$pv==0,"<1e-324" ,sprintf("%.1e",res$pv))
  return(res)
}

pearson <- function(X,Y){
  library(broom)
  res = cor.test(x = X, y=Y , method = "pearson",use='pairwise.complete',exact=F) |>
    broom::tidy() |>
    mutate( pv = p.value,
            p.value= ifelse(pv==0,"<1e-324" ,sprintf("%.1e",pv))) |>
    dplyr::rename(r=estimate,p=p.value)
  return(res)
}

spearman <- function(X,Y){
  library(broom)
  res = cor.test(x = X, y=Y , method = "spearman",use='pairwise.complete',exact=F) |>
        broom::tidy() |>
    mutate( pv = p.value,
            p.value= ifelse(pv==0,"<1e-324" ,sprintf("%.1e",pv))) |>
    dplyr::rename(r=estimate,p=p.value)
  return(res)
}

# Get spearman correlation parameters ready to plot
spearman.toplot = function(X,Y,rm.slope=T){
  s   = spearman(X,Y)
  s$slope = slope(X,Y)
  s$N = sum(complete.cases(X,Y))
  s$toshow = sprintf("%.3f\n%.3f\n%4s\n%4s",s$slope,s$r,s$p,s$N)
  if(rm.slope){
    s$toshow = sprintf("%.3f\n%4s\n%4s",s$r,s$p,s$N)
  }
  s$xmax = max_(X)
  s$ymax = max_(Y)
  s$xmin = min_(X)
  s$ymin = min_(Y)
  return(s)
}

pearson.toplot = function(X,Y,rm.slope=T){
  p   = pearson(X,Y)
  p$slope = slope(X,Y)
  p$N = sum(complete.cases(X,Y))
  p$toshow = sprintf("%.3f\n%.3f\n%4s\n%4s",p$slope,p$r,p$p,p$N)
  if(rm.slope){
    p$toshow = sprintf("%.3f\n%4s\n%4s",p$r,p$p,p$N)
  }
  p$xmax = max_(X)
  p$ymax = max_(Y)
  p$xmin = min_(X)
  p$ymin = min_(Y)
  return(p)
}

# Extract correlation parameters as dataframe
get.cor.param = function(x,y,...){ as.data.frame(scor(x,y,...)[c('estimate','pv','p.value')]) }

# precomputed data -------------------------------------------------------------
get.AA1 = function(){ unlist(strsplit("ACDEFGHIKLMNPQRSTVWY","")) }
#aaa=c('ALA',"TYR","K","ser","Thr")

regex_AA1 = function(extra=c('\\*','\\-'),negate=F,any=T){
# Amino acid scores used in Dubreuil et al. 2019
  bracket1="["
  bracket2="]"
  if(negate){ bracket1="[^" }
  if(any){ bracket2="]+" }
  aas = c(bracket1,get.AA1(),extra,bracket2)
  return(concat(aas))
}

# Checks wether character corresponds to 1-letter amino acid code
is.aa1    <- function(x,.add.names=T,.ignore.case=T){
  return( is.in(x,get.AA1(), case.sensitive=!.ignore.case, withNames=.add.names) )
}

# Checks wether character corresponds to 3-letter amino acid code
is.aa3    <- function(x,.add.names=T,.ignore.case=T){
  return( is.in(x,get.AA3(), case.sensitive=!.ignore.case, withNames=.add.names) )
}

get.AAA = function(){ c("ALA","CYS","ASP","GLU","PHE",
                        "GLY","HIS","ILE","LYS","LEU",
                        "MET","ASN","PRO","GLN","ARG",
                        "SER","THR","VAL","TRP","TYR") }

get.AA3 = function(){ c('A'="ALA",'C'="CYS",'D'="ASP",'E'="GLU",'F'="PHE",
                        'G'="GLY",'H'="HIS",'I'="ILE",'K'="LYS",'L'="LEU",
                        'M'="MET",'N'="ASN",'P'="PRO",'Q'="GLN",'R'="ARG",
                        'S'="SER",'T'="THR",'V'="VAL",'W'="TRP",'Y'="TYR") }

get.AA.df = function(){
  aaa = get.AA3()
  return( data.frame( a=names(aaa), aaa ) )
}

# Amino acid properties based on seqinr::SEQINR.UTIL$AA.PROPERTY
get.aa.poperties = function(){
  aa_prop = list(
    tiny      = c('A','C','G','S','T'),
    small     = c('A','B','C','D','G','N','P','S','T','V'),
    aliphatic = c('I','L','V'),
    aromatic  = c('F','H','W','Y'),
    nonpolar  = c('A','C','F','G','I','L','M','P','V','W','Y'),
    polar     = c('D','E','H','K','N','Q','R','S','T','Z'),
    charged   = c('B','D','E','H','K','R','Z'),
    basic     = c('H','K','R'),
    acidic    = c('B','D','E','Z'),
    alcohol   = c('S','T'),
    turnlike  = c('A','C','D','E','G','H','K','N','Q','R','S','T')
  )

  aa_grouped = sapply(aa_prop,paste0,collapse="")
  names(aa_prop) = paste0(names(aa_prop),"_",aa_grouped)
  return(aa_prop)
}

# Amino acid scores used in Dubreuil et al. 2019
get.aggrescan = function(){
  setNames(
    object=c(-0.036, 0.604, -1.836, -1.412, 1.754, -0.535, -1.033, 1.822, -0.931, 1.38,
        0.91, -1.302, -0.334, -1.231, -1.24, -0.294, -0.159, 1.594, 1.037, 1.159),
    nm=get.AA1()
  )
}

get.camsol = function(){
  setNames(
    object=c(-0.4533509211, -3.039164569,   3.806164055,  4.112172118, -3.922859742,
         0.809691392,  -0.2864452149, -2.986464178,  3.64837007,  -2.150991033,
        -1.649718164,   1.353162125,   1.405662227,  0.2038292233, 4.093823461,
         0.5372203047, -0.9718594221, -3.593600781, -3.759140236, -2.931244491),
    nm = get.AA1()
  )
}

get.foldamyloid = function(){
  setNames(
    object = c(0.086, 0.568, -0.776, -0.632, 0.958, -1.088, 0.025, 1.217, -0.565, 1.015, 0.725,
  -0.713, -2.303, -0.271, 0.032, -0.73, -0.349, 0.92, 1.027, 0.851),
    nm = get.AA1()
  )
}

get.kytedoolittle = function(){
  setNames(
    object=c(1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5,
  -4.5, -0.8, -0.7, 4.2, -0.9, -1.3),
    nm=get.AA1()
  )
}

get.pawar_ph7 = function(){
  setNames(
    object=c(-3.31, 1.61, -9.42, -10.38, 2.8, -3.96, -4.31, 0.93, -9.55, -0.25, -1.06, -6.02,
  -11.96, -6, -11.93, -5.08, -2.12, 0.49, 2.92, 1.03),
    nm=get.AA1()
  )
}

get.roseman = function(){
  setNames(
    object=c(0.39, 0.25, -3.81, -2.91, 2.27, 0, -0.64, 1.82, -2.77, 1.82,  0.96, -1.91, 0.99,
  -1.3, -3.95, -1.24, -1, 1.3, 2.13, 1.47),
    nm=get.AA1()
  )
}

get.stickiness = function(){
  setNames(
    object=c(0.0062, 1.0372, -0.7485, -0.7893, 1.2727, -0.1771, 0.1204, 1.1109, -1.1806,
  0.9138, 1.0124, -0.2693, -0.1799, -0.4114, -0.0876, 0.1376, 0.1031, 0.7599,
  0.7925, 0.8806),
    nm=get.AA1()
  )
}

get.voronoi_stickiness = function(){
  setNames(
    object=c(-0.040638565,  0.739811788, -0.281053408, -0.392609711, 0.855804199,
      -0.015564908,  0.189364177,  0.593731643, -0.546173776, 0.560032835,
       0.569250951, -0.163722093, -0.089014045, -0.142636678, 0.027334912,
      -0.101035135, -0.024028749,  0.372857295,  0.913761361, 0.79208657),
    nm=get.AA1()
  )

}

get.wimleywhite = function(){
  setNames(
     object=c(4.08, 4.49, 3.02, 2.23, 5.38, 4.24, 4.08, 4.52, 3.77, 4.81, 4.48, 3.83, 3.8,
     3.67, 3.91, 4.12, 4.11, 4.18, 6.1, 5.19),
     nm=get.AA1()
  )
}

sticky = get.stickiness()

# Calculate the delta score between amino acids
get.score.mutation = function(wt,mut,aascore=get.stickiness(), na.score=0){
  if( length(wt) != length(mut)){
    stop("Inputs have different length.")
  }
  validAA = all(c(wt,mut) %in% get.AA1())
  nonAA = setdiff(c(wt,mut),get.AA1())

  if( !validAA ){ warning(sprintf("Some amino acids are not recognized. (%s)",toString(nonAA))) }

  # USE NORMALIZED SCORE (ONLY POSITIVE VALUES)
  # scale(aascore,center=min(aascore),scale=1)
  aascore.norm = normin(aascore)

  # REPLACE NA VALUES BY FIXED SCORE
  wt.score= aascore.norm[wt]
  wt.score[is.na(wt.score)] = na.score

  mut.score= aascore.norm[mut]
  mut.score[is.na(mut.score)] = na.score

  deltascore = (mut.score - wt.score)

  names(deltascore) = paste0(wt,"->",mut)

  return(deltascore)
}

