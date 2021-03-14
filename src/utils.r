# author: Benjamin Dubreuil
# created: 20/12/2020
# comments: general purpose functions
# ```{r}Sys.time()```
#

# I/O --------------------------------------------------------------------------

# Retrieve text content from URL
open.url <- function(file_url) {
  con <- gzcon(url(file_url))
  txt <- readLines(con)
  closeAllConnections()
  return(textConnection(txt))
}

# Execute loading call and save data unless saved data already exists
preload = function(saved.file,loading.call,doing='create data...'){
  library(tictoc)
  if( !file.exists(saved.file) ){
    cat(doing)
    tic(doing)
    res = eval(substitute(loading.call))
    saveRDS(res,saved.file)
    toc()
  }else{
    res = readRDS(saved.file)
  }
  return(res)
}

# Testing/Subsetting -----------------------------------------------------------
is.whole  <- function(x){ all(floor(x) == x) }      # checks if a value has decimal part (not necessarily integer e.g. 1.0 is whole)
is.string <- function(x){ is.character(x) && length(x) == 1 } # cheks if a value is a single string of letters
is.binary <- function(x){ all(1*x %in% 0:1) }       # checks if a value is 0/1 (binary/logical)
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
normin  <- function(x){ return( scale(x,center=min(x),scale=1)[,1] ) } # adds minimum to values (keeps original range)
norm01  <- function(x){ return( scale(x,center=min(x),scale=diff(range(x)))[,1] ) } # scales values between 0 and 1
maxrep  <- function(x,n){ sort(table(x), decreasing = TRUE)[1:n] } # gets the Nth most repeated values
sum.na  <- function (x,notNA=FALSE){ sum(is.na(x) == !notNA) } # returns sum of NA or not-NA values
geomean <- function(x) {  exp(mean(log(x[x != 0 & !is.na(x)]))) } # returns geometrical mean
geosd   <- function(x) {  exp(sd(log(x[x != 0 & !is.na(x)]))) } # returns geometrical standard deviation

rowsNA = function(m){ # counts how many rows have NAs
  nas = apply(as.matrix(m),1,sum.na)
  if(sum(nas)>0){ message(sprintf('%s rows contain NA values.',sum(nas))) }
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
concat   <- function(x){ return( paste(x,collapse='') ) }        # concatenate characters to string

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
  b4sep = sprintf("([^%s]+).+",sep)
  part1 = sub(b4sep, "\\1", x=name)
  if(lc){ tolower(part1) }
  return(part1)
}

strfind = function(strings, patterns){ # search multiple patterns in character vectors
  sapply(patterns,  function(p){ grep(x = strings, pattern = p, value = T) })
}

file_ext <- function (fn) {
  # remove a path
  splitted    <- strsplit(x=fn, split='/')[[1]]
  # or use .Platform$file.sep in stead of '/'
  fn          <- splitted [length(splitted)]
  ext         <- ''
  splitted    <- strsplit(x=fn, split='\\.')[[1]]
  l           <-length (splitted)
  if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l]
  # the extention must be the suffix of a non-empty name
  ext
}

# Format values ----------------------------------------------------------------

RoundUpToNearest <- function(nb, roundto=1){ # Round up numbers to unit or decimal specified
  if (roundto){ return(nb) }
  else{ return( ceiling(nb / roundto) * roundto ) }
}

RoundDownToNearest <- function(nb, roundto=1){ # Round down numbers to unit or decimal specified
  if (roundto){ return (nb) }
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
charge.asym <- function(p,n){ return( (p-n)**2 / (p+n) ) } #p=positives n=negatives (sums)

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

# Colors -----------------------------------------------------------------------
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

# shorthands -------------------------------------------------------------------
lsos <- function(..., n=10) { # checks top10 memory consumption from R objects
  ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

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

# Correlation with spearman method (ranks) and pairwise value
scor <- function(x,y,met='spearman',use='pairwise.complete.obs'){ cor.test(x,y,method = met, use=use, exact=F) }
spearman <- function(X,Y){
  library(broom)
  res = cor.test(x = X, y=Y , method = "spearman",use='pairwise.complete',exact=F) %>% broom::tidy()
  return(res)
}

# Get spearman correlation parameters ready to plot
spearman.toplot = function(X,Y){
  s   = spearman(X,Y)
  s$N = sum(complete.cases(X,Y))
  s$toshow = sprintf(" r %.3f \n p %.1e \n N %s \n",s$estimate,s$p.value,s$N)
  s$xmax = max_(X)
  s$ymax = max_(Y)
  s$xmin = min_(X)
  s$ymin = min_(Y)
  return(s)
}
# Extract correlation parameters as dataframe
get.cor.param = function(x,y,...){ as.data.frame(scor(x,y,...)[c('estimate','p.value')]) }

# precomputed data -------------------------------------------------------------
# Amino acid scores used in Dubreuil et al. 2019
get.AA1 = function(){ unlist(strsplit("ACDEFGHIKLMNPQRSTVWY","")) }


#aaa=c('ALA',"TYR","K","ser","Thr")

# Checks wether character corresponds to 1-letter amino acid code
is.aa1    <- function(x,.add.names=T,.ignore.case=T){
  return( is.in(x,get.AA1(),
                case.sensitive=!.ignore.case,
                withNames=.add.names)
        )
}

# Checks wether character corresponds to 3-letter amino acid code
is.aa3    <- function(x,.add.names=T,.ignore.case=T){
  return( is.in(x,get.AA3(),
                case.sensitive=!.ignore.case,
                withNames=.add.names)
  )
}

get.AA3 = function(){ c('A'="ALA",'C'="CYS",'D'="ASP",'E'="GLU",'F'="PHE",
                        'G'="GLY",'H'="HIS",'I'="ILE",'K'="LYS",'L'="LEU",
                        'M'="MET",'N'="ASN",'P'="PRO",'Q'="GLN",'R'="ARG",
                        'S'="SER",'T'="THR",'V'="VAL",'W'="TRP",'Y'="TYR") }

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

