make.bins <- function(tobin, nbin = 5, mode=c('equals','distrib'),
                      force0=T,
                      lowest=50, highest=50){
  # Generate bins of a vector based on quantile or with equal number of observations
  # First and last bin can be adjusted manually
  library(hablar)
  library(dplyr)
  N = length(tobin)
  # Define the lowest/highest bin according to number or fraction of values
  if( between(lowest,0,1) ){
    Plow = lowest
    idx_low = round(lowest*N,0)
  }else if( lowest > 1 ){
    Plow = (lowest/N)
    idx_low = lowest
  }

  if( between(highest,0,1) ){
    Phigh = highest
    idx_high = round(highest*N,0)+1
  }else if( highest > 1 ){
    Phigh = (highest/N)
    idx_high = highest+1
  }

  # Extract lowest/highest bin and find the breakpoint value for each
  low = sort(tobin)[1:idx_low]
  lowval = tail(low,n=1)
  message(sprintf("Lowest  bin: p=%02.0f%% n=%.0f (value = %.1f)\n",Plow*100,length(low),lowval))
  high = sort(tobin,dec=T)[1:idx_high]
  highval = tail(high,n=1)
  message(sprintf("Highest  bin: p=%02.0f%% n=%.0f (value = %.1f)\n",Phigh*100,length(high),highval))

  mid = tobin[ between(tobin,lowval,highval) ]
  message(sprintf("Remaining values to bin: %s\n",length(mid)))

  # Finding breaking points for middle values
  mode = match.arg(arg = mode, choices = c("equals", "distrib"))
  if (mode == "equals") {
    midbreaks = quantile(mid, probs = seq(0, 1, len = (nbin-2+1)), na.rm = TRUE)
  }else if (mode == "distrib") {
    midbreaks = seq(min_(mid), max_(mid), len = (nbin-2+1) )
  }

  # Finalize binning with manual and auto-adjusted bins
  if(force0){ initval = 0 }else{ initval=min_(tobin) }
  breaks=unique(c(initval,lowval,midbreaks,highval,max_(tobin)))
  binned = cut(tobin, breaks, include.lowest = T, labels = paste.pair(sprintf("%.0f",breaks)),ordered_result = T)

  return(binned)
}



# Calculate the spearman correlation of two variables by group
cor.sub.by = function(DATA,  XX, YY, BY, ID=NULL,na.rm=T){
  BYCOL=DATA[[BY]]
  if( !is.factor(BYCOL) ){ BINS = na.exclude(unique(BYCOL)) }
  BINS = levels(BYCOL)
  nby = length(BINS)
  message(sprintf("Number of groups: %s\n",nby))
  bins = sprintf("[%s] %s",seq_along(BINS),as.character(levels(BINS)))
  cat(str_wrap(toString(bins), width = 90))
  cat("\n")
  CC = DATA %>% dplyr::select(XX,YY,BY,ID) %>%
    group_by(!!sym(BY),.drop = T) %>%
    mutate( R = scor(!!sym(XX),!!sym(YY))$estimate, P=scor(!!sym(XX),!!sym(YY))$'p.value', n=n()) %>%
    summarise( r=unique(R) , p=unique(P),
               N=n(),
               na.xy= sum( is.na(!!sym(XX)) | is.na(!!sym(YY))),
               nax=sum(is.na(!!sym(XX))), nay=sum(is.na(!!sym(YY))),
               n=N-na.xy
    )
  if(na.rm){
    message('removing NAs...')
    return( CC %>% dplyr::filter( !is.na(!!sym(BY)) ) )
  }
  return(CC)
}


# Measure centrality of a network
network.centrality = function(fromTo){
  if(missing(fromTo)){ stop("A from-to matrix is required to build the network!") }
  library(tidygraph)
  net %>%
  mutate(cent_alpha = centrality_authority())

  }

# Find the specific quantiles
q25 = function(x){ quantile(x,0.25) }
q75 = function(x){ quantile(x,0.75) }
d10 = function(x){ quantile(x,0.1) }
d40 = function(x){ quantile(x,0.4) }
d60 = function(x){ quantile(x,0.6) }
d90 = function(x){ quantile(x,0.9) }

# Get the values below/above a specified quantile
get.d10 = function(x,with_ties=T){ as.logical( with_ties*(x==d10(x)) + (x < d10(x)) ) }
get.q25 = function(x,with_ties=T){ as.logical( with_ties*(x==q25(x)) + (x < q25(x)) ) }
get.q75 = function(x,with_ties=T){ as.logical( with_ties*(x==q75(x)) + (x > q75(x)) ) }
get.d90 = function(x,with_ties=T){ as.logical( with_ties*(x==d90(x)) + (x > d90(x)) ) }
# xx= 0:100
# xx[ get.d10(xx,with_ties = F) ]
# xx[ get.q25(xx,with_ties = F) ]

# Get values in interquantile range
get.iqr = function(x,lower,upper,with_ties=T){
    if(with_ties){
      x >= quantile(x,lower) & x <= quantile(x,upper)
    }else{
      x > quantile(x,lower) & x < quantile(x,upper)
    }
}

# Slice tibble data for specicic quantiles or interquantile range
slice.d10  = function(x,...){ dplyr::slice_min(x,prop=0.1,with_ties = T,...) }
slice.d90  = function(x,...){ dplyr::slice_max(x,prop=0.1,with_ties = T,...) }
slice.q25  = function(x,...){ dplyr::slice_min(x,prop=0.25,with_ties = T,...) }
slice.q75  = function(x,...){ dplyr::slice_max(x,prop=0.25,with_ties = T,...) }
slice.iqr  = function(.data,x,lower,upper,...){
  col <- enquo(x)
  dplyr::filter(.data,between(!!col,quantile(!!col,lower),quantile(!!col,upper)))
}

# Generate five bins arbitrarily defined at quantiles: 10, 25, 40-60, 75, 90
fivebins = function(x,
                    binnames=c('lowest','lower','mid','high','highest'),
                    b1=get.d10,
                    b2=get.q25,
                    b3=get.iqr,
                    b4=get.q75,
                    b5=get.d90){
  BINS = list(B1 = b1(x), B2  = b2(x), B3  = b3(x,lower=0.4,upper=0.6), B4  = b4(x), B5  = b5(x) )
  names(BINS) = binnames
  return(BINS)
}

# TESTING
# fivebins(x = xx) %>% lapply(FUN = function(B){ xx[B] } )
# df = data.frame(id=paste0('id',0:100), num=0:100, stringsAsFactors = F)
# fivebins(x = df$num) %>% lapply(FUN = function(B){ df$id[B] } )
