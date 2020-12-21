make.bins <- function(tobin, nbin = 5, mode=c('equals','distrib'),
                      force0=T,
                      lowest=50, highest=50){
  # Generate bins of a vector based on quantile or with equal number of observations
  # First and last bin can be adjusted manually
  require(hablar)
  require(dplyr)
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
  message(sprintf("Lowest  bin: p=%02.0f%% n=%.0f (value = %.1f)",Plow*100,length(low),lowval))
  high = sort(tobin,dec=T)[1:idx_high]
  highval = tail(high,n=1)
  message(sprintf("Highest  bin: p=%02.0f%% n=%.0f (value = %.1f)",Phigh*100,length(high),highval))

  mid = tobin[ between(tobin,lowval,highval) ]
  message(sprintf("Remaining values to bin: %s",length(mid)))

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
  message(sprintf("Number of groups: %s",nby))
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
    return( CC %>% filter( !is.na(!!sym(BY)) ) )
  }
  return(CC)
}
