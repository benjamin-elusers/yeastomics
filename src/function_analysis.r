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
  ft_string = STRING %>%
              dplyr::filter(experiments>100) %>%
              dplyr::select(ORF1,ORF2)
  fromTo=ft_string
  if(missing(fromTo)){ stop("A from-to matrix is required to build the network!") }
  library(tidygraph)
  library(igraph)
  library(netrankr)

  fullnet = as_tbl_graph(fromTo)
  ranknet = fullnet %>% activate(nodes) %>%
        mutate(
          # FAILED (out of memory issue)
          #cen.alpha                       = centrality_alpha(),
          #cen.power                       = centrality_power(),
          # IRRELEVANT
          #cen.manual                      = centrality_manual(),
          #cen.communicability_odd         = centrality_communicability_odd(),
          #cen.communicability_even        = centrality_communicability_even(),
          #cen.subgraph_odd                = centrality_subgraph_odd(),
          #cen.subgraph_even               = centrality_subgraph_even(),
          # CENTRALITY MEASURES TO USE
          cen.degree                      = centrality_degree(),      # degree = regulatory relevance (signaling hubs)
          cen.betweenness                 = centrality_betweenness(cutoff = 20), # number of shortest path going through a node
          cen.closeness                   = centrality_closeness(cutoff = 20), # inverse average length of shortest path
          cen.eigen                       = centrality_eigen(),       # node centrality index (important proteins interact with other important proteins)
          cen.authority                   = centrality_authority(),   # authority = receive from hubs (authority and hub can overlap)
          cen.hub                         = centrality_hub(),         # hub = connected to authorities (authority and hub can overlap)
          cen.pagerank                    = centrality_pagerank(),    # similar to eigen centrality
          cen.subgraph                    = centrality_subgraph(),    # number of subgraphs connected (closed loops from this node)
          cen.integration                 = centrality_integration(), # walk counts
          cen.communicability             = centrality_communicability(),
          cen.katz                        = centrality_katz()
          #cen.betweenness_network         = centrality_betweenness_network(),
          #cen.betweenness_current         = centrality_betweenness_current(),
          #cen.betweenness_communicability = centrality_betweenness_communicability(),
          #cen.betweenness_rsp_simple      = centrality_betweenness_rsp_simple(), # SP = capability of a protein to bring in communication distant proteins
          #cen.betweenness_rsp_net         = centrality_betweenness_rsp_net(),
          #cen.information                 = centrality_information(),
          #cen.decay                       = centrality_decay(),
          #cen.random_walk                 = centrality_random_walk(),
          #cen.expected                    = centrality_expected()
        )


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
fivebins = function(x, applyto=NULL,
                    binnames=c('lowest','lower','mid','high','highest'),
                    b1=get.d10,
                    b2=get.q25,
                    b3=get.iqr,
                    b4=get.q75,
                    b5=get.d90){
  BINS = list(B1 = b1(x), B2  = b2(x), B3  = b3(x,lower=0.4,upper=0.6), B4  = b4(x), B5  = b5(x) )
  names(BINS) = binnames

  if( is.vector(applyto) && length(applyto) == length(x) ){
      BINS = lapply(BINS,function(B){ applyto[ B ]})
  }
  return(BINS)
}

# TESTING
# fivebins(x = xx) %>% lapply(FUN = function(B){ xx[B] } )
# df = data.frame(id=paste0('id',0:100), num=0:100, stringsAsFactors = F)
# fivebins(x = df$num,applyto=df$id)
