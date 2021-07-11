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
  # if( length(BY) == 1){
  #   BYCOL=DATA[[BY]]
  #   if( !is.factor(BYCOL) ){ BINS = na.exclude(unique(BYCOL)) }
  #   BINS = levels(BYCOL)
  #   nby = length(BINS)
  #   message(sprintf("Number of groups: %s\n",nby))
  #   bins = sprintf("[%s] %s",seq_along(BINS),as.character(levels(BINS)))
  #   cat(str_wrap(toString(bins), width = 90))
  #   cat("\n")
  # }

  CC = DATA %>% dplyr::select(XX,YY,BY,ID) %>%
    group_by(across(all_of(BY)),.drop = T) %>%
    drop_na(!!sym(XX),!!sym(YY)) %>%
    mutate( R = scor(!!sym(XX),!!sym(YY))$estimate, P=scor(!!sym(XX),!!sym(YY))$'p.value', n=n()) %>%
    summarise( r=unique(R) , p=unique(P),
               N=n(),
               na.xy= sum( is.na(!!sym(XX)) | is.na(!!sym(YY))),
               nax=sum(is.na(!!sym(XX))), nay=sum(is.na(!!sym(YY))),
               n=N-na.xy
    )
  if(na.rm){
    message('removing NAs...')
    #return( CC %>% dplyr::filter(  !is.na(across(all_of(BY))) )
  }
  return(CC)
}



# Measure centrality of a network
network.centrality = function(fromTo,namenet=''){
  if(missing(fromTo)){ stop("A from-to matrix is required to build the network!") }
  library(igraph)
  library(tictoc)
  #library(tidygraph)
  #library(netrankr)
  tic(" - Build graph...")
  fullnet = graph_from_data_frame(fromTo,directed = F)
  message("Number of nodes: ",length(V(fullnet)))
  comp = components(fullnet)
  message("Number of components: ",comp$no)
  comp.sizes = table(comp$csize)
  cat("Components sizes (# count) :\n")
  cat("----------------------------\n")
  cat(sprintf("\n%-5s (# %-3s)",names(comp.sizes),comp.sizes))
  cat("\n\n")
  toc()

  node.centrality = tibble( ids = as_ids(V(fullnet)),
          cent_deg = igraph::degree(fullnet)
         )

  # Get the largest connected component
  tic(" - Find largest connected component...")
  connet = decompose(fullnet)[[which.max(comp$csize)]]
  nnodes= length(V(connet))
  message("Maximum number of nodes: ",nnodes)
  toc()

  # Compute centrality measures for nodes of largest connected component
  tic(" - Compute centrality measures in connected component...")
  path.centrality = tibble( ids = as_ids(V(connet)) ) %>%
    mutate(
      #cent_alpha = igraph::alpha_centrality(connet,loops = T),
      cent_betweenness = igraph::betweenness(connet,normalized=T),
      cent_closeness = igraph::closeness(connet,normalized = T),
      cent_eccentricity = igraph::eccentricity(connet),
      cent_pagerank = igraph::page_rank(connet)$vector,
      cent_eigen = igraph::eigen_centrality(connet,scale = T)$vector,
      cent_authority = igraph::authority_score(connet,scale = T)$vector,
      cent_hub = igraph::hub.score(connet,scale = T)$vector,
      cent_subgraph= igraph::subgraph_centrality(connet,diag = F)
  )
  toc()

  centrality = left_join(node.centrality,path.centrality)

  # Benchmark correlation for the centrality measures
  tic(" - Find correlation between centrality measures...")
  corrplot::corrplot(cor(centrality[,-1],use = 'complete',method ='spearman'),
                     method='circle', addCoef.col = 'violet',
                     number.cex = 0.7, number.font = 1, number.digits = 2,
                     title=paste0("\n",namenet," network"))
  toc()
  return(centrality)
}

# Find the specific quantiles
q25 = function(x){ quantile(x,0.25) }
q75 = function(x){ quantile(x,0.75) }
d10 = function(x){ quantile(x,0.1) }
d40 = function(x){ quantile(x,0.4) }
d60 = function(x){ quantile(x,0.6) }
d90 = function(x){ quantile(x,0.9) }

# Get the values below/above a specified quantile
get.d10 = function(x,with_ties=T){ (with_ties & (x==d10(x)) + (x < d10(x)) ) }
get.q25 = function(x,with_ties=T){ (with_ties & (x==q25(x)) + (x < q25(x)) ) }
get.q75 = function(x,with_ties=T){ (with_ties & (x==q75(x)) + (x > q75(x)) ) }
get.d90 = function(x,with_ties=T){ (with_ties & (x==d90(x)) + (x > d90(x)) ) }
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


slice_min_max <- function(df, order_by = value, n = 1) {
  order_by = enquo(order_by)
  min <- slice_min(df, !!order_by, n = n) %>%
    mutate(type = "min")
  max <- slice_max(df, !!order_by, n = n) %>%
    mutate(type = "max")
  df <- bind_rows(min, max) %>%
    as_tibble()
  return(df)
}

# TESTING
# fivebins(x = xx) %>% lapply(FUN = function(B){ xx[B] } )
# df = data.frame(id=paste0('id',0:100), num=0:100, stringsAsFactors = F)
# fivebins(x = df$num,applyto=df$id)

# COUNT AMINO ACIDS BY PROTEIN
#dubres = load.dubreuil2019.data(3)
#fullAA=dubres$resname
# APPLY AMINO ACID SCORES COUNT AMINO ACIDS BY PROTEIN
#data(AAindex,package = 'protr')

AACOUNT2SCORE = function(COUNT,SCORE, opposite=F){
  if( is.null(dim(COUNT)) ){ stop("COUNT must be a 2D object") }
  if( is.null(names(SCORE)) ){ stop("SCORE must be a 1D vector with names corresponding to row or col names of COUNT") }
  #nms1 = setdiff(names(SCORE),colnames(COUNT))
  #nms2 = setdiff(names(SCORE),rownames(COUNT))
  if(opposite){
    SCORE = SCORE * (-1)
    print("TAKING OPPOSITE VALUES ")
  }
  if( length(SCORE) == nrow(COUNT) ){ # AA SCORE BY ROW, PROT BY COL
    #print("AMINO ACID BY ROW, PROTEIN BY COLUMN")
    score.m = matrix(as.numeric(SCORE),nrow=nrow(COUNT),ncol=ncol(COUNT),byrow=F)
    rownames(score.m) = rownames(COUNT)
    S = colSums(COUNT*score.m,na.rm = T)/colSums(COUNT,na.rm = T)
  }else if( length(SCORE) == ncol(COUNT) ){ # AA SCORE BY COL, PROT BY ROW
    #print("PROTEIN BY ROW, AMINO ACID BY COLUMN")
    score.m = matrix(data = as.numeric(SCORE),nrow=nrow(COUNT),ncol=ncol(COUNT),byrow=T)
    colnames(score.m) = colnames(COUNT)
    S = rowSums(COUNT*score.m,na.rm = T)/rowSums(COUNT,na.rm = T)
  }else if( !(length(SCORE) %in% dim(COUNT) ) ){
    stop("SCORE must be a vector of length equal to number of col OR number of rows of COUNT")
  }
  return(S)
}


coalesce_join <- function(x, y,
                          by = NULL, suffix = c(".x", ".y"),
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))

  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce,
    1,
    nchar(to_coalesce) - nchar(suffix_used)
  ))

  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]],
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce

  dplyr::bind_cols(joined, coalesced)[cols]
}

nearest=function(x,y,n=1){
  d=abs(x-y)
  ord = order(d)
  return( ord[1:n] )
}
