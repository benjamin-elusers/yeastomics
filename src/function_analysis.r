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
  library(network)
  library(CINNA)
  library(centiserve)
  library(linkcomm)
  library(sna)
  library(tictoc)
  #library(tidygraph)
  #library(netrankr)
  tic(" - Build graph...")
  fullnet = graph_from_data_frame(fromTo,directed = F)
  full = network::network(as_edgelist(fullnet))

  message("Number of nodes: ",length(V(fullnet)))
  comp = igraph::components(fullnet)
  message("Number of components: ",comp$no)
  comp.sizes = table(comp$csize)
  cat("Components sizes (# count) :\n")
  cat("----------------------------\n")
  cat(sprintf("\n%-5s (# %-3s)",names(comp.sizes),comp.sizes))
  cat("\n\n")
  toc()

  message("Compute global network centralities (even if not fully connected)...")
  tic('  - Global centrality...')
  node.centrality = tibble( ids = as_ids(V(fullnet)) )

  tic("      * CINNA centralities...")
  CINNA.centrality = node.centrality %>%
    add_column(
      cent_closeness_dangalchev = CINNA::dangalchev_closeness_centrality(fullnet),
      #cent_wienerindex = CINNA::wiener_index_centrality(fullnet), # returns Inf for not connected graph
      cent_group = CINNA::group_centrality(fullnet),
      cent_harmonic = CINNA::harmonic_centrality(fullnet),
      cent_localbridging = CINNA::local_bridging_centrality(fullnet)
    )
  toc()

  tic("      * sna centralities...")
  sna.centrality = node.centrality %>%
    add_column(
      cent_stress = sna::stresscent(full),
      cent_load = sna::loadcent(full),
      #cent_flowbet = sna::flowbet(full), # quite slow
      #cent_info = sna::infocent(full) # quite slow
    )
  toc()

  tic("     * igraph centralities...")
  igraph.centrality = node.centrality %>%
    add_column(
      cent_deg = igraph::degree(fullnet),
      cent_betweenness = igraph::betweenness(fullnet,normalized=T),
      cent_eccentricity = igraph::eccentricity(fullnet),
      cent_pagerank = igraph::page_rank(fullnet)$vector,
      cent_eigen = igraph::eigen_centrality(fullnet,scale = T)$vector,
      cent_authority = igraph::authority_score(fullnet,scale = T)$vector,
      cent_hub = igraph::hub_score(fullnet,scale = T)$vector,
      #cent_subgraph= igraph::subgraph_centrality(fullnet,diag = F), # quite slow
      cent_coreness = igraph::coreness(fullnet)
    )
  toc()

  tic("      * centiserve centralities...")
  centiserve.centrality = node.centrality %>%
    add_column(
      cent_topcoef = centiserve::topocoefficient(fullnet), # quite slow -- finished in 2mn
      #cent_bottleneck = centiserve::bottleneck(fullnet),  # quite slow
      cent_clustrank  = centiserve::clusterrank(fullnet), # may be incorrect -- error expected
      cent_diffusion = centiserve::diffusion.degree(fullnet),
      cent_dmnc = centiserve::dmnc(fullnet),
      cent_mnc = centiserve::mnc(fullnet),
      cent_geodesic = centiserve::geokpath(fullnet),
      #cent_katz = centiserve::katzcent(fullnet), # quite slow
      #cent_markov = centiserve::markovcent(fullnet), # quite slow
      cent_laplacian = centiserve::laplacian(fullnet), # quite slow -- finished in 1mn
      #cent_epc = centiserve::epc(fullnet), # quite slow
      cent_leverage = centiserve::leverage(fullnet),
      #cent_entropy = centiserve::entropy(fullnet), # quite slow
      #cent_crossclique = centiserve::crossclique(fullnet), # quite slow
      #cent_community = centiserve::communitycent(fullnet), # quite slow -- requires linkcomm package
      cent_closeness_latora =  centiserve::closeness.latora(fullnet,normalized = T),
      cent_closeness_resid  = centiserve::closeness.residual(fullnet),
      #cent_semilocal = centiserve::semilocal(fullnet) # too slow
   )
  toc()

  global_centrality = node.centrality %>%
                      left_join(igraph.centrality) %>%
                      left_join(centiserve.centrality) %>%
                      left_join(CINNA.centrality) %>%
                      left_join(sna.centrality) %>%
                      dplyr::slice(gtools::mixedorder(ids))
  toc()

  # Get the largest connected component
  tic(" - Find largest connected component...")
  connet = decompose(fullnet)[[which.max(comp$csize)]]
  nnodes= length(V(connet))
  message("Maximum number of nodes: ",nnodes)
  toc()

  # Compute centrality measures for nodes of largest connected component
  tic("Compute centrality measures in connected component...")
  connected.centrality = tibble( ids = as_ids(V(connet)) ) %>%
    mutate(
      #cent_alpha = igraph::alpha_centrality(connet,loops = T),
      cent_closeness = igraph::closeness(connet,normalized = T),
      #cent_closeness_current = centiserve::closeness.currentflow(connet), # quite slow
      #cent_closeness_vitality = centiserve::closeness.vitality(connet), # ERROR Subgraph of graph is not strongly connected
      #cent_hubbell = centiserve::hubbell(connet), # ERROR Hubbell index centrality is not solvable for this graph
      cent_avgdist = centiserve::averagedis(connet),
      cent_barycenter = centiserve::barycenter(connet),
      #cent_centroid = centiserve::centroid(connet), # quite slow
      cent_decay = centiserve::decay(connet),
      cent_radiality = centiserve::radiality(connet) # quite slow -- finished in 1mn
  )
  toc()

  centrality = left_join(global_centrality,connected.centrality)
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
q25 = function(x){ quantile(x,0.25,na.rm=T) }
q75 = function(x){ quantile(x,0.75,na.rm=T) }
d10 = function(x){ quantile(x,0.1,na.rm=T) }
d40 = function(x){ quantile(x,0.4,na.rm=T) }
d60 = function(x){ quantile(x,0.6,na.rm=T) }
d90 = function(x){ quantile(x,0.9,na.rm=T) }

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
slice.iqr  = function(.data,x,lower,upper,negate=F,...){
  col <- enquo(x)
  if(!negate){
    dplyr::filter(.data,between(!!col,quantile(!!col,lower),quantile(!!col,upper)))
  }else{
    dplyr::filter(.data,!between(!!col,quantile(!!col,lower),quantile(!!col,upper)))
  }
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
  min <- slice_min(df, !!order_by, n = n) %>% mutate(type = "min")
  max <- slice_max(df, !!order_by, n = n) %>% mutate(type = "max")
  df <- bind_rows(min, max) %>% as_tibble()
  return(df)
}

slice_min_max <- function(df, order_by = value, n = 1) {
  # Find extreme records using quoted column name
  order_by = enquo(order_by)
  min <- slice_min(df, !!order_by, n = n) %>% mutate(type = "min")
  max <- slice_max(df, !!order_by, n = n) %>% mutate(type = "max")
  df <- bind_rows(min, max) %>% as_tibble()
  return(df)
}

get_outliers_index = function(v,nout=10){
  # Get index of extreme (outliers) values from a vector
  ord=order(v,decreasing=F)
  rev_ord=order(v,decreasing=T)
  out = c(ord[1:nout],rev_ord[nout:1])
  return(out)
}

get_outliers = function(v,nout){
  # Get outliers values (extremes) from a vector
  return( v[get_outliers_index(v,nout)] )
}

get_outliers_boundary = function(v,nout){
  val = v[get_outliers_index(v,nout)]
  return(c('min'=head(val,n=1),'out_min'=val[nout], 'out_max'=val[nout+1], max=last(val)))
}



get_extremes = function(df,column,n=3){
  # Find extreme records using string as column names
  extremes = df[ get_outliers_index(df[[column]],nout = n), ]
  return(extremes)
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

find_na_rows = function(df,as.indices=F){
  NA_codon = rowsNA(df) > 0
  if(as.indices){ return(which(NA_codon)) }
  return(df[which(NA_codon),])
}

coalesce_join <- function(x, y,
                          by = NULL, suffix = c(".x", ".y"),
                          join = dplyr::full_join, ...) {
  # Coalesce: To combine two datasets without duplicating rows on NAs
  # Missing values in x will be replaced by defined values from y for the same key
  # And conversely between NAs of y and defined values of x

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

  coalesced <- purrr::map_dfc(set_names(to_coalesce), ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]],
    joined[[paste0(.x, suffix[2])]]
  )) %>% setNames(to_coalesce)

  dplyr::bind_cols(joined, coalesced)[cols]
}

nearest=function(x,y,n=1,value=F){
  if(length(x)>1){
    warning('x must be a scalar (single value). Using the first value of x...')
    x=x[1]
  }
  d=abs(x-y)
  ord = order(d)
  if(value){ return(y[ord[1:n]]) }
  return( ord[1:n] )
}

text2corner = function(where){
  corners=list()
  corners$bottomleft = c(X=-Inf,Y=-Inf)
  corners$topleft = c(X=-Inf,Y=Inf)
  corners$topright = c(X=Inf,Y=Inf)
  corners$bottomright = c(X=Inf,Y=-Inf)
  which=match.arg(where, choices = c('bottomleft','topleft','topright','bottomright'),several.ok = F)
  return(corners[[which]])
}


make_scatterplot  = function(data2plot,xvar,yvar, grp=NULL,
                             pt.size = 0.7, pt.shape=19, pt.alpha=0.8, pt.col='black', pt.fill='black',
                             labx='',laby='',
                             pal='Spectral',paldir=-1,
                             theme2use = NULL,
                             txtcorner='bottomleft'){
  library(ggthemes)
  library(ggplot2)
  theme_set(theme2use)

  print(head(data2plot))
  corner=text2corner(txtcorner)
  rho=spearman.toplot(data2plot[[xvar]],data2plot[[yvar]])
  p = ggplot(data=data2plot,aes_string(x=xvar,y=yvar)) +
      geom_point(size=pt.size, alpha=pt.alpha,shape=pt.shape,col=pt.col,fill=pt.fill)

  if(!is.null(grp)){
    p = ggplot(data=data2plot,aes_string(x=xvar,y=yvar,fill=grp)) +
        geom_point(size=pt.size, alpha=pt.alpha,shape=pt.shape)
  }
  p = p +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour=NA,size=0.25,alpha=0.25)+
    geom_smooth(method = 'lm',col='gray50',se=F,size=1) +
    geom_text(data=rho,aes(x=corner['X'],y=corner['Y'],label=toshow),hjust='inward',vjust='inward',size=5) +
    xlab(labx) +
    ylab(laby) +
    scale_fill_distiller(palette = pal,direction=paldir)
  plot(p)
  return(p)
}

get_clade_residual_evorate = function(cladedata){
  cat("==> Get residuals of clade-specific evolutionary rate from abundance <==\n")
  if(missing(cladedata)){ stop("Run get_clade_data() with the two clades you want to compare...") }

  fit.evo = cladedata %>%
    broom::augment_columns(x=lm(XX~log10(ppm1),data=.)) %>%
    dplyr::rename(.resid.1 = .resid, .fitted.1=.fitted, .se.fit.1=.se.fit, .hat.1=.hat, .sigma.1=.sigma, .cooksd.1=.cooksd,.std.resid.1=.std.resid)  %>%
    broom::augment_columns(x=lm(YY~log10(ppm2),data=.)) %>%
    dplyr::rename(.resid.2 = .resid, .fitted.2=.fitted, .se.fit.2=.se.fit, .hat.2=.hat, .sigma.2=.sigma, .cooksd.2=.cooksd,.std.resid.2=.std.resid) %>%
    mutate( rY=.resid.2, rX=.resid.1)
  return(fit.evo)
}

#### _analysis function ####
show_density = function(input, # Input data
                        var,   # Variable name
                        nsd=2  # Number of standard deviation to show
){

  V = input[[var]]

  # Statistics
  mu = mean_(V) # Y mean
  md = median_(V) # Y median
  s  = sd_(V) # Y standard deviation
  rg = range_(V)
  symmetry = function(x){ return( sort(c(x,x)) * (-1)^(1:length(x)) ) }
  sd_x = mu + s* symmetry(x=1:nsd)
  # Density
  D=density(V,bw=0.1,na.rm=T)
  dy.ind = sapply(sd_x,function(X){ nearest(X,D$x) })
  sd_y  = D$y[dy.ind]

  df_sd =tibble( xx = sd_x, yy = sd_y)

  A = ggplot(input) +
    ggpubr::grids() +
    stat_density(bw = 0.1, aes_string(x=V)) +
    geom_vline(xintercept = mu, col='white', linetype=1,size=1) + # mean
    geom_segment(data=df_sd,aes(x=xx,xend=xx,y=rep(0,nsd*2),yend=yy), col='white', linetype=2,size=1) + # median
    geom_segment(data=df_sd[1:2,],aes(x=xx,xend=xx,y=yy,yend=c(1,1)), col='red', linetype=2,size=0.25) + # standard deviation (outside)
    geom_errorbarh(data=df_sd,aes(xmin=xx[1],xmax=xx[2],y=0.8,height=0.05),col='red',linetype=1,size=1) + # standard deviation
    ylim(0,1) + scale_y_continuous(position = "left") + xlab("Mean Evolutionary Rate") +
    theme( axis.line.x = element_blank(),axis.ticks = element_blank()) + #aspect.ratio = 3/2,
    coord_flip(xlim=rg, ylim=c(1,0))
  return(A)
}

show_sample = function(input, # Input data
                       pop.mean,pop.range,
                       name='property', # column for sample name
                       id='ORF', # column for id of single observation
                       value=c('MF_nucleotide_binding','MF_molecular_function'), # sample names
                       var){

  which.sample = input[[name]] %in% value
  selected = input %>% dplyr::filter(which.sample) %>% dplyr::select(!!sym(name), !!sym(id), !!sym(var))
  mu.sample = mean_(selected[[var]])

  df.mu = group_by(selected,!!sym(name)) %>%
    summarise(n=n(),mu.sample=mean_(!!sym(var)),sd.sample=sd_(!!sym(var)),
              sdmin=mu.sample-sd.sample, sdmax=mu.sample+sd.sample) %>%
    mutate(MU = pop.mean)


  #selected.prop = c('MF_nucleotide_binding','MF_molecular_function','essential_core','essential_dispensable','pangenome_rare','pangenome_cloud')
  #ER.prop = propfit %>% dplyr::filter(property %in% selected.prop) %>% dplyr::select('property','EVO.FULL','.resid.evo')
  library(ggbeeswarm)
  library(see)

  B = ggplot(selected,aes(y=!!sym(var),fill=!!sym(name),col=!!sym(name))) +
    geom_violinhalf(aes_string(x=name),color=NA,alpha=0.9,show.legend = F) +
    geom_beeswarm(aes_string(x=name),size = 2,  shape = 21, stroke = 0, alpha=0.5,groupOnX = T,dodge.width=0) +
    geom_pointrange_borderless(data=df.mu,mapping = aes(x=!!sym(name),y=mu.sample,ymin=sdmin,ymax=sdmax),size=1,fatten=6, position=position_dodge2(width=1)) +
    #geom_crossbar(data=df.mu,mapping = aes(x=!!sym(name),y=mu.sample,ymin=sdmin,ymax=sdmax),alpha=0.1,size=0.5,fatten=0) +
    geom_hline(df.mu,mapping = aes(yintercept = MU), col='black',linetype=1,size=1) + # mean
    geom_text(data=df.mu,mapping = aes(label=n,x=!!sym(name),y=2.5,vjust='inward'),show.legend = F) +
    xlab(name) + scale_y_continuous(name='',limits = pop.range, expand = c(0,0)) + ylab('') +
    #theme(legend.position = 'none',legend.direction = 'vertical',axis.text.x = element_text(angle = 45)) +
    ggpubr::grids()+
    #scale_fill_material_d(palette = "ice") + scale_color_material_d(palette = "ice")
    scale_fill_material_d(palette="full") + scale_color_material_d(palette="full")
  B
  return(B)
}

get_XY_data = function(input,x=X,y=Y,noNA=T){
  res=dplyr::lst(XX=input[[x]], YY=input[[y]],
                 varnames=c('x'=x,'y'=y),
                 n = c('x'=sum(!is.na(XX)), 'y' = sum(!is.na(YY)), 'xy' = sum(complete.cases(YY,XX)) ),
                 mu = c('x'=mean_(XX),'y'=mean_(YY)),
                 md = c('x'=median_(XX),'y'=median_(YY)),
                 var = c('x'=var_(XX),'y'=var_(YY))
  )

  res$df = input %>% ungroup()
  if( noNA ){
    df.noNA = res$df %>% dplyr::filter(complete.cases(!!sym(y),!!sym(x)))
    res$df = df.noNA
    res$YY=res$df[[y]]
    res$XX=res$df[[x]]

  }
  return(res)
}

get_model_params = function(m,x,y){

  nXY=sum( complete.cases(x,y) )
  mu.y=mean_(y)
  m.params = dplyr::lst(
    fit  = m,
    xx = x,
    yy = y,
    yfit = fitted(fit), # Y-Fitted
    yres = residuals(fit), # Y-Residual
    pfit = coef(fit), # Fitted parameters (intercept, PPM)
    TSS = sum( (yy-mu.y)^2 ),
    ESS = sum( (yfit-mu.y)^2 ), # Explained variance
    RSS = TSS-ESS, # Deviance (Unexplained variance)
    RSE = sqrt( (1 / (nXY-2)) * RSS ), # Residual standard error
  )
  return(m.params)
}

get_fit_data = function(d,m){
  mu.y = d$mu['y']
  var.y=d$var['y']
  nXY = d$n['xy']
  xname=d$varnames['x']
  yname=d$varnames['y']
  fit = d$df %>%
    broom::augment_columns(x=m) %>%
    mutate(ESS=sum_( (.fitted-mu.y)^2 ), TSS=sum_( (!!sym(yname)-mu.y)^2 ), RSS=TSS-ESS,
           s2=TSS/(nXY-1), s2.y = var.y, RS=sum(.resid),
           RSE=sqrt( (1 / (nXY-2)) * RSS), AIC = AIC(m), BIC=BIC(m), LL = readr::parse_number(as.character(logLik(m))))
  return(fit)
}


make_linear_fit = function(input,    # Input data
                           x=X, y=Y, # X/Y Variables to fit
                           only.params=T){

  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']
  # model #
  model.name='Linear'
  f=as.formula(paste0(y,"~",x))
  m = lm(formula =f , data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model = model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

make_logistic_fit = function(input,    # Input data
                             y=Y, x=X, # X/Y Variables to fit
                             only.params=T){
  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']

  # model #
  model.name='Sigmoid\n(logisitic)'
  f=as.formula(paste0(y,"~","SSlogis(",x,",Asym,xmid,scal)"))
  init.params = getInitial( f,data=xydata$df)
  m = nls(formula = f, start = init.params, data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model =model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

make_poly_fit = function(input,    # Input data
                         y=Y, x=X, # X/Y Variables to fit
                         deg=3,
                         only.params=T){
  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']

  # model #
  model.name=sprintf('Polynomial\n(d = %s)',deg)
  f=as.formula(paste0(y,"~poly(",x,",degree=",deg,",raw=T)"))
  m = glm(formula =f , data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model = model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

make_expo_fit = function(input,    # Input data
                         y=Y, x=X, # X/Y Variables to fit
                         only.params=T){
  xydata = get_XY_data(input,x,y,noNA=T)
  mu.y = xydata$mu['y']
  var.y=xydata$var['y']
  nXY = xydata$n['xy']

  fexpo = function(alpha=1,beta=0,x){ return(alpha * exp(beta*x)) }
  # model #
  model.name='Exponential'
  f=as.formula(paste0(y,"~ expo(alpha,beta,",x,")"))
  init.params=list(alpha=1,beta=-1)
  m=nls(f,start = init.params,data=xydata$df)
  m.params = get_model_params(m,xydata$XX,xydata$YY) %>%
    purrr::list_modify(model = model.name, xname=xydata$varnames['x'],yname=xydata$varnames['y'])
  if(only.params){ return(m.params) }
  # data with model #
  fit = get_fit_data(xydata,m) %>% mutate(model=model.name)
  return(fit)
}

decompose_variance = function(LM,to.df=F){
  # DECOMPOSE VARIANCE FROM LINEAR REGRESSION
  N=sum(complete.cases(LM$model))
  TSS = var( LM$model[,1] ) * (N-1)
  df.var = summary(aov(LM))[[1]]
  nvar = nrow(df.var)
  RSS = sum(df.var$`Sum Sq`[nvar])
  rss.pc =100*RSS/TSS
  ESS = sum(df.var$`Sum Sq`[-nvar])

  ess.max = sum( (fitted(LM) - mean_(LM$model[,1]))^2 )
  ess = TSS-RSS
  ess.pc =100*ess/TSS
  #RSS =  deviance(LM)
  #rss = sum( residuals(LM)^2)
  #ESS =  TSS - RSS
  one_line_formula =  paste(deparse1(formula(LM))) %>%
    str_trunc(side = 'center', width = 80)
  nterms = n_distinct(labels(LM))

  cat(sprintf("%s\n",one_line_formula))
  cat(sprintf("(%s predictor variables)\n",nterms))
  cat(sprintf("TSS %.1f (n=%s)\n",TSS,N))
  cat(sprintf("--> ESS %.1f (%.0f%%)\n",ess, ess.pc))
  cat(sprintf("--> RSS %.1f (%.0f%%)\n", RSS,  rss.pc))
  if(to.df){
    res=tibble(N=N,nterms=nterms,TSS=TSS,ESS=ESS,RSS=RSS,
               RSS_rel = rss.pc, ESS_rel=ess.pc
    )
    return(res)
  }
}
