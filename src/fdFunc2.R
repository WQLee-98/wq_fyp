library(ade4)


# creating FD function that returns functional dispersion of all communities
funcdiv = function (x, a, original = FALSE){
  # Arguments
  #   x = species x traits matrix/data frame. all traits need to be numeric and have no NA
  #   a = site x species matrix/data frame
  #   original = whether matrices inputted are original. Default is 'FALSE'
  
  # get tolerance value specifying small values that should be 0
  tol = 1e-07
  
  x = data.frame(x)
  a = as.matrix(a)
  
  ## checks on validity of arguments if using original species compositions=====
  if(original){
    # check presence of row names in 'x' i.e. species names
    if (is.null(row.names(x) ) ) stop("'x' must have row names.","\n")
    
    # check if one community has total abundance of zero (no species) *MAY NOT NEED
    abun.sum <- apply(a, 1, sum)
    if (any(abun.sum == 0) ) stop("At least one community has zero-sum 
                                  abundances (no species).","\n")
    
    # check if one species has total abundance of zero (never occurs) *MAY NOT NEED
    abun.sum2 <- apply(a, 2, sum)
    if (any(abun.sum2 == 0) ) stop("At least one species does not occur in any 
                                   community (zero total abundance across all 
                                   communities).","\n")
    
    # check if no. of species are same for 'x' and 'a' and in same order
    x.spec = rownames(x)
    a.spec = colnames(a)
    if (length(x.spec) != length(a.spec)) stop("Number of species for 'x' and 
                                               'a' do not match.","\n")
    if (any(x.spec != a.spec)) stop("Order of species in 'x' and 'a' do not 
                                    match","\n")
    
    # check whether any species has no traits
    no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]) )
    if (any(no.traits == 0) ) stop("At least one species has no trait data.","\n")
  }
  
  
  
  ## creating distance matrix 'x.dist'==========================================
  # distance measure = euclidean
  x.s <- apply(x, 2, scale, center = TRUE, scale = TRUE)
  x.dist <- dist(x.s)
  
  x.rn <- row.names(x)
  attr(x.dist, "Labels" ) <- x.rn
  
 
  
  ## Functional Dispersion =====================================================
  n <- attr(x.dist, "Size")
  
  # create a matrix of distances with transformations (not sure why)
  A <- matrix(0, ncol = n, nrow = n)
  
  A[row(A) > col(A)] <- -0.5 * x.dist^2 #why?
  
  A <- A + t(A) # t() transposes the dataframe
  
  # gower's double-centering
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE) # i assume this does the ordination
  
  vectors <- e$vectors
  eig <- e$values
  
  # check if 'x.dist' is Euclidean or not
  w0 <- eig[n] / eig[1]
  if (w0 > -tol) r <- sum(eig > (eig[1] * tol)) else r <- length(eig)
  
  # PCoA axes; coordinates of all species in morphospace
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), r)
  dimnames(vectors) <- list(colnames(a), NULL)
  
  pos <- eig > 0
  
  
  # create vector to store functional dispersion values
  avg.dist.cent <- rep(NA, nrow(a)) ; names(avg.dist.cent) <- row.names(a)
  
  # create list to store euclidean distances of each species to each community centroid
  dist.cent = list()
  
  # create vector to store number of unique species
  uniq.sp <- rep(NA, nrow(a)) ; names(uniq.sp) <- row.names(a)
  
  com <- nrow(a)
  
  for (i in 1:com){
    # index of species present in community i
    pres <- which(a[i ,] > 0)
    # coordinates of all species present in community i
    vec <- vectors[pres, , drop = F] 
    # number of unique species in community i (2 diff species may have the exact same traits)
    nb.sp <- nrow(unique(vec))
    uniq.sp[i] = nb.sp
    
    # if community has at least 2 unique species, calculate distances to centroid
    # and functional dispersion
    if (nb.sp >= 2){
      # community i centroid coordinates
      centroid <- apply(vec, 2, mean)
      
      # distance of species from community i centroid
      dist.pos <- sweep(vec[, pos , drop = F], 2, centroid[pos]) 
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos) ){
        dist.neg <- sweep(vec[, !pos , drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      
      # euclidian distance of all species in community i to centroid through Pythagoras' Theorem
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- mean(zij)
      dist.cent[[i]] = zij
    }
    else if (nb.sp == 1){
      avg.dist.cent[i] = 0
      sp = c(0)
      names(sp) = rownames(vec)
      dist.cent[[i]] = sp
    }
    else {
      avg.dist.cent[i] = NA
      dist.cent[[i]] = NA
    }
  }
  
  # changing element names in list 'dist.cent'
  names(dist.cent) = rownames(a)
  
  # returns
  #   FDis: Average distance to centre of centroid (Functional Dispersion) of each cell
  #   eig: Eigen values
  #   vectors: PCoA axes coordinates
  #   dist.cent: Distance of all species to centroid of each cell
  #   uniq.sp: Number of unique species in each cell
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors, dist.cent = dist.cent, uniq.sp = uniq.sp))
}


