library(ade4)


# creating FD function that returns functional dispersion of all communities
funcdiv = function (vectors, eig, a){
  # Arguments
  #   vectors = coordinates of all species on all 4 axis of PCoA
  #   eig = eigen values of PCoA
  #   a = site x species matrix
  
  pos = eig > 0
  
  # create vector to store functional dispersion values
  avg.dist.cent = rep(NA, nrow(a)) ; names(avg.dist.cent) = row.names(a)
  
  # create list to store euclidean distances of each species to each community centroid
  dist.cent = list()
  
  # number of communities/cells
  com = nrow(a)
  
  for (i in 1:com){
    # index of species present in community i
    pres = which(a[i ,] > 0)
    # coordinates of all species present in community i
    vec = vectors[pres, , drop = F]
    
    # number of unique species in community i (2 diff species may have the exact same traits)
    nb.sp = nrow(unique(vec))
    
    # if community has at least 2 unique species, calculate distances to centroid
    # and functional dispersion
    if (nb.sp >= 2){
      # community i centroid coordinates
      centroid = apply(vec, 2, mean)
      
      # distance of species from community i centroid
      dist.pos = sweep(vec[, pos , drop = F], 2, centroid[pos]) 
      dist.pos = rowSums(dist.pos^2)
      if (any(!pos) ){
        dist.neg = sweep(vec[, !pos , drop = F], 2, centroid[!pos])
        dist.neg = rowSums(dist.neg^2)
      }
      else dist.neg = 0
      
      # euclidian distance of all species in community i to centroid through Pythagoras' Theorem
      zij = sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] = mean(zij)
      dist.cent[[i]] = zij
    } 
    else if (nb.sp == 1){
      # functional dispersion of a cell with 1 species is 0
      avg.dist.cent[i] = 0
      sp = c(0)
      names(sp) = rownames(vec)
      dist.cent[[i]] = sp
    }
    else {
      # functional dispersion of a cell with 0 species is NA
      avg.dist.cent[i] = NA
      dist.cent[[i]] = NA
    }
  }
  
  # changing element names in list 'dist.cent'
  names(dist.cent) = rownames(a)
  
  # returns
  #   FDis: Average distance to centre of centroid (Functional Dispersion) of each cell
  #   dist.cent: Distance of all species to centroid of each cell
  return(list(FDis = avg.dist.cent, dist.cent = dist.cent))
}
