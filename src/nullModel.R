
# Create function that calculates mean fdisp values of null models of each cell
null_mod = function(orig_mat, nb.sp, vectors, eig){
  # arguments
  #   orig_mat: site x species matrix of original assemblage
  #   nb.sp: vector storing the number of species present in each cell after a 
  #     defaunation scenario
  #   vectors = coordinates of all species on all 4 axis of PCoA
  #   eig = eigen values of PCoA

  pos = eig > 0
  
  # create list to store randomised fdisp values for each cell
  nulls = list()
  
  for (i in 1:nrow(orig_mat)){
    if (nb.sp[i] >= 2){
      # identities of species in a cell
      sp = colnames(orig_mat)[as.logical(orig_mat[i,])]
      
      # creating list of randomised combinations of species in the cell
      combin = lapply(1:999, FUN = function(x){sample(sp,nb.sp[i])})
      
      avg.dist = lapply(combin, FUN = function(x, vectors, pos){
        #calculate centroidd of the species combination
        pres = x
        vec = vectors[pres,,drop = F]
        centroid = apply(vec, 2, mean)
        
        # calculate fdisp of cell
        dist.pos <- sweep(vec[, pos , drop = F], 2, centroid[pos]) 
        dist.pos <- rowSums(dist.pos^2)
        if (any(!pos) ){
          dist.neg <- sweep(vec[, !pos , drop = F], 2, centroid[!pos])
          dist.neg <- rowSums(dist.neg^2)
        }
        else dist.neg <- 0
        
        # euclidian distance of all species in community i to centroid through Pythagoras' Theorem
        zij <- sqrt(abs(dist.pos - dist.neg))
        return (mean(zij))
      }, vectors = vectors, pos = pos)
      
      avg.dist = unlist(avg.dist)
      nulls[[i]] = avg.dist
    } 
    else if (nb.sp[i] == 1){
      # functional dispersion of a cell with 1 species is 0
      avg.dist = rep(0, 999)
      nulls[[i]] = avg.dist
    } 
    else {
      # functional dispersion of a cell with 0 species is NA
      avg.dist = rep(NA, 999)
      nulls[[i]] = avg.dist
    }
  }
  # replace names of nulls with grid ids
  names(nulls) = rownames(orig_mat)
  
  # returns list of all null values of each cell
  return (nulls)
}
