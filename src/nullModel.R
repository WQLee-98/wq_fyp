library(utils)
library(parallel)
library(foreach)
library(doParallel)

# Create function that calculates mean fdisp values of null models of each cell
# and their z-scores
null_mod = function(orig_mat, defaun_mat, orig_fd, defaun_fd){
  # arguments
  #   orig_mat: site x species matrix of original assemblage
  #   defaun_mat: site x species matrix of defaunated assemblage
  #   orig_fd: results of functional dispersion analysis on original assemblage
  #   defaun_fd: results of functional dispersion analysis on defaunated assemblage
  
  num_sp_def = orig_fd$uniq.sp - defaun_fd$uniq.sp
  
  
  mean_fdisp = rep(NA, nrow(orig_mat)) ; names(mean_fdisp) = row.names(orig_mat)
  stan_dev = rep(NA, nrow(orig_mat)) ; names(stan_dev) = row.names(orig_mat)
  z_scores = rep(NA, nrow(orig_mat)) ; names(z_scores) = row.names(orig_mat)
  orig_fdisp = orig_fd$FDis
  
  vectors = orig_fd$vectors
  pos = orig_fd$eig > 0

  
  for (i in 1:length(num_sp_def)){
    if (num_sp_def[i] != 0){
      # identities of species in a cell
      sp = colnames(orig_mat)[as.logical(orig_mat[i,])]

      # combination of species defaunated
      combin = combn(sp, num_sp_def[i])
      
      # creating vector to store all functional dispersion values for all null models
      fdisp = rep(NA,ncol(combin))
      
      for (j in 1:ncol(combin)){
        pres = sp
        pres = pres[!pres %in% combin[,j]]
        
        vec = vectors[pres,,drop = F]
        centroid = apply(vec, 2, mean)
        
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
        fdisp[j] <- mean(zij)
      }
      mean_fdisp[i] = mean(fdisp)
      stan_dev[i] = sd(fdisp)
    }
    stopCluster(cl)
  }
  z_scores = (orig_fdisp - mean_fdisp)/stan_dev
  
  return(list(z_scores = z_scores, mean_fdisp = mean_fdisp))
}

# num_cores = detectCores()
# cl = makeCluster(num_cores-1, type = "PSOCK")
# print(cl)
# registerDoParallel(cl)
# getDoParRegistered()
# getDoParWorkers()
# 
# 
# 
# system.time(
#   x <- foreach(
#     i = 1:10, 
#     .combine = 'c'
#   ) %dopar% {
#     sqrt(i)
#   })
# 
# stopCluster(cl)
# 
# 
# system.time({
#   x = vector()
#   for(i in 1:100){
#     x[i] = sqrt(i)
#   }
# })

# clusterExport(cl, list("site_species_orig","site_species_HL","res_orig",
#                        "res_HL", "null_mod"))