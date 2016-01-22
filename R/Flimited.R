#'Compute a limited F-measure
#'
#'A limited version of F-measure for
#'
#'@name Flimited
#'
#'@param n_small_clst the size of the 
#'
#'@param partition_ref
#'
#'@param partition_est
#'
#'
#'
#'@export
Flimited <- function(n_small_clst, pred, ref){
  
  partition_ref <- as.numeric(ref)
  partition_est <- as.numeric(pred)
  
  label_obs_in_small_class_ref <- which(table(partition_ref)<=n_small_clst)
  if(!length(label_obs_in_small_class_ref)){
    stop('n_small_clst should be increased')
  }
  
  label_obs_in_small_class_ref <- as.numeric(names(label_obs_in_small_class_ref))
  index_obs_in_small_class_ref <- which(partition_ref %in% label_obs_in_small_class_ref)
  label_obs_in_small_class_est <- unique(partition_est[index_obs_in_small_class_ref]) 
  
  index_obs_restricted <- which(partition_est%in%label_obs_in_small_class_est)
  partition_est_restricted <- partition_est[index_obs_restricted]
  partition_ref_restricted <- partition_ref[index_obs_restricted]
  
  return(FmeasureC(ref=partition_ref_restricted, pred=partition_est_restricted))
}
