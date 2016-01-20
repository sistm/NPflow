#'Compute a limited F-measure
#'
#'@name Flimited
#'
#'@param n_small_clst
#'
#'@param partition_ref
#'
#'@param partition_est
#'
#'
#'
#'@export Flimited
#'
#'
#'

Flimited<-function(n_small_clst,partition_ref,partition_est){
  
  partition_ref<-as.numeric(partition_ref)
  partition_est<-as.numeric(partition_est)
  
  label_obs_in_small_class_ref<-which(table(partition_ref)<=n_small_clst)
  if(!length(label_obs_in_small_class_ref)){
    stop('Error : n_small_clst should be increased')
  }
  
  label_obs_in_small_class_ref <- as.numeric(names(label_obs_in_small_class_ref))
  index_obs_in_small_clust_ref <- which(partition_ref %in% label_obs_in_small_class_ref)
  label_obs_in_small_class_est <- unique(partition_est[index_obs_in_small_clust_ref]) 
  
  partition_est_restricted<-partition_est[partition_est%in%label_obs_in_small_class_est]
  
  
  find_index_obs_coclust<-partition_est%in%label_obs_in_small_class_est
  find_index_obs_coclust[index_obs_in_small_clust_ref]<-FALSE
  index_obs_in_small_clust_ref_2<- which(find_index_obs_coclust==TRUE)
  
  partition_ref_restricted<-c(partition_ref[sort(c(index_obs_in_small_clust_ref,
                                                   index_obs_in_small_clust_ref_2))])
  
  return(FmeasureC(ref=partition_ref_restricted, pred=partition_est_restricted))
}
