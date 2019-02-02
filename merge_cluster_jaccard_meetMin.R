

merge_cluster <- function(cluster_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=TRUE){
  
  
  intraGene = unique(unlist(cluster_genelist))
  
  cluster_gene2clust = sapply( cluster_genelist, function(x) intraGene%in%x )
  rownames(cluster_gene2clust) = intraGene
  
  cluster_score = get_cluster_ovlp_score(cluster_gene2clust,clusterSizeDecreasing)
  
  cluster_gene2clust = cluster_score$cluster_gene2clust
  cluster_size = cluster_score$cluster_size
  cluster_jaccard = cluster_score$cluster_jaccard
  cluster_meetMin = cluster_score$cluster_meetMin
  
  # cluster_gene2clust = cluster_gene2clust[ , order( colSums(cluster_gene2clust) ,decreasing = T)  ]
  # cluster_size = colSums(cluster_gene2clust)
  # cluster_ovlp = t(cluster_gene2clust) %*% cluster_gene2clust
  # 
  # temp1 = !cluster_gene2clust
  # temp2 = t(temp1) %*% temp1
  # cluster_union = length(intraGenes)-temp2
  # cluster_jaccard = cluster_ovlp/cluster_union
  # diag(cluster_jaccard) = 0
  # 
  # temp1 = apply( as.matrix(cluster_size), 1, function(x) rep(x, nrow(cluster_size)) )
  # cluster_minSize=apply( cbind( t(temp1),temp1 ), 1, function(x) pmin( x[1:(length(x)/2) ] ,x[ (length(x)/2+1):length(x) ] )  )
  # cluster_meetMin = cluster_ovlp / cluster_minSize
  # diag(cluster_meetMin) = 0
  
  merge_log=c()
  while( (max(cluster_jaccard)>=jaccard_cutoff) | (max(cluster_meetMin)>=meetMin_cutoff) ){
    
    if(max(cluster_jaccard)>=jaccard_cutoff){
      temp_index = which(cluster_jaccard==max(cluster_jaccard), arr.ind = T)[1,]
      cluster_gene2clust[ ,temp_index[2] ]= (cluster_gene2clust[ ,temp_index[2] ] | cluster_gene2clust[ ,temp_index[1] ])        
      cluster_gene2clust=cluster_gene2clust[,-temp_index[1]]
      
      merge_log=rbind( merge_log,temp_index )
      
      cluster_score = get_cluster_ovlp_score(cluster_gene2clust,clusterSizeDecreasing)
      
      cluster_gene2clust = cluster_score$cluster_gene2clust
      cluster_size = cluster_score$cluster_size
      cluster_jaccard = cluster_score$cluster_jaccard
      cluster_meetMin = cluster_score$cluster_meetMin
      
      # cluster_ovlp = t(cluster_gene2clust) %*% cluster_gene2clust
      # 
      # temp1 = !cluster_gene2clust
      # temp2 = t(temp1) %*% temp1
      # cluster_union = length(intraGenes)-temp2
      # cluster_jaccard = cluster_ovlp/cluster_union
      # diag(cluster_jaccard) = 0
      
      next
    }
    
    
    if(max(cluster_meetMin)>=meetMin_cutoff){
      temp_index = which(cluster_meetMin==max(cluster_meetMin), arr.ind = T)[1,]
      cluster_gene2clust[ ,temp_index[2] ]= (cluster_gene2clust[ ,temp_index[2] ] | cluster_gene2clust[ ,temp_index[1] ])        
      cluster_gene2clust=cluster_gene2clust[,-temp_index[1]]
      
      merge_log=rbind( merge_log,temp_index )
      
      cluster_score = get_cluster_ovlp_score(cluster_gene2clust,clusterSizeDecreasing)
      
      cluster_gene2clust = cluster_score$cluster_gene2clust
      cluster_size = cluster_score$cluster_size
      cluster_jaccard = cluster_score$cluster_jaccard
      cluster_meetMin = cluster_score$cluster_meetMin
      
      next
    }
    
    
  }
  
  result = list(intraGenes = intraGene,
                mergedCluster_gene2clust = cluster_gene2clust,
                mergedCluster_jaccard = cluster_jaccard,
                mergedCluster_meetmin = cluster_meetMin,
                merge_log = merge_log)
  
}




get_cluster_ovlp_score <- function(cluster_gene2clust, clusterSizeDecreasing){

  
  cluster_gene2clust = cluster_gene2clust[ , order( colSums(cluster_gene2clust) ,decreasing = clusterSizeDecreasing)  ]
  
  cluster_size = colSums(cluster_gene2clust)
  
  cluster_ovlp = t(cluster_gene2clust) %*% cluster_gene2clust
  
  temp1 = !cluster_gene2clust
  temp2 = t(temp1) %*% temp1
  cluster_union = nrow(cluster_gene2clust) - temp2
  cluster_jaccard = cluster_ovlp/cluster_union
  diag(cluster_jaccard) = 0
  
  temp1 = apply( as.matrix(cluster_size), 1, function(x) rep(x, length(cluster_size)) )
  cluster_minSize=apply( cbind( t(temp1),temp1 ), 1, function(x) pmin( x[1:(length(x)/2) ] ,x[ (length(x)/2+1):length(x) ] )  )
  cluster_meetMin = cluster_ovlp / cluster_minSize
  diag(cluster_meetMin) = 0
  
  result = list( cluster_gene2clust = cluster_gene2clust,
                 cluster_size = cluster_size,
                 cluster_jaccard = cluster_jaccard,
                 cluster_meetMin = cluster_meetMin)
  
}

