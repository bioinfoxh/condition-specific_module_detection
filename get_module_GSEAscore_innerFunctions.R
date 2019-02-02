


# get_module_GSEAscore <- function(moduleList, exp, label, time_series, degContrast){
# 
#   design <- model.matrix(~0+factor(label, level=time_series))
#   colnames(design) <- levels( factor(label,level=time_series) )
#   fit <- lmFit(exp, design)
#   # To make all pair-wise comparisons between the three groups the appropriate contrast matrix can be created by
#   contrast.matrix <- makeContrasts(contrasts = degContrast, levels=design)
#   fit2 <- contrasts.fit(fit, contrast.matrix)
#   fit2 <- eBayes(fit2)
#   # A list of top genes differential expressed in group2 versus group1 can be obtained from
#   deg_group1 = topTable(fit2, coef=1, adjust="BH",n=Inf)
#   deg_group2 = topTable(fit2, coef=2, adjust="BH",n=Inf)
#   deg_group3 = topTable(fit2, coef=3, adjust="BH",n=Inf)
#   
#   deg_group1 = deg_group1[ order(deg_group1$t, decreasing = TRUE)  ,]
#   deg_group2 = deg_group2[ order(deg_group2$t, decreasing = TRUE)  ,]
#   deg_group3 = deg_group3[ order(deg_group3$t, decreasing = TRUE)  ,]
#   
#   gsea_obs_group1= sapply(moduleList, function(x) GSEA.EnrichmentScore( rownames(deg_group1) , x, 1, abs(deg_group1$logFC) )$ES )
#   gsea_obs_group2= sapply(moduleList, function(x) GSEA.EnrichmentScore( rownames(deg_group2) , x, 1, abs(deg_group2$logFC) )$ES )
#   gsea_obs_group3= sapply(moduleList, function(x) GSEA.EnrichmentScore( rownames(deg_group3) , x, 1, abs(deg_group3$logFC) )$ES )
#   
#   gsea_obs = data.frame( gsea_obs_group1,gsea_obs_group2,gsea_obs_group3, stringsAsFactors = FALSE )
#   colnames(gsea_obs) = degContrast
#   
#   return(gsea_obs)
# 
# }



GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}





empirical_p_value <- function( score_obs, score_random){
  if(score_obs>=0){
    p = signif(sum(score_random > score_obs)/sum(score_random>=0)   , digits=5)
  }else{
    p = signif(sum(score_random < score_obs)/sum(score_random<0)   , digits=5)
  }
  return(p)
} 


normlize_GSEAscore <- function( score_obs, score_random){
  
  mean_pos = abs( mean( score_random[ score_random>=0  ]  ) )
  mean_neg = abs( mean( score_random[ score_random<0  ]  ) )
  
  score_random_norm_pos = score_random[score_random>=0] / mean_pos
  score_random_norm_neg = score_random[score_random<0] / mean_neg

  if(score_obs>=0){
    score_obs_norm = score_obs/mean_pos
  }else{
    score_obs_norm = score_obs/mean_neg
  }

  return( list(NormScore_obs = score_obs_norm, 
               NormScore_random_pos = score_random_norm_pos, 
               NormScore_random_neg = score_random_norm_neg ) )
}




# normlize_GSEAscore_Z <- function( score_obs, score_random){
#   
#   mean_pos =  mean( score_random[ score_random>=0  ]  ) 
#   mean_neg =  mean( score_random[ score_random<0  ]  ) 
#   
#   sd_pos = sd( score_random[ score_random>=0  ]  )
#   sd_neg = sd( score_random[ score_random<0  ]  )
#   
#   score_random_norm_pos = (score_random[score_random>=0] - mean_pos)/sd_pos
#   score_random_norm_neg = (score_random[score_random<0] - mean_neg)/sd_neg
#   
#   if(score_obs>=0){
#     score_obs_norm = (score_obs - mean_pos)/sd_pos
#   }else{
#     score_obs_norm = (score_obs - mean_neg)/sd_neg
#   }
#   return( list(NormScore_obs = score_obs_norm, 
#                NormScore_random_pos = score_random_norm_pos, 
#                NormScore_random_neg = score_random_norm_neg ) )
# }



get_GSEA_fdr <- function( score_obs_vec, score_random_matrix ){
  
  normGSEA_obs = c()
  normGSEA_random_pos = c()
  normGSEA_random_neg = c()
  
  for (j in 1:length(score_obs_vec) ){
    temp_normGSEA = normlize_GSEAscore( score_obs_vec[j], score_random_matrix[j,]  )
    
    normGSEA_obs = c( normGSEA_obs,  temp_normGSEA$NormScore_obs)
    normGSEA_random_pos = c( normGSEA_random_pos, temp_normGSEA$NormScore_random_pos )
    normGSEA_random_neg = c( normGSEA_random_neg, temp_normGSEA$NormScore_random_neg )
  }
  
  normGSEA_obs_pos = normGSEA_obs[ score_obs_vec>=0  ]
  normGSEA_obs_neg = normGSEA_obs[ score_obs_vec<0  ]
  
  GSEA_fdr = c()
  for (k in 1:length(score_obs_vec) ){
    if( score_obs_vec[k]>=0  ){
      temp_ratio_obs = sum( normGSEA_obs_pos >= normGSEA_obs[k] )/length(normGSEA_obs_pos)
      temp_ratio_perm = sum( normGSEA_random_pos >= normGSEA_obs[k] ) / length(normGSEA_random_pos)
      temp_fdr = temp_ratio_perm/temp_ratio_obs
    }else{
      temp_ratio_obs = sum( normGSEA_obs_neg <= normGSEA_obs[k] )/length(normGSEA_obs_neg)
      temp_ratio_perm = sum( normGSEA_random_neg <= normGSEA_obs[k] ) / length(normGSEA_random_neg)
      temp_fdr = temp_ratio_perm/temp_ratio_obs   
    }
    GSEA_fdr = c(GSEA_fdr, temp_fdr)
  }
  return(GSEA_fdr)
}



