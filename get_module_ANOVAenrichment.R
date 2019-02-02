
source("get_module_GSEAscore_innerFunctions.R")


get_geneANOVA <- function(x, exp_label){
  x=as.numeric(x)
  temp = data.frame(x=x, stage = exp_label,stringsAsFactors = F)
  temp1 = summary( aov(x ~ stage, temp) )
  stat = c(temp1[[1]][1,4], temp1[[1]][1,5])
  return(stat)
}

#gene_anova = t( apply( exp_gene, 1, function(x) get_geneANOVA(x,exp_label) ) )

get_module_ANOVAenrichment <- function(moduleList, exp, label, permTimes){
  
  if( is.null(moduleList) ){ names(moduleList)=1:length(moduleList)  }
  
  gene_anova = t( apply( exp, 1, function(x) get_geneANOVA(x,label) ) )
  rownames(gene_anova) = rownames(exp)
  temp_stat = gene_anova[ order(gene_anova[,1],decreasing = TRUE) , ]
  gene_order = rownames(temp_stat)
  gene_weight = temp_stat[,1]
  module_GSEAscore = sapply(moduleList, function(x) GSEA.EnrichmentScore( gene_order, x, 1, correl.vector = gene_weight )$ES )
  

  perm_module_GSEAscore=c()
  for (j in 1:permTimes){
    print(paste("GSEA permutation times ",j,sep="" ) )
    perm_label = sample(label)
    
    perm_anova = t( apply( exp, 1, function(x) get_geneANOVA(x,perm_label) ) )
    rownames(perm_anova) = rownames(exp)
    temp_stat = perm_anova[ order(perm_anova[,1],decreasing = TRUE) ,]
    temp_gene_order = rownames(temp_stat)
    temp_gene_weight = temp_stat[,1]
    temp_GSEAscore = sapply(moduleList, function(x) GSEA.EnrichmentScore( temp_gene_order, x, 1, correl.vector = temp_gene_weight )$ES )
    
    perm_module_GSEAscore = cbind(perm_module_GSEAscore, temp_GSEAscore)
    

  }
  
  temp = cbind( module_GSEAscore, perm_module_GSEAscore )
  
  module_GSEAp = apply(temp, 1, function(x) empirical_p_value(x[1],x[2:ncol(temp)]  ) )
  
  module_GSEAfdr = get_GSEA_fdr(module_GSEAscore, perm_module_GSEAscore)
  
  
  GSEAresult = list( gene_anova = gene_anova,
                     GSEA_score = module_GSEAscore,
                     GSEA_p = module_GSEAp,
                     GSEA_fdr = module_GSEAfdr )
  
  return(GSEAresult)
  
}




