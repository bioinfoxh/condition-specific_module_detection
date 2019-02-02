

source("get_deg_list.R")
source("get_module_GSEAscore_innerFunctions.R")
#source("get_GSEA_enrichment_score.R")

library(limma, quietly = T)





get_module_GSEAscore <- function(moduleList, exp, label, level, contrast, permTimes){
  
  if( is.null(moduleList) ){ names(moduleList)=1:length(moduleList)  }


deg_list = get_deg_list( exp,label, level, contrast  )

module_GSEAscore = c()
for (i in 1:length(deg_list)){
  temp_deg = deg_list[[i]]
  temp_deg = temp_deg[ order( temp_deg$t, decreasing=TRUE ) , ]
  temp_gene_order = rownames(temp_deg)
  temp_corVec = abs(temp_deg$t)
  temp_GSEAscore = sapply(moduleList, function(x) GSEA.EnrichmentScore( temp_gene_order, x, 1, correl.vector = temp_corVec )$ES )
  
  module_GSEAscore = cbind(module_GSEAscore, temp_GSEAscore)
}
colnames(module_GSEAscore) = names( deg_list)
rownames(module_GSEAscore) = names(moduleList)



perm_module_GSEAscore_list = list()
for (j in 1:permTimes){
  print(paste("GSEA permutation times ",j,sep="" ) )
  perm_label = sample(label)
  perm_deg_list = get_deg_list( exp, perm_label, level, contrast)
  
  perm_module_GSEAscore = c()
  for (i in 1:length(perm_deg_list)){
    temp_deg = perm_deg_list[[i]]
    temp_deg = temp_deg[ order( temp_deg$t, decreasing=TRUE ) , ]
    temp_gene_order = rownames(temp_deg)
    temp_corVec = abs(temp_deg$t)
    temp_GSEAscore = sapply(moduleList, function(x) GSEA.EnrichmentScore( temp_gene_order, x, 1, correl.vector = temp_corVec)$ES )
    
    perm_module_GSEAscore = cbind(perm_module_GSEAscore, temp_GSEAscore)
  }
  colnames(perm_module_GSEAscore) =names( perm_deg_list)
  perm_module_GSEAscore_list[[j]] = perm_module_GSEAscore
}



module_GSEAp=c()
module_GSEAfdr=c()
for( i in 1:ncol(module_GSEAscore)){
  temp_gsea_obs = module_GSEAscore[,i]
  temp_gsea_perm = sapply( perm_module_GSEAscore_list, function(x) x[,i])
  
  temp = cbind( temp_gsea_obs, temp_gsea_perm )
  temp_module_gseaP = apply(temp, 1, function(x) empirical_p_value(x[1],x[2:ncol(temp)]  ) )
  module_GSEAp = cbind(module_GSEAp, temp_module_gseaP)
  
  temp_GSEAfdr = get_GSEA_fdr(temp_gsea_obs, temp_gsea_perm)
  module_GSEAfdr = cbind(module_GSEAfdr, temp_GSEAfdr)
}
colnames(module_GSEAp) = colnames(module_GSEAscore)
rownames(module_GSEAp) = names(moduleList)

colnames(module_GSEAfdr) = colnames(module_GSEAscore)
rownames(module_GSEAfdr) = names(moduleList)

GSEAresult = list( DEG_list = deg_list,
               GSEA_score = module_GSEAscore,
               GSEA_p = module_GSEAp,
               GSEA_fdr = module_GSEAfdr )

return(GSEAresult)

}


# 
# write.table(module_GSEAp,"module_GSEAp.txt", sep="\t", quote=F, row.names=F, col.names=T)
# write.table(module_GSEAfdr,"module_GSEAfdr.txt", sep="\t", quote=F, row.names=F, col.names=T)
# write.table(module_GSEAscore,"module_GSEAscore.txt", sep="\t", quote=F, row.names=F, col.names=T)
# for( i in 1:ncol(module_GSEAscore)){
#   temp_gsea_perm = sapply( perm_module_GSEAscore_list, function(x) x[,i])
#   temp_filename = paste("perm_module_GSEAscore_contrast_",i,".txt",sep="")
#   write.table(temp_gsea_perm, temp_filename, sep="\t",quote=F, row.names=F, col.names=F)
# }


