

source("get_GSEA_enrichment_score.R")

dir.create("pivot_gseaEnrichedModule/")
setwd("pivot_gseaEnrichedModule/")


gene2mod = c()

for (i in 1:length(moduleList)){
  temp_gene = moduleList[[i]]
  temp_size = length(temp_gene)
  temp_index = rep( i,temp_size )
  temp_gene2mod = data.frame(GeneID=temp_gene,moduleID=temp_index,stringsAsFactors = FALSE)
  gene2mod = rbind(gene2mod,temp_gene2mod)
}
gene2mod = gene2mod[  is.element( gene2mod[,2],pathwayID )  , ]
gene2mod = gene2mod[  order(gene2mod[,1]) , ]
gene2mod = data.frame( gene = gene2mod[,1], mod=gene2mod[,2], stringsAsFactors = F )

gene_all = genelist
gene_intra = unique(gene2mod[,1])
gene_inter = setdiff( gene_all,gene_intra )


temp = ppi[, c(2,1,3)]
colnames(temp) = colnames(ppi)[1:3]
network_double = rbind( ppi[,1:3] , temp   )
network_double[,3] = abs(network_double[,3])
temp = table( network_double[,1]  )
degree = data.frame(gene=names(temp) ,temp,stringsAsFactors = F)[,c(1,3)]


hub_cutoff = quantile( degree[,2], 0.9  )


hub = degree[ degree[,2]>=hub_cutoff  ,]
hub_inter = hub[ hub[,1] %in% gene_inter   ,]
hub_intra = hub[ hub[,1] %in% gene_intra   ,]


dir.create("pivot_hub")
setwd("pivot_hub")


temp_deg = deg_list$`X8cell-X4cell`

hub2mod_p = c()
hub2mod_ES = c()
hub_pivot = c()
for( i in 1:nrow(hub_inter)  ){
 # print(i)
  
  temp_hub = hub_inter[i,1]
  temp_hub_degree = hub_inter[i,2]
  
  temp_hub_ppi = network_double[ network_double[,1] %in% temp_hub  ,]
  temp_hub_deg = temp_deg[  match(temp_hub_ppi[,2],rownames(temp_deg) ) ,]
  temp_hub_ppi$t = temp_hub_deg$t
  temp_hub_ppi = temp_hub_ppi[ order( abs(temp_hub_ppi[,4]), decreasing = T )  , ]
  
  temp_hub_nb2mod = gene2mod[ gene2mod[,1] %in% temp_hub_ppi[,2], ]
  if( nrow(temp_hub_nb2mod) == 0 ){next}
  
  temp_hub_mod = as.data.frame( table( temp_hub_nb2mod[,2] ), stringsAsFactors = FALSE )
  if( sum( temp_hub_mod[,2]>=3 ) ==0  ){
    next
  }else{  
      temp_hub_mod_ge3 = as.numeric(  temp_hub_mod[ temp_hub_mod[,2] >=3 ,1 ]  )  
      
  } 
  
  temp_ESscore = c()
  for( j in 1:length(temp_hub_mod_ge3)  ){
    temp_mod_gene = temp_hub_nb2mod[temp_hub_nb2mod[,2]==temp_hub_mod_ge3[j]   ,1]
#    temp_mod_score = GSEA.EnrichmentScore( temp_hub_ppi[,2], temp_mod_gene, 1, temp_hub_ppi[,3] )
    temp_mod_score = GSEA.EnrichmentScore( temp_hub_ppi[,2], temp_mod_gene, 1, temp_hub_ppi$corR )
    temp_ESscore = c(temp_ESscore, temp_mod_score$ES)
  }
  ESscore_obs = data.frame( mod = temp_hub_mod_ge3, ES= temp_ESscore, stringsAsFactors = F  )
  
  
########  random  
  
  # exp_hub = exp_gene[ exp_gene[,1]==temp_hub ,   ]
  # exp_nb = exp_gene[ match( temp_hub_ppi[,2] , exp_gene[,1] ),    ]  
  # 
  # rand_exp_nb = exp_nb[,  2:ncol(exp_nb)]
  # rownames(rand_exp_nb) = exp_nb[,1]
  
  ESscore_rand = c()
  
  for(k in 1:1000){
    if(k%%100 == 0)print(paste( i, "gene, perm", k, sep = " "   )  )
    
#    rand_exp_hub = sample( exp_hub[2:length(exp_hub)]  )
    
#    rand_exp_nb = apply( exp_nb[,  2:ncol(exp_nb)] , 1, function(x) sample(x)  )
  
#    rand_cor = abs( cor( as.numeric(rand_exp_hub) , t(as.matrix(rand_exp_nb) ) )[1,] )
    
  
    rand_hub_ppi = cbind( temp_hub_ppi[,1:3], t=sample(temp_hub_ppi$t))
    rand_hub_ppi = rand_hub_ppi[ order(abs(rand_hub_ppi$t),decreasing = TRUE)  ,]
    
    temp_ESscore = c()
    for( m in 1:length(temp_hub_mod_ge3)  ){
      temp_mod_gene = temp_hub_nb2mod[temp_hub_nb2mod[,2]==temp_hub_mod_ge3[m]   ,1]
#      temp_mod_score = GSEA.EnrichmentScore( rand_hub_ppi[,2], temp_mod_gene, 1, rand_hub_ppi[,3] )
      temp_mod_score = GSEA.EnrichmentScore( rand_hub_ppi[,2], temp_mod_gene, 1, rand_hub_ppi$corR )
      temp_ESscore = c(temp_ESscore, temp_mod_score$ES)
    }
    
    ESscore_rand = cbind(ESscore_rand, as.matrix(temp_ESscore) )
    rownames(ESscore_rand) = temp_hub_mod_ge3
    
   }
  ESscore_temp = cbind( ESscore_obs[,2], ESscore_rand )
  p = apply(ESscore_temp, 1, function(x) empirical_p_value( x[1],x[2:length(x)]   )  )
  p = cbind(ESscore_obs,p  )
  
  temp_p = p[  match( enriched_module_id, p[,1] )  ,  3]
  temp_ES = p[  match( enriched_module_id, p[,1] )  ,  2]
  
  hub2mod_p = rbind( hub2mod_p, temp_p  )
  hub2mod_ES = rbind( hub2mod_ES, temp_ES  )
  hub_pivot = c( hub_pivot, temp_hub  )
  
  
  rownames(hub2mod_ES) = hub_pivot
  rownames(hub2mod_p) = hub_pivot
  colnames(hub2mod_ES) = enriched_module_id
  colnames(hub2mod_p) = enriched_module_id
  
  temp_output = data.frame( module=rownames(ESscore_temp), ESscore_temp, stringsAsFactors = F)
  colnames(temp_output) = c("module","obs",paste("random_", 1:(ncol(temp_output)-2 ), sep="") )
  
  temp_output_file = paste("ES_score_hub_",temp_hub,".txt",sep="")
  write.table(temp_output, temp_output_file, quote=F, sep="\t", col.names = T, row.names = F)
  }





empirical_p_value <- function( score_obs, score_random){
  if(score_obs>=0){
    p = signif(sum(score_random >= score_obs)/sum(score_random>=0)   , digits=5)
  }else{
    p = signif(sum(score_random <= score_obs)/sum(score_random<0)   , digits=5)
  }
  return(p)
} 



#### output the hub2gene result

temp = data.frame( gene = rownames(hub2mod_ES), hub2mod_ES, stringsAsFactors = F )
write.table( temp, "hub2mod_ES.txt", sep="\t", quote=F, row.names=F, col.names=T  )
#
temp = data.frame( gene = rownames(hub2mod_p), hub2mod_p, stringsAsFactors = F )
write.table( temp, "hub2mod_p.txt", sep="\t", quote=F, row.names=F, col.names=T  )


######## get pivot gene with ES score greater than 0 (overrepresented)


temp1 = hub2mod_ES
temp2 = hub2mod_ES>0
temp2[temp2==FALSE] = NA

hub2mod_ES_pos = hub2mod_ES * temp2
hub2mod_p_pos = hub2mod_p * temp2

temp = data.frame( gene = rownames(hub2mod_ES_pos), hub2mod_ES_pos, stringsAsFactors = F )
write.table( temp, "hub2mod_ES_pos.txt", sep="\t", quote=F, row.names=F, col.names=T  )
temp = data.frame( gene = rownames(hub2mod_p_pos), hub2mod_p_pos, stringsAsFactors = F )
write.table( temp, "hub2mod_p_pos.txt", sep="\t", quote=F, row.names=F, col.names=T  )



temp = hub2mod_p_pos
temp[is.na(temp)]=2
pivot = rownames(temp[(rowSums(temp<=0.05)>0),])
pivot2 = rownames(temp[(rowSums(temp<=0.05)>1),])


temp1 = which(temp<=0.05, arr.ind = T, useNames = T)
temp2 = rownames(temp)[temp1[,1]]
temp3 = colnames(temp)[temp1[,2]]
temp4 = gene2symbol[ match( temp2, gene2symbol[,1] ),  ]

hub2mod_pivot = data.frame( temp4,module=temp3, stringsAsFactors = F )
hub2mod_pivot2 = hub2mod_pivot[ hub2mod_pivot[,1]%in%pivot2  ,]
hub2mod_pivot2 = hub2mod_pivot2[ order(hub2mod_pivot2[,1])  ,]

setwd("..")
write.table(hub2mod_pivot, "hub2mod_pivot.txt",sep="\t", quote=F, row.names=F, col.names=T)
write.table(hub2mod_pivot2, "hub2mod_pivot2.txt",sep="\t", quote=F, row.names=F, col.names=T)






