find_hub_module_crosstalk = function ( network, pathway_genelist ,hub_cutoff , hub2mod_interaction_cutoff) {

# module_enriched=module_enriched[order(module_enriched)]
  
network_double = rbind( network[,1:2] , cbind(network[,2],network[,1] )   )
temp = table(network_double[,1] )
degree = as.data.frame(temp)

gene_all = unique( network_double[,1]  )
gene_intra = unique(unlist(pathway_genelist))
gene_inter = setdiff( gene_all,gene_intra )
#module_ge06_enriched = module_ge06[  module_enriched , ]

gene2mod=c()
for(i in 1:length(pathway_genelist)){
  temp_size = length( pathway_genelist[[i]] )
  temp_idx = rep( names(pathway_genelist)[i], temp_size )
  temp_gene2mod = data.frame( pathway_genelist[[i]], temp_idx, stringsAsFactors = F )
  gene2mod = rbind(gene2mod, temp_gene2mod)
}
colnames(gene2mod) = c("GeneID","PathwayName")

# gene2mod = sapply( pathway_genelist, function(x) gene_intra%in%x )
# rownames(gene2mod) = gene_intra
# colnames(gene2mod) = names(pathway_genelist)

# gene2mod = c()
# for (i in 1:nrow(module_ge06)){
# 	temp_gene = module_ge06[ i, !is.na(  module_ge06[i,] )]
# 	temp_size = length(temp_gene)
# 	temp_index = rep( i,temp_size )
# 	temp_gene2mod = cbind(temp_gene,temp_index)
# 	gene2mod = rbind(gene2mod,temp_gene2mod)
# }
# gene2mod = gene2mod[  is.element( gene2mod[,2],module_enriched )  , ]
# gene2mod = gene2mod[  order(gene2mod[,1]) , ]


# temp = degree[ match(gene_intra,degree[,1] ) , 2 ]
# gene2mod_degree = cbind( gene2mod, temp )

# degree_module = tapply( gene2mod_degree[,3], gene2mod_degree[,2], sum  )
# temp_mid = as.integer (names(degree_module) )
# edge_module = cbind( temp_mid,degree_module )
# edge_module = edge_module[ order( edge_module[,1] )  ,]

edge_module_out=c()
for (i in 1:length(pathway_genelist)){
  temp_module_gene = pathway_genelist[[i]]
  temp_module_ppi = network_double[ is.element( network_double[,1] , temp_module_gene  ) , ]
  temp_module_ppi_out = temp_module_ppi[  !is.element( temp_module_ppi[,2],temp_module_gene  )     ,]
  if( is.vector(temp_module_ppi_out)  ){
    temp_module_ppi_out_num=1
  }else{
    temp_module_ppi_out_num = nrow( temp_module_ppi_out )
  }
  edge_module_out = c( edge_module_out, temp_module_ppi_out_num     ) 
}
names(edge_module_out) = names(pathway_genelist)



edge_all = nrow( network )

hub10 = degree[ degree[,2]>=hub_cutoff  ,]
hub10_inter = intersect( gene_inter,hub10[,1]   ) 
hub10_intra = setdiff( hub10[,1],hub10_inter  )

p_hub2mod = c()

for (i in 1:length(hub10_inter) ){
	print(paste("processing hub ",i," / ",length(hub10_inter) ) )
	
	temp_hub = hub10_inter[i]
	temp_hub_ppi = network_double[ network_double[,1]==temp_hub , ]
	
	edge_hub = nrow( temp_hub_ppi )
	edge_nonhub = edge_all - edge_hub

	edge_hub = rep(edge_hub, length(pathway_genelist)  )
	edge_nonhub = rep( edge_nonhub, length(pathway_genelist) )

	edge_hub2mod = gene2mod[ is.element( gene2mod[,1],temp_hub_ppi[,2]  ), ]
	if( length(dim(edge_hub2mod))==0 ){
		edge_hub2mod = t( as.matrix(edge_hub2mod)  )
	}
	temp_edge_hub2mod = table( factor(edge_hub2mod[,2],level=names(pathway_genelist))   )
	#temp_mid = as.integer ( names(  temp_edge_hub2mod ) )
	edge_hub2mod = as.data.frame( temp_edge_hub2mod   )
	
	# edge_hub2mod = temp_edge_hub2mod_num
	# edge_hub2mod = edge_hub2mod[ order(edge_hub2mod[,1])  ,]
	
	num_matrix = cbind( edge_hub2mod[,2],edge_hub,edge_nonhub,edge_module_out )
  
  num_matrix[ num_matrix[,1]<hub2mod_interaction_cutoff , 1] = 0
	
	pvalue = apply( num_matrix, 1, function(x) phyper(x[1]-1, x[2], x[3], x[4],lower.tail=F) ) 
	p_hub2mod = rbind(p_hub2mod, pvalue)

}

rownames(p_hub2mod)=hub10_inter
colnames(p_hub2mod)=names(pathway_genelist)

crosstalk=c()
crosstalk$p_hub2mod = p_hub2mod
crosstalk$degree = degree
crosstalk$hub_inter = hub10_inter
crosstalk$hub_intra = hub10_intra

return(crosstalk)

}
