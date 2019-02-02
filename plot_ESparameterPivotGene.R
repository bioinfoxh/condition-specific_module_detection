



pathway_ppi = lapply( pathway_genelist, function(x) ppi[ ppi[,1]%in%x & ppi[,2]%in%x, ]  )

setwd("featureModuleGene/ESparameterFeatureGene_figure/")

ESgene_path = "/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/parameter_analysis/gene_parameter_association/GBLasso_inf_2018_allNET_e3k12/featureGenes_selected/featureGene_in_coexpNetwork/"

temp_ESgene = read.table( paste(ESgene_path,"featureGenePivot_TE.txt", sep=""), header=TRUE, stringsAsFactors=FALSE )

if( nrow(temp_ESgene)>0 ) {
  temp_gene = unique( temp_ESgene[,1] )
  
  for(i in 1:length(temp_gene) ){
    temp_pivot = temp_gene[i]
    temp_pivot_symbol = unique(temp_ESgene[ temp_ESgene[,1]==temp_pivot ,2])
    temp_mid = as.character( temp_ESgene[ temp_ESgene[,1]==temp_pivot ,3] )
    
    temp_module_ppi = c()
    temp_module_gene = c()
    for( j in 1:length(temp_mid)){
      temp_module_ppi = rbind( temp_module_ppi, pathway_ppi[ names(pathway_ppi)==temp_mid[j] ][[1]] )
      temp_module_gene = c(temp_module_gene,  pathway_genelist[ names(pathway_genelist)==temp_mid[j] ][[1]] )
    }
    
    temp_module_ppi = unique(temp_module_ppi)
    temp_module_gene = unique(temp_module_gene)
    
    temp_pivot_ppi = ppi_double[ ppi_double[,1]==temp_pivot & ppi_double[,2]%in%temp_module_gene  ,]
    
    temp_dataframe = rbind(temp_module_ppi, temp_pivot_ppi)
    temp_vertice = c(temp_pivot, temp_module_gene)
    temp_verticeShape = c("square", rep("circle",length(temp_module_gene)))
    
    temp_symbol = gene2symbol[ match( temp_vertice, gene2symbol[,1] )  , 2]
    
    library(igraph)
    
    temp_g = graph.data.frame( temp_dataframe[,1:2], directed = FALSE,vertices = temp_vertice  )
    
    V(temp_g)$symbol = temp_symbol
    V(temp_g)$shape = temp_verticeShape
    
    E(temp_g)$weight = abs( temp_dataframe[,3] )
    E(temp_g)$direction = temp_dataframe[,5]
    
    e_lty =  E(temp_g)$direction 
    e_lty[e_lty==0] = 2
    e_lty[e_lty==1] = 1
    
    e_width = E(temp_g)$weight
    e_width = ((e_width)^2) * 2
    
    pdf( paste("TE_pivot_",temp_pivot_symbol,".pdf", sep=""),width = 10, height = 10 )  
    plot(temp_g,
         layout = layout.fruchterman.reingold,
         vertex.shape=V(temp_g)$shape,
         vertex.size=20*10/length(V(temp_g)),
         vertex.color="dark gray",
         vertex.frame.color="dark gray",
         vertex.label=V(temp_g)$symbol,
         vertex.label.dist=0,
         vertex.label.font=3,
         edge.lty=e_lty,
         edge.width=e_width   )
    dev.off()
    
  }
  
  
}






