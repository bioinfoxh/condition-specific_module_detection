
setwd("/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_network_analysis")


ppi_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/art_exp_development_pathwayCommon/pathwayCommon_v8/ppi.entrez.human.uniq"

exp_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_series_data/gene_expression_RMA_call_PM1/gene_expression_RMA_call_PM1.txt"

exp_sample_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_series_data/gene_expression_RMA_call_PM1/sample_RMA_call_PM1.txt"


exp_gene = as.matrix( read.table(exp_file,header=T,stringsAsFactors=F ) )
exp_sample = read.table(exp_sample_file,header=T,stringsAsFactors = F)



####### construct co-expression ppi network 


ppi = as.matrix( read.table( ppi_file,sep="\t",header=F )  )

genelist = exp_gene[,1]

ppi_in_exp = ppi[(ppi[,1]%in%genelist) & (ppi[,2]%in%genelist),  ]

pcc = c()
for (i in 1:nrow(ppi_in_exp)){
  print(i)
  temp_pcc = c()
  
  temp_gene1 = ppi_in_exp[i,1]
  temp_gene2 = ppi_in_exp[i,2]
  temp_exp1 = exp_gene[ exp_gene[,1]==temp_gene1 ,2:ncol(exp_gene)]
  temp_exp2 = exp_gene[ exp_gene[,1]==temp_gene2 ,2:ncol(exp_gene)]
  temp_pcc_test = cor.test( temp_exp1, temp_exp2  )
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1
  }else{
    temp_pcc[5] = 0
  }
  pcc = rbind(pcc,abs(temp_pcc) ) 
}

rm( list=ls(pat="temp_") )
rownames(pcc)=NULL
colnames(pcc)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

pcc_sig=pcc[ pcc[,4]<=0.05    ,]
write.table(pcc,"coexp_network_all.txt",sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc_sig,"coexp_network_significant.txt",sep="\t",quote=F,col.names=T,row.names=F)

save.image("network_pcc.RData")




