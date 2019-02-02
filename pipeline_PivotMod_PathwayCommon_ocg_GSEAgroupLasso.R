

setwd("/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_network_analysis/OCG/")


################################
### process exp and ppi data, construct the network
################################

get_corR <- function(x){
  n = length(x)/2
  res = cor.test( x[1:n], x[(n+1):(n*2)] )
  corR = res$estimate
  p = res$p.value
  return(c(corR,p))
}



####### load gene expression data (Helen)

ppi_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/art_exp_development_pathwayCommon/pathwayCommon_v8/ppi.entrez.human.uniq"

exp_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_series_data/gene_expression_RMA_call_PM1/gene_expression_RMA_call_PM1.txt"

exp_sample_file="/Users/Hui/myFiles/cambridge/projects/epihealthnet/EmbryoScope/development_series_data/gene_expression_RMA_call_PM1/sample_RMA_call_PM1.txt"

exp_gene = read.table(exp_file,header=T,stringsAsFactors=F )
temp = as.character(exp_gene[,1])
exp_gene = exp_gene[,-1]
rownames(exp_gene)=temp

sampleLabel = read.table(exp_sample_file,header=T,stringsAsFactors = F)
exp_label = sampleLabel[ match( colnames(exp_gene), sampleLabel[,1] ) , 2]




####### construct co-expression ppi network 

ppi =  read.table( ppi_file,sep="\t",header=F,stringsAsFactors = F )  
ppi[,1] = as.character(ppi[,1])
ppi[,2] = as.character(ppi[,2])


ppi_in_exp = ppi[ (ppi[,1]%in%rownames(exp_gene)) & (ppi[,2]%in%rownames(exp_gene)),  ]
ppi = ppi_in_exp
rm(ppi_in_exp)

genelist = unique( c(ppi[,1], ppi[,2]) )
exp_gene = exp_gene[ rownames(exp_gene)%in%genelist  , ]



pcc = c()
nBlock = nrow(ppi) %/% 10000 + 1
for (i in 1:nBlock){
  print(i)
  if(i<nBlock){
    temp_ppi = ppi[ (10000*(i-1)+1):(10000*i) , ]
  }else{
    temp_ppi = ppi[ (10000*(i-1)+1):nrow(ppi) , ]
  }
  temp_exp1 = exp_gene[ match(temp_ppi[,1], rownames(exp_gene)) , ]
  temp_exp2 = exp_gene[ match(temp_ppi[,2], rownames(exp_gene)) , ]
  temp_cor = apply( as.matrix( cbind(temp_exp1,temp_exp2) ),1, function(x) get_corR(x) )
  temp_cor = t(temp_cor)
  pcc = rbind( pcc, temp_cor )
}
ppi_pcc = cbind(ppi, pcc)

temp = rep(1, nrow(ppi_pcc))
temp[ ppi_pcc[,3]<0 ] = 0
ppi_pcc = cbind(ppi_pcc,temp)
rownames(ppi_pcc)=NULL
colnames(ppi_pcc)=c("GeneID1","GeneID2","corR","pvalue","pcc_direction")

ppi = ppi_pcc
rm(ppi_pcc,pcc)


ppi_all = ppi
ppi_sig=ppi[ ppi$pvalue<=0.05    ,]
ppi = ppi_sig
rm(ppi_sig)

genelist = unique( c(ppi[,1],ppi[,2])  )
exp_gene = exp_gene[ rownames(exp_gene) %in% genelist  ,]

#######

library(org.Hs.eg.db, quietly = TRUE)

#uniKeys = keys(org.Hs.eg.db, keytype="ENTREZID")
uniKeys = rownames(exp_gene)
cols = c("SYMBOL")
gene2symbol = select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="ENTREZID")

for (i in 1:nrow(gene2symbol)){
  if( is.na(gene2symbol[i,2]) ){
    gene2symbol[i,2] = gene2symbol[i,1]
  }
}




###########################
### clustering
###########################




########### clustering the co-expression network by using MCL

# system("mcl coexp_network_significant --abc -I 2 -o mcl.output.I2")

library(igraph)



# library(MCL)
# temp_ppi = ppi[,1:2]
# colnames(temp_ppi)=NULL
# ppi_g = graph.data.frame( temp_ppi, directed=F, vertices=genelist)
# ppi_adj = get.adjacency( ppi_g ,sparse=T)
# 
# mcl_I5=mcl(ppi_adj, addLoops = TRUE, expansion = 2, inflation = 5, allow1 = TRUE, max.iter = 100, ESM = FALSE)


#load("/media/hx239/Disk1/projects/epihealthnet/art_exp_development_pathwayCommons/test/co_exp_network_sig.RData")

library(linkcomm)

temp_ppi = ppi[,1:3]
temp_ppi[,3] = abs(temp_ppi[,3])

# t1 = Sys.time()
# 
# lc <- getLinkCommunities(temp_ppi, hcmethod = "single", 
#                          use.all.edges = FALSE,
#                          edglim = 5e4, 
#                          directed = FALSE,
#                          bipartite = FALSE, 
#                          dist = NULL, 
#                          plot = TRUE,
#                          check.duplicates = TRUE, 
#                          removetrivial = TRUE,
#                          verbose = TRUE )
# 
# t2 = Sys.time()


t3 = Sys.time()
ocg <- getOCG.clusters(temp_ppi[,1:2])
t4 = Sys.time()


save.image("clusters_OCG_networkSig.RData")










####### calculated co-expression modularity


ocg_nodecluster = ocg$nodeclusters

temp=table(ocg_nodecluster$cluster)
ocg_clustersize = data.frame( cluster=as.integer(names(temp)), size=as.integer(temp) )

temp_cluster_idx = ocg_clustersize[ ocg_clustersize[,2]>=5 ,1]
ocg_nodecluster_selected = ocg_nodecluster[ ocg_nodecluster[,2]%in%temp_cluster_idx, ] 

cluster_genelist = split( ocg_nodecluster_selected$node, ocg_nodecluster_selected$cluster )
names(cluster_genelist) = 1:length(cluster_genelist)

cluster_genelist = lapply( cluster_genelist, function(x) as.character(x) )



################################

# source("get_module_GSEAscore.R")
# 
# 
# get_geneANOVA <- function(x, exp_label){
#   x=as.numeric(x)
#   temp = data.frame(x=x, stage = exp_label,stringsAsFactors = F)
#   temp1 = summary( aov(x ~ stage, temp) )
#   stat = c(temp1[[1]][1,4], temp1[[1]][1,5])
#   return(stat)
# }
# 
# gene_anova = t( apply( exp_gene, 1, function(x) get_geneANOVA(x,exp_label) ) )


source("get_module_ANOVAenrichment.R")

cluster_anova = get_module_ANOVAenrichment(cluster_genelist, exp_gene, exp_label, permTimes=100)


temp = cbind(cluster_anova$GSEA_score, cluster_anova$GSEA_p)
module_idx = as.numeric( rownames( temp[temp[,1]>0 & temp[,2]<=0.05 ,]))

#temp_p = cluster_anova$GSEA_p[module_idx]
#module_idx = order(temp_p)

module_genelist = cluster_genelist[module_idx]
names(module_genelist) = 1:length(module_genelist)

##################### gruop lasso gel-penalty pathway selection

library(grpregOverlap)

label_level = c("Oocyte", "FourCell", "EightCell","Blastocyst")


temp = exp_label
temp[temp!="Oocyte"]="Other"
Y_oocyte = factor( temp, levels=c("Other","Oocyte") )
temp = exp_label
temp[temp!="FourCell"]="Other"
Y_4cell = factor( temp, levels=c("Other","FourCell") )
temp = exp_label
temp[temp!="EightCell"]="Other"
Y_8cell = factor( temp, levels=c("Other","EightCell") )
temp = exp_label
temp[temp!="Blastocyst"]="Other"
Y_blastocyst = factor( temp, levels=c("Other","Blastocyst") )


X = t(exp_gene)
colnames(X) = paste("g",colnames(X),sep="")

group = lapply( module_genelist,function(x) paste("g",x,sep=""))

# set.seed(123)
# cvfit = cv.grpregOverlap(X, Y_oocyte , group = group, nfolds=5,alpha=0.5,nlambda=1000,family="binomial", penalty="grLasso",trace=T)
# lambda = cvfit$lambda[ which( cvfit$cve==min(cvfit$cve))[1] ]
# fit = grpregOverlap(X, Y_oocyte , group = group, lambda=lambda,family="binomial", penalty="grMCP")



get_cveList <- function( X, Y, group, alphaList, nlambda=1000, nfolds=5, family="binomial", penalty="gel",seed=2018){
  cveList = c()
  for(i in 1:length(alphaList)){
    print( paste("alpha = ",alphaList[i],sep="") )
    temp_alpha = alphaList[i]
    temp_cvfit = cv.grpregOverlap(X, Y , group = group, nfolds=nfolds,alpha=temp_alpha,nlambda=nlambda,family=family, penalty=penalty,seed=seed)
    temp_cve = temp_cvfit$cve
    temp_cvse = temp_cvfit$cvse
    temp_pe = temp_cvfit$pe
    temp_lambda = temp_cvfit$lambda
    temp_res = cbind( rep(temp_alpha,length(temp_lambda)), temp_lambda, temp_cve,temp_cvse,temp_pe )
    cveList = rbind(cveList,temp_res)
  }
  colnames(cveList)=c("alpha","lambda","cve","cvse","pe")
  return(cveList)
}

alphaList = 0.05
cveList_oocyte_1 = get_cveList(X, Y_oocyte, group, alphaList, nlambda=1000, nfolds=5, family="binomial", penalty="gel")
cveList_4cell_1 = get_cveList(X, Y_4cell, group, alphaList, nlambda=1000, nfolds=5, family="binomial", penalty="gel")
cveList_8cell_1 = get_cveList(X, Y_8cell, group, alphaList, nlambda=1000, nfolds=5, family="binomial", penalty="gel")
cveList_blastocyst_1 = get_cveList(X, Y_blastocyst, group, alphaList, nlambda=1000, nfolds=5, family="binomial", penalty="gel")


temp = cveList_oocyte[ cveList_oocyte[,3]== min(cveList_oocyte[,3]) ,]
alpha_oocyte = temp[1]
lambda_oocyte = temp[2]
fit_oocyte = grpregOverlap(X, Y_oocyte, group = group, alpha=alpha_oocyte, lambda=lambda_oocyte, family="binomial", penalty="gel")

temp = cveList_4cell[ cveList_4cell[,3]== min(cveList_4cell[,3]) ,]
alpha_4cell = temp[1]
lambda_4cell = temp[2]
fit_4cell = grpregOverlap(X, Y_4cell, group = group, alpha=alpha_4cell, lambda=lambda_4cell, family="binomial", penalty="gel")


temp = cveList_8cell[ cveList_8cell[,3]== min(cveList_8cell[,3]) ,]
alpha_8cell = temp[1]
lambda_8cell = temp[2]
fit_8cell = grpregOverlap(X, Y_8cell, group = group, alpha=alpha_8cell, lambda=lambda_8cell, family="binomial", penalty="gel")


temp = cveList_blastocyst[ cveList_blastocyst[,3]== min(cveList_blastocyst[,3]) ,]
alpha_blastocyst = temp[1]
lambda_blastocyst = temp[2]
fit_blastocyst = grpregOverlap(X, Y_blastocyst, group = group, alpha=alpha_blastocyst, lambda=lambda_blastocyst, family="binomial", penalty="gel")


sum(fit_oocyte$beta.latent!=0)
sum(fit_4cell$beta.latent!=0)
sum(fit_8cell$beta.latent!=0)
sum(fit_blastocyst$beta.latent!=0)

temp = fit_oocyte$beta.latent
beta_oocyte = temp[temp!=0,]

temp = fit_4cell$beta.latent
beta_4cell = temp[temp!=0,]

temp = fit_8cell$beta.latent
beta_8cell = temp[temp!=0,]

temp = fit_blastocyst$beta.latent
beta_blastocyst = temp[temp!=0,]

save.image("clusters_OCG_networkSig_cluster5_ANOVAsig_grLasso.RData")



temp = beta_oocyte
temp = temp[-1]
temp1 = sapply( strsplit(names(temp),"_"), function(x) x[[1]] )
temp2 = sub("grp","",temp1)
temp3 = sapply( strsplit(names(temp),"_"), function(x) x[[2]] )
temp4 = sub("g","",temp3)
temp5 = gene2symbol[ match(temp4,gene2symbol[,1])  ,2]
featureGene_oocyte = data.frame( gene = temp4,symbol=temp5, module=temp2, beta=temp, stringsAsFactors = FALSE)

temp = beta_4cell
temp = temp[-1]
temp1 = sapply( strsplit(names(temp),"_"), function(x) x[[1]] )
temp2 = sub("grp","",temp1)
temp3 = sapply( strsplit(names(temp),"_"), function(x) x[[2]] )
temp4 = sub("g","",temp3)
temp5 = gene2symbol[ match(temp4,gene2symbol[,1])  ,2]
featureGene_4cell = data.frame( gene = temp4,symbol=temp5, module=temp2, beta=temp, stringsAsFactors = FALSE)

temp = beta_8cell
temp = temp[-1]
temp1 = sapply( strsplit(names(temp),"_"), function(x) x[[1]] )
temp2 = sub("grp","",temp1)
temp3 = sapply( strsplit(names(temp),"_"), function(x) x[[2]] )
temp4 = sub("g","",temp3)
temp5 = gene2symbol[ match(temp4,gene2symbol[,1])  ,2]
featureGene_8cell = data.frame( gene = temp4,symbol=temp5, module=temp2, beta=temp, stringsAsFactors = FALSE)

temp = beta_blastocyst
temp = temp[-1]
temp1 = sapply( strsplit(names(temp),"_"), function(x) x[[1]] )
temp2 = sub("grp","",temp1)
temp3 = sapply( strsplit(names(temp),"_"), function(x) x[[2]] )
temp4 = sub("g","",temp3)
temp5 = gene2symbol[ match(temp4,gene2symbol[,1])  ,2]
featureGene_blastocyst = data.frame( gene = temp4,symbol=temp5, module=temp2, beta=temp, stringsAsFactors = FALSE)

module_oocyte = unique( featureGene_oocyte[,3] )
module_4cell = unique( featureGene_4cell[,3] )
module_8cell = unique( featureGene_8cell[,3] )
module_blastocyst = unique( featureGene_blastocyst[,3] )

module_selected = unique( c(module_oocyte,module_4cell,module_8cell,module_blastocyst))

pathway_genelist = module_genelist[ as.numeric( module_selected ) ]

save.image("clusters_OCG_networkSig_cluster5_ANOVAsig_grLasso.RData")


pathway_size = sapply(pathway_genelist, function(x) length(x) )
temp1 = unlist(pathway_genelist)
temp2 = rep(names(pathway_genelist), pathway_size)
gene2pathway = data.frame( gene=temp1, module=temp2 ,stringsAsFactors = FALSE)


########################
################################
label_level = c("Oocyte", "FourCell", "EightCell","Blastocyst")
limmaContrast = c("FourCell-Oocyte", "EightCell-FourCell", "Blastocyst-EightCell")
permTimes = 1000

pathway_degGSEA = get_module_GSEAscore(pathway_genelist, exp_gene, exp_label, label_level, limmaContrast, permTimes)

temp_score = pathway_degGSEA$GSEA_score
temp_score[temp_score>0] = 1
temp_score[temp_score<0] = -1

temp_p = pathway_degGSEA$GSEA_p
temp_p[temp_p<=0.05] = -1
temp_p[temp_p>0.05] = 0
temp_p = abs(temp_p)

pathway_gsea_flag = temp_score * temp_p



######### get pivot genes

temp1 = ppi
temp2 = ppi[, c(2,1,3:ncol(ppi))]
colnames(temp2) = colnames(ppi)
ppi_double  = rbind(temp1,temp2)

gene_intra = unique(gene2pathway[,1])

temp_ppi_intra = ppi_double[ ppi_double[,1]%in%gene_intra ,]
gene_intra_neighbor = setdiff(temp_ppi_intra[,2], gene_intra  )

interGene2pathway = matrix(0, nrow=length(gene_intra_neighbor), ncol=length(pathway_genelist) )
for (i in 1:length(gene_intra_neighbor)){
  print( paste("gene # ",i,sep="")  )
  temp_gene = gene_intra_neighbor[i]
  temp_ppi = ppi_double[ ppi_double[,1]==temp_gene, ]
  for (j in 1:length(pathway_genelist)){
    temp_module = pathway_genelist[[j]]
    temp_ovlp = length(intersect(temp_ppi[,2], temp_module) )
    interGene2pathway[i,j]=temp_ovlp
  } 
}
rownames(interGene2pathway) = gene_intra_neighbor
colnames(interGene2pathway) = names(pathway_genelist)


pivot = rownames( interGene2pathway )[rowSums( interGene2pathway>=2 )>=2 ]
temp1 = interGene2pathway[ rownames(interGene2pathway)%in%pivot  ,]
temp2 = which( temp1>=2, arr.ind = TRUE, useNames = FALSE)
temp3 = apply( temp2, 1, function(x)temp1[x[1],x[2]] )

pivot2pathway = data.frame( gene= rownames(temp1)[temp2[,1]], module=colnames(temp1)[temp2[,2]], num_ppi=temp3, stringsAsFactors = FALSE )
pivot2pathway = pivot2pathway[ order(pivot2pathway[,1]) ,]

temp = gene2symbol[ match(gene2pathway[,1],gene2symbol[,1]) ,2]
gene2pathway_symbol = data.frame( gene=gene2pathway[,1], symbol=temp, module=gene2pathway[,2] )

temp = gene2symbol[ match(pivot2pathway[,1],gene2symbol[,1]) ,2]
pivot2pathway_symbol = data.frame( gene=pivot2pathway[,1], symbol=temp, module=pivot2pathway[,2], num_ppi=pivot2pathway[,3] )



##########################

dir.create("featureModuleGene")

setwd("featureModuleGene")

write.table( featureGene_oocyte, "featureModuleGene_oocyte.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table( featureGene_4cell, "featureModuleGene_4cell.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table( featureGene_8cell, "featureModuleGene_8cell.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table( featureGene_blastocyst, "featureModuleGene_blastocyst.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

write.table(gene2pathway_symbol, "gene2pathway.txt",quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
write.table(pivot2pathway_symbol, "pivot2pathway_2edges.txt",quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )



source("/Users/Hui/Dropbox/my/Rscripts/enrichmentGO/FDRoptional/enrichmentGO_human_EntrezGeneID.R")

dir.create("pathway_enrichment")
enrichmentGO( pathway_genelist,  category="BP", minDEG=3, cutoff=0.1, FDRcorrection=F, DelRedundency=F, "pathway_enrichment" )




source("../getEnrichmentBarplot.R")

dir.create("pathway_enrichment_figure")

tempFiles = list.files("pathway_enrichment",pattern = ".txt")
for( i in 1:length(tempFiles)){
    temp_inputFile = paste("pathway_enrichment/",tempFiles[i],sep="")
    print(temp_inputFile)
    temp_outputFile = paste("pathway_enrichment_figure/",tempFiles[i],".jpg",sep="")
    getEnrichmentBarplot(temp_inputFile,20,temp_outputFile)
}




# setwd("module_enrichmentGO")
# 
# enriched_module_name = c()
# files = list.files()
# for (i in 1:length(files)){
#   if( file.info(files[i])[1] > 1 ) {
#     temp = unlist(strsplit(files[i],"_"))[4]
#     temp_name = unlist( strsplit(temp,".",fixed = T) )[1]
#     enriched_module_name = c(enriched_module_name,temp_name )  
#   }
# }
# 
# setwd("..")



################## module expression plot



temp_exp_ooc= rowMeans( exp_gene[ , grep( "Oocyte",exp_label )    ]   )
temp_exp_4cell= rowMeans( exp_gene[ , grep( "FourCell",exp_label )    ]   )
temp_exp_8cell= rowMeans( exp_gene[ , grep( "EightCell",exp_label )    ]   )
temp_exp_blast= rowMeans( exp_gene[ , grep( "Blastocyst",exp_label )    ]   )

exp_mean=data.frame( rownames(exp_gene),temp_exp_ooc,temp_exp_4cell,temp_exp_8cell,temp_exp_blast,stringsAsFactors = F )
colnames(exp_mean)=c("GeneID","oocyte","4cell","8cell","blastocyst")

# max_exp=max(exp_mean[,2:ncol(exp_mean)])
# min_exp=min(exp_mean[,2:ncol(exp_mean)])

# enriched_module_id = read.table("enriched_module_id",header=F)$V1



dir.create("pathway_expression_figure")

setwd("pathway_expression_figure")

for (mid in names(pathway_genelist)){ 
  temp1 = pathway_genelist[ names(pathway_genelist)==mid ][[1]]
  #  temp1 = temp1[ !is.na(temp1)  ]
  temp1_exp = exp_mean[ exp_mean[,1]%in%temp1 , ]
  max_exp=max(temp1_exp[,2:ncol(temp1_exp)])
  min_exp=min(temp1_exp[,2:ncol(temp1_exp)])
  
  temp_figure_name=paste("figure_exp_of_module_",mid,".jpg",sep="")
  jpeg(temp_figure_name)
  plot( as.numeric(exp_mean[  exp_mean[,1]==temp1[1] , 2:ncol(exp_mean) ] ) ,type="o",pch=20,col="grey",ylim=c(min_exp,max_exp),xlab="Development Stage",ylab="Gene Expression", xaxt="n")
  axis(side=1,at=1:4,labels=c("oocyte","4cell","8cell","blastocyst"))
  for (i in 2:length(temp1) ){
    lines( as.numeric(exp_mean[  exp_mean[,1]==temp1[i] , 2:ncol(exp_mean) ]) ,type="o",pch=20,col="grey"   )
  }
  dev.off()
}

setwd("..")
###################################

setwd("..")
save.image("clusters_OCG_networkSig_cluster5_ANOVAsig_grLasso.RData")


###############

plot_featureModule.R

#################

temp1= unique(featureGene_oocyte[,3])
temp_oocyte = data.frame( temp1, stage=rep("oocyte",length(temp1)) )
temp1= unique(featureGene_4cell[,3])
temp_4cell = data.frame( temp1, stage=rep("4cell",length(temp1)) )
temp1= unique(featureGene_8cell[,3])
temp_8cell = data.frame( temp1, stage=rep("8cell",length(temp1)) )
temp1= unique(featureGene_blastocyst[,3])
temp_blast = data.frame( temp1, stage=rep("blast",length(temp1)) )

pathway2stage = rbind( temp_oocyte,temp_4cell,temp_8cell,temp_blast )

# gene2pathway_gseaSig = gene2pathway[gene2pathway[,3]%in%pathway_gseaSig,]
# 
# write.table(gene2pathway,"gene2pathway.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
# write.table(gene2pathway_gseaSig,"gene2pathway_gseaSig.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
# 





# cluster5_modularity = sapply( cluster5_nodelist, function(x){mean( temp_ppi[ temp_ppi$GeneID1%in%x & temp_ppi$GeneID2%in%x, 3] )} )

# permTimes = 10000
# cluster5_modularity_perm = c()
# for (i in 1:permTimes){
#   print(paste("modularity permutation times ",i,sep="")  )
#   temp_perm_cor = sample( temp_ppi$corR )
#   temp_perm_ppi = cbind( temp_ppi[,1:2], temp_perm_cor )
#   temp_modularity = sapply( cluster5_nodelist, function(x) mean( temp_perm_ppi[ temp_perm_ppi$GeneID1%in%x & temp_perm_ppi$GeneID2%in%x, 3] ) )
#   cluster5_modularity_perm = cbind(cluster5_modularity_perm, temp_modularity)
# }
# cluster5_modularity_p = apply( cbind(cluster5_modularity,cluster5_modularity_perm),1, function(x){sum(x[2:(permTimes+1)]>=x[1])/permTimes} )

# cluster5_modularity_p=c()
# for(i in 1:nrow(cluster5_size)){
#   print( paste("permutation modularity score for cluster ",i,sep="") )
#   
#   temp_clusterSize = cluster5_size[i,2]
#   temp_modularity = cluster5_modularity[i]
#   
#   if(temp_clusterSize<10){
#     permTimes=1000
#   }else{
#     permTimes = 100
#     }
#   
#   temp_perm_modularity = c()
#   for(j in 1:permTimes){
#     temp_perm_score = mean( sample(temp_ppi$corR, temp_clusterSize )  )
#     temp_perm_modularity = c(temp_perm_modularity, temp_perm_score)
#   }
#   
#   temp_p = sum(temp_perm_modularity>temp_modularity) / permTimes
#   cluster5_modularity_p = c(cluster5_modularity_p, temp_p )
#   
# }
# 
# names(cluster5_modularity_p) = cluster5_size$cluster
# 
# 
# 
# cluster5_sig = names( cluster5_modularity_p[cluster5_modularity_p<=0.05]  )
# 
# cluster5_sig_genelist = cluster5_nodelist[ names(cluster5_nodelist)%in%cluster5_sig   ]
# 
# 
# save.image("clusters_OCG_networkSig_cluster5sig.RData")
# 


######### calculate cluster GSEA enrichment of DEG

# source("get_module_GSEAscore.R")
# 
# label_level = c("Oocyte", "FourCell", "EightCell","Blastocyst")
# limmaContrast = c("FourCell-Oocyte", "EightCell-FourCell", "Blastocyst-EightCell")
# permTimes = 100

# cluster5_sig_gsea_z = get_module_GSEAscore(cluster5_sig_genelist, exp_gene, exp_label, label_level, limmaContrast, permTimes)
# 
# cluster5_gsea_z = get_module_GSEAscore(cluster5_nodelist, exp_gene, exp_label, label_level, limmaContrast, permTimes)
# 
# 
# 
# 
# cluster5_sig_gsea_genelist = cluster5_sig_genelist[ rowSums(cluster5_sig_gsea$GSEA_p<=0.05)>0 ]
# cluster5_gsea_genelist = cluster5_nodelist[ rowSums(cluster5_gsea$GSEA_p<=0.05)>0 ]
# 
# cluster5_sig_gseaZ_genelist = cluster5_sig_genelist[ rowSums(cluster5_sig_gsea_z$GSEA_p<=0.05)>0 ]
# cluster5_gseaZ_genelist = cluster5_nodelist[ rowSums(cluster5_gsea_z$GSEA_p<=0.05)>0 ]
# 
# save.image("clusters_OCG_networkSig_cluster5sig_gsea.RData")
# 
# #############
# 
# 
# source("merge_cluster_jaccard_meetMin.R")
# 
# merged_cluster5_sig_gsea = merge_cluster(cluster5_sig_gsea_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=TRUE)
# merged_cluster5_gsea = merge_cluster(cluster5_gsea_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=TRUE)
# 
# merged_cluster5_sig_gseaZ = merge_cluster(cluster5_sig_gseaZ_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=TRUE)
# merged_cluster5_gseaZ = merge_cluster(cluster5_gseaZ_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=TRUE)
# 
# 
# merged_cluster5_sig_gsea_asc = merge_cluster(cluster5_sig_gsea_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=FALSE)
# merged_cluster5_gsea_asc = merge_cluster(cluster5_gsea_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=FALSE)
# 
# merged_cluster5_sig_gseaZ_asc = merge_cluster(cluster5_sig_gseaZ_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=FALSE)
# merged_cluster5_gseaZ_asc = merge_cluster(cluster5_gseaZ_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=FALSE)
# 
# save.image("clusters_OCG_networkSig_cluster5sig_gsea_merged.RData")


#########


# source("merge_cluster_jaccard_meetMin.R")

#mergedCluster_j10m50_dec = merge_cluster(cluster5_sig_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=TRUE)
#mergedCluster_j10m50_asc = merge_cluster(cluster5_sig_genelist, jaccard_cutoff=0.1, meetMin_cutoff=0.5, clusterSizeDecreasing=FALSE)

#mergedCluster_j20m50_dec = merge_cluster(cluster5_sig_genelist, jaccard_cutoff=0.2, meetMin_cutoff=0.5, clusterSizeDecreasing=TRUE)
#mergedCluster_j20m50_asc = merge_cluster(cluster5_sig_genelist, jaccard_cutoff=0.2, meetMin_cutoff=0.5, clusterSizeDecreasing=FALSE)


#save.image("clusters_OCG_networkSig_cluster5sig_mergedCluster.RData")



# source("merge_cluster_jaccard_asc.R")
# mergedCluster_jacd10p_asc = merge_cluster(cluster3_sig_genelist,0.1  )
# mergedCluster_jacd20p_asc = merge_cluster(cluster3_sig_genelist,0.2  )


# source("merge_cluster_jaccard_dec.R")
# mergedCluster_jacd10p_dec = merge_cluster(cluster3_sig_genelist,0.1  )
# mergedCluster_jacd20p_dec = merge_cluster(cluster5_sig_genelist,0.2  )

# source("merge_cluster_dice_dec.R")
# mergedCluster_dice10p_dec = merge_cluster(cluster5_sig_genelist,0.1  )
# mergedCluster_dice20p_dec = merge_cluster(cluster5_sig_genelist,0.2 )




##############################
######## continue with the mergedCluster_j20m50_asc




pathway_genelist = module_genelist[ names(module_genelist)%in% enriched_module_name   ]

intraGene = unique(unlist(pathway_genelist))
pathway_gene2clust = sapply( pathway_genelist, function(x) intraGene%in%x )
rownames(pathway_gene2clust) = intraGene

temp = which(pathway_gene2clust==TRUE,arr.ind = TRUE,useNames = TRUE)
temp1 = gene2symbol[ match(rownames(temp),gene2symbol[,1]) ,]
gene2pathway = data.frame( temp1, PathwayName = colnames(pathway_gene2clust)[temp[,2]] ,stringsAsFactors = F )

# dir.create("mergedCluster_jacd10p_asc_enrichment")
# 
# temp_gene2clust = mergedCluster_jacd10p_asc$mergedCluster_gene2clust
# temp_genelist = apply( temp_gene2clust, 2, function(x) names(x[x==TRUE])  )
# enrichmentGO( temp_genelist,  category="BP", minDEG=3, cutoff=0.05, FDRcorrection=T, DelRedundency=T, "mergedCluster_jacd10p_asc_enrichment" )
# 
# 
# dir.create("mergedCluster_jacd10p_dec_enrichment")
# 
# temp_gene2clust = mergedCluster_jacd10p_dec$mergedCluster_gene2clust
# temp_genelist = apply( temp_gene2clust, 2, function(x) names(x[x==TRUE])  )
# enrichmentGO( temp_genelist,  category="BP", minDEG=3, cutoff=0.05, FDRcorrection=T, DelRedundency=T, "mergedCluster_jacd10p_dec_enrichment" )


# library("GSA")
# 
# temp_label_group = c("EightCell","Blastocyst")
# temp_exp = exp_gene[ , exp_label %in% temp_label_group ,]
# temp_label = exp_label[ exp_label %in% temp_label_group ]
# 
# temp_label[ grep(temp_label_group[1],temp_label) ] = 1
# temp_label[ grep(temp_label_group[2],temp_label) ] = 2
# temp_label = as.numeric(temp_label)
# 
# temp_gsa3 = GSA(temp_exp, temp_label, module_genelist, rownames(temp_exp), resp.type="Two class unpaired", minsize=1,maxsize=1000,nperms=1000 )
# 
# temp1 = cbind(temp_gsa1$pvalues.lo,temp_gsa1$pvalues.hi)
# temp2 = cbind(temp_gsa2$pvalues.lo,temp_gsa2$pvalues.hi)
# temp3 = cbind(temp_gsa3$pvalues.lo,temp_gsa3$pvalues.hi)


#########################

# network=as.matrix(ppi[,1:2])
# hub_cutoff =quantile( table( c(network[,1],network[,2])), 0.85)
hub_cutoff = 40

source("find_hub_pathway_crosstalk.R")
crosstalk = find_hub_module_crosstalk(network,pathway_genelist,hub_cutoff,3)
p_hub2mod = crosstalk$p_hub2mod

pivot=as.matrix( unique( rownames(which(p_hub2mod<=0.05,arr.ind=T,useNames=T)))   )  
pivot = gene2symbol[ match(pivot, gene2symbol[,1]) , ]

pivot2= as.matrix(names(which(rowSums(p_hub2mod<=0.05)>1) ) )
pivot2 = gene2symbol[ match(pivot2, gene2symbol[,1]) , ]

temp=which(p_hub2mod<=0.05,arr.ind=T,useNames=T)
pivot2mod = cbind( rownames(temp), names(pathway_genelist)[ temp[,2] ]  )
temp1 = gene2symbol[ match(pivot2mod[,1], gene2symbol[,1]) , ]
pivot2mod = data.frame( temp1, pivot2mod[,2] )
colnames(pivot2mod)[3] = "PathwayName"
# pivot2mod = t( apply(temp1, 1 , function(temp1) c(temp1[1],colnames(p_hub2mod)[temp1[2]] ) )  )

pivot2mod_pivot2=pivot2mod[  is.element( pivot2mod[,1],pivot2[,1]  ) ,]
# pivot2mod_pivot2=temp[order(temp[,1]),]


pivot2mod_pathway_gseaSig = pivot2mod[ pivot2mod[,3]%in%pathway_gseaSig  ,]

temp = table( pivot2mod_pivot2[ pivot2mod_pivot2[,3]%in%pathway_gseaSig  , 1] )
temp1 = names(temp)[temp>1]
pivot2mod_pivot2_pathway_gseaSig = pivot2mod_pivot2[ pivot2mod_pivot2[,3]%in%pathway_gseaSig & pivot2mod_pivot2[,1]%in%temp1 ,]
pivot2mod_pivot2_pathway_gseaSig = pivot2mod_pivot2_pathway_gseaSig[ order(pivot2mod_pivot2_pathway_gseaSig[,1]) , ]



multiFunGene2clust = pathway_gene2clust[ rowSums(pathway_gene2clust)>1, ]
temp = which(multiFunGene2clust==TRUE, arr.ind = TRUE, useNames = TRUE)
multiFunGene = data.frame(GeneID=rownames(temp), PathwayName=colnames(pathway_gene2clust)[temp[,2]] ,stringsAsFactors = F )
temp1 = gene2symbol[ match(multiFunGene[,1], gene2symbol[,1]) ,]
multiFunGene = data.frame(temp1, PathwayName=multiFunGene[,2])

temp = multiFunGene[ multiFunGene[,3]%in%pathway_gseaSig  ,]
temp1=table(multiFunGene[,1])
temp2 = names(temp1)[temp1>1]
multiFunGene_pathway_gseaSig = temp[ temp[,1]%in%temp2 ,]

dir.create("crosstalk_hubcutoff_40")

setwd("crosstalk_hubcutoff_40")
write.table(pivot,"pivot.hub10.ppi3.GeneID.symbol.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(pivot2,"pivot2.hub10.ppi3.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
write.table(pivot2mod,"pivot2mod.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
write.table(pivot2mod_pivot2,"pivot2mod_pivot2.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
write.table(pivot2mod_pathway_gseaSig,"pivot2mod_pathway_gseaSig.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
write.table(pivot2mod_pivot2_pathway_gseaSig,"pivot2mod_pathway_gseaSig_pivot2.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
write.table(multiFunGene,"multiFunGene.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
write.table(multiFunGene_pathway_gseaSig,"multiFunGene_pathway_gseaSig.GeneID.symbol.txt", quote=F, row.names=F, col.names=F,sep="\t")
setwd("..")

#######################
######### result visualization

####################################  enrichment barplot

dir.create("figure_enrichment")

setwd("figure_enrichment")

for (i in names(pathway_genelist)){
  print(i)
  inputfile=paste( "../module_enrichmentGO/enriched_GOterm_list_",i,".txt",sep=""  )
 
    temp1=read.delim(inputfile,sep="\t", header=T)
    outputfile = paste("enrichment_module_",i,".jpg",sep="")
    jpeg(outputfile,width=800,height=480)
    par(mar=c(10, 42, 10, 0.5) )
    if (nrow(temp1)<=10){
      barplot(rev(-log(temp1$Pvalues)),main="",xlab="-lg(Pvalue)",ylab="",horiz=TRUE,names.arg=rev(temp1$GOtermName),las=1 ,space=1,col="steel blue",width=rep(1,nrow(temp1)))
    }else{
      barplot(rev(-log(temp1$FDR[1:10])),main="",xlab="-lg(Pvalue)",ylab="",horiz=TRUE,names.arg=rev(temp1$GOtermName[1:10]),las=1 ,space=1,col="steel blue",width=rep(1,nrow(temp1)))
    }
    dev.off()
  
}

setwd("..")




########## get ppi symbol, output ppi for each module

temp1 = gene2symbol[match(ppi[,1],gene2symbol[,1]), 2]
temp2 = gene2symbol[match(ppi[,2],gene2symbol[,1]), 2]

ppi_symbol = data.frame( GeneID_1=ppi[,1], Symbol_1=temp1, GeneID_2=ppi[,2], Symbol_2=temp2, ppi[3:ncol(ppi)],stringsAsFactors = FALSE)

dir.create("module_links") 

setwd("module_links")

for( i in names(pathway_genelist) ){
  temp_gene = pathway_genelist[ names(pathway_genelist)==i ][[1]]
  temp_link = ppi_symbol[ ppi_symbol[,1]%in%temp_gene & ppi_symbol[,3]%in%temp_gene  , ]
  filename=paste("ppi_module_",i,".txt",sep="")
  write.table(temp_link,filename,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
}
setwd("..")

# module_ge10_ppi = c()
# for (i in 1:nrow(module_ge10)){
#   temp_md_gene = as.numeric(module_ge10[i, !is.na( module_ge10[i,]  )])
#   temp_md_ppi = pcc[ is.element( pcc[,1],temp_md_gene ) & is.element( pcc[,2], temp_md_gene ),   ]
#   temp_index= rep( i, nrow(temp_md_ppi) )  
#   temp_md_ppi = cbind(temp_md_ppi, temp_index)
#   filename=paste("ppi_module.",i,sep="")
#   write.table(temp_md_ppi,filename,quote=F,sep="\t",row.names=F,col.names=F)
#   #  module_ge10_ppi= rbind( module_ge10_ppi, cbind(temp_md_ppi, temp_index)   )
# }





