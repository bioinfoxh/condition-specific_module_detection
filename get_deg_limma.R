


options(scipen=999)

library(limma,quietly = TRUE)

# exp=as.matrix(read.table( "/Users/Hui/myFiles/cambridge/projects/epihealthnet/art_exp_development_pathwayCommon/exp_data_helen/exp_mas5.geneid",sep="\t",header=T   ))

exp=exp_gene

temp1 = exp[,1]
exp = exp[, 2:ncol(exp) ]
rownames(exp) = temp1

label = colnames(exp)
label = sapply(strsplit(label,".",fixed=T), function(x) x[[1]] )


########

design <- model.matrix(~0+factor(label))
colnames(design) <- levels( factor(label) )
fit <- lmFit(exp, design)
# To make all pair-wise comparisons between the three groups the appropriate contrast matrix can be created by
contrast.matrix <- makeContrasts(X4cell-oocyte,
                                 X8cell-X4cell, 
                                 Blast-X8cell,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# A list of top genes differential expressed in group2 versus group1 can be obtained from
deg_0_4 = topTable(fit2, coef=1, adjust="BH",n=Inf)
#deg_0_8 = topTable(fit2, coef=2, adjust="BH",n=Inf)
deg_4_8 = topTable(fit2, coef=2, adjust="BH",n=Inf)
#deg_0_b = topTable(fit2, coef=4, adjust="BH",n=Inf)
#deg_4_b = topTable(fit2, coef=5, adjust="BH",n=Inf)
deg_8_b = topTable(fit2, coef=3, adjust="BH",n=Inf)


# The outcome of each hypothesis test can be assigned using
# results <- decideTests(fit2)

sig_0_4 = deg_0_4[ deg_0_4[,4]<=0.01 ,  ]
#sig_0_8 = deg_0_8[ deg_0_8[,5]<=0.1 ,  ]
#sig_0_b = deg_0_b[ deg_0_b[,5]<=0.1 ,  ]
sig_4_8 = deg_4_8[ deg_4_8[,4]<=0.01 ,  ]
#sig_4_b = deg_4_b[ deg_4_b[,5]<=0.1 ,  ]
sig_8_b = deg_8_b[ deg_8_b[,4]<=0.01 ,  ]


sig_0_4_up = rownames( sig_0_4[ sig_0_4[,1]>0 ,] )
sig_0_4_down = rownames( sig_0_4[sig_0_4[,1]<0,] )
sig_4_8_up = rownames( sig_4_8[ sig_4_8[,1]>0 , ] )
sig_4_8_down = rownames( sig_4_8[ sig_4_8[,1]<0 , ] )
sig_8_b_up = rownames( sig_8_b[ sig_8_b[,1]>0 , ] )
sig_8_b_down = rownames( sig_8_b[ sig_8_b[,1]<0 , ] )


sig_0_4 = data.frame(GeneID=rownames(sig_0_4), sig_0_4, stringsAsFactors = F)
write.table(sig_0_4,"deg_sig_0_4",quote=F,sep="\t",col.names=T,row.names=F)
sig_0_8 = data.frame(GeneID=rownames(sig_0_8), sig_0_8, stringsAsFactors = F)
write.table(sig_0_8,"deg_sig_0_8",quote=F,sep="\t",col.names=T,row.names=F)
sig_0_b = data.frame(GeneID=rownames(sig_0_b), sig_0_b, stringsAsFactors = F)
write.table(sig_0_b,"deg_sig_0_b",quote=F,sep="\t",col.names=T,row.names=F)
sig_4_8 = data.frame(GeneID=rownames(sig_4_8), sig_4_8, stringsAsFactors = F)
write.table(sig_4_8,"deg_sig_4_8",quote=F,sep="\t",col.names=T,row.names=F)
sig_4_b = data.frame(GeneID=rownames(sig_4_b), sig_4_b, stringsAsFactors = F)
write.table(sig_4_b,"deg_sig_4_b",quote=F,sep="\t",col.names=T,row.names=F)
sig_8_b = data.frame(GeneID=rownames(sig_8_b), sig_8_b, stringsAsFactors = F)
write.table(sig_8_b,"deg_sig_8_b",quote=F,sep="\t",col.names=T,row.names=F)









