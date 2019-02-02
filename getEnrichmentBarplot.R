getEnrichmentBarplot <- function(enrichmentResultFile, topN, outputFileName){
  if( file.info(enrichmentResultFile)$size > 1   ){
    temp1=read.delim(enrichmentResultFile,sep="\t", header=T, stringsAsFactors = F)
    #outputfile = paste("enrichment_module_",i,".jpg",sep="")
    jpeg(outputFileName,width=800,height=800)
    par(mar=c(10, 42, 10, 0.5) )
    if (nrow(temp1)<=topN){
      barplot(rev(-log(temp1$FDR)),main="",xlab="- log ( p.value )",ylab="",horiz=TRUE,names.arg=rev(temp1$GOtermName),las=1 ,space=1,col="light blue",width=rep(1,nrow(temp1)))
    }else{
      barplot(rev(-log(temp1$FDR[1:topN])),main="",xlab="- log ( p.value )",ylab="",horiz=TRUE,names.arg=rev(temp1$GOtermName[1:topN]),las=1 ,space=1,col="light blue",width=rep(1,nrow(temp1)))
    }
    dev.off()
  }
}
