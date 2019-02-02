

get_deg_list <- function( exp, label, level, contrast){
  
  require(limma, quietly = TRUE)

  design <- model.matrix(~0+factor(label, level=level ))
  colnames(design) <- levels( factor(label,level=level) )
  fit <- lmFit(exp, design)
  # To make all pair-wise comparisons between the three groups the appropriate contrast matrix can be created by
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  # A list of top genes differential expressed in group2 versus group1 can be obtained from
  
  deg_list = list()
  for (i in 1:length(contrast) ){
    temp_deg = topTable(fit2, coef=i, adjust="BH",n=Inf)
    
    temp_deg = temp_deg[ match( rownames(exp), rownames(temp_deg) ) ,]
    
    deg_list[[i]] = temp_deg
  }
                       
  names(deg_list) = contrast
                       
  return( deg_list  )
                       
}