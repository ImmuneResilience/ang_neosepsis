# DESeq2 Helper functions

# Using the group terms, generate any comparison of interest
# @param {S4 object} dds = DESeq2 dds object
# @param {character} terms = "group"
# @param {character} numberator = first term to contrast
# @param {character} denominator = second term to contrast
# This filters for FC +/- 1.5
group.contrast <- function(dds, terms, numerator, denominator, tissue){
  
  resdds <- results(dds, contrast = c(terms, numerator, denominator))
  resOrdered <- resdds[order(resdds$padj),] #sort results by padj value
  resOrdered$ABSLFC <- abs(resOrdered$log2FoldChange) #add ABS log2FoldChange column
  resOrdered$FC <- (sign(resOrdered$log2FoldChange))*(2^(resOrdered$ABSLFC)) #add FC column
  filename1 <- paste0("../results/DESeq_", tissue, "_", numerator, "_vs_", denominator, ".csv")
  #write.csv(as.data.frame(resOrdered), file = filename1) #write results to .csv file
  
  DEgenes <- subset(resOrdered, padj <= 0.05 & ABSLFC >= 0.5849625007) #subset results for only DE genes >= +/- 1.5 FC, <= 0.05 padj
  filename2 <- paste0("../results/", tissue, "_", numerator, "_vs_", denominator, "_DEgenes.csv")
  write.csv(as.data.frame(DEgenes), file = filename2) #write DE genes to .csv file
  
  as.data.frame(DEgenes)
}
