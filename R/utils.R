correct.genes <- function (gene, reg.matrix){
  genes <- colnames(reg.matrix)
  
  gene <- sapply(gene, function(g) {
    if (g %in% genes ){
      return(g)
    } else {
      pattern <- paste("\\b", g, "\\b", sep = "")
      in_complex <- grepl(pattern, genes)
      
      if(sum(in_complex) == 0) stop ("One or more genes are not present in the knowledge graph")
      
      return(genes[in_complex][1])
    }
    
  })
  
  return(gene)
}

correct.file.name <- function(gene.name){
  if(grepl(":", gene.name)){
    return(paste(strsplit(gene.name, ":")[[1]][1], "complex", sep="_"))
  } else{
    return(gene.name)
  }
}

