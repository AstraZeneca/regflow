#' Write table with regulation of intermediate genes in a pathway
#' 
#' Writes a CSV table with the amount of regulation that flows through 
#' intermediate genes or proteins in a pathway. Values are given as absolute 
#' values and relative values with respect to the overall regulation from 
#' `gene.from` to `gene.to`. 
#' 
#' The call is ignored if `gene.from` or `gene.to`, after removing overlapping 
#' values from `gene.to`, is `NULL` or has no elements. 
#' 
#' @param regulation.flow A matrix with all-step regulation flow values, 
#'     usually obtained with [calculate.regulation.matrices()]. 
#' @param gene.from Character vector with gene names as source of regulation. 
#' @param gene.to Character vector with gene names as target of regulation. 
#'     Elements also in `gene.from` are ignored.
#' @param write.dir Character string with the directory to write the table.  
#' @param write.file.name Character string with the table file name. If `NULL`, 
#'     default values are generated depending on the number of genes in 
#'     `gene.from` and `gene.to`. 
#' 
#' @return `NULL`, invisibly. 
#' 
#' @export

write.pathway <- function( regulation.flow, gene.from, gene.to, write.dir, 
    write.file.name = NULL )
{
    gene.from <- correct.genes(gene.from, regulation.flow)
    gene.to <- correct.genes(gene.to, regulation.flow)
    
    gene.to.ex <- setdiff( gene.to, gene.from )
    
    gene.from.n <- length( gene.from )
    gene.to.n <- length( gene.to.ex )
    
    if ( gene.from.n == 0 || gene.to.n == 0 )
        invisible()
    
    if ( gene.from.n > 1 )
        regulation.flow.from <- colMeans( as.matrix( abs( 
            regulation.flow[ gene.from, ] ) ) )
    else
        regulation.flow.from <-  abs( regulation.flow[ gene.from, ] )
    
    if ( gene.to.n > 1 )
        regulation.flow.to <- rowMeans( as.matrix( abs( 
            regulation.flow[ , gene.to.ex ] ) ) )
    else
        regulation.flow.to <- abs( regulation.flow[ , gene.to.ex ] )
    
    regulation.intermediate <- regulation.flow.from * regulation.flow.to
    
    gene.to.net <- gene.to.ex[ regulation.intermediate[ gene.to.ex ] != 0 ]
    gene.to.net.n <- length( gene.to.net )
    
    gene.from.to.net <- union( gene.from, gene.to.net )
    
    regulation.intermediate <- regulation.intermediate[ 
        regulation.intermediate != 0 ]
    
    if ( length( regulation.intermediate ) == 0 )
        return( invisible() )
    
    regulation.intermediate.gene <- names( sort( 
        regulation.intermediate, decreasing = TRUE ) )
    
    regulation.intermediate.gene <- regulation.intermediate.gene[ 
        ! regulation.intermediate.gene %in% gene.from.to.net ]
    
    regulation.intermediate.table <- data.frame( 
        gene = regulation.intermediate.gene, 
        absolute.regulation = 
            regulation.intermediate[ regulation.intermediate.gene ], 
        relative.regulation = 
            regulation.intermediate[ regulation.intermediate.gene ] / 
            mean( regulation.flow.from[ gene.to.net ] )
    )
    
    if ( ! dir.exists( write.dir ) )
        dir.create( write.dir, recursive = TRUE )
    
    if ( ! is.null( write.file.name ) )
        write.pw.file.name <- write.file.name
    else if ( gene.from.n == 1 && gene.to.n == 1 )
        write.pw.file.name <- 
            sprintf( "pathway_%s_2_%s.csv", correct.file.name(gene.from), correct.file.name(gene.to) )
    else if ( gene.from.n == 1 )
        write.pw.file.name <- sprintf( "pathway_%s_2_gene_set.csv", correct.file.name(gene.from) )
    else if ( gene.to.n == 1 )
        write.pw.file.name <- sprintf( "pathway_gene_set_2_%s.csv", correct.file.name(gene.to) )
    else
        write.pw.file.name <- "pathway_gene_set_2_gene_set.csv"
    
    write.pathway.file.path <- file.path( write.dir, write.pw.file.name )
    
    write.csv( regulation.intermediate.table, write.pathway.file.path, 
        row.names = FALSE )
    
    invisible()
}

