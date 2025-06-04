#' Write table with regulation flow from or to a gene 
#' 
#' Writes a CSV table with the amount of regulation flow from or to a gene or 
#' protein. 
#' 
#' The call is ignored if `gene` is `NULL` or has no elements. 
#' 
#' @param regulation.flow A matrix with regulation flow values, obtained with 
#'     [calculate.regulation.matrices()]. 
#' @param gene Character vector with gene names as source or target of 
#'     regulation. 
#' @param is.from Boolean indicating: 
#' * `TRUE`: Downstream regulation from `gene`. 
#' * `FALSE`: Upstream regulation to `gene`. 
#' @param write.dir Character string with the directory to write the table.  
#' 
#' @return `NULL`, invisibly. 
#' 
#' @export

write.regulation.flow <- function( regulation.flow, gene, is.from, write.dir )
{
    gene.n <- length( gene )
    
    if ( gene.n == 0 )
        return( invisible() )
    
    gene <- correct.genes(gene, regulation.flow)
    
    if ( is.from )
    {
        downstream.gene <- colnames( regulation.flow )
        
        if ( gene.n > 1 )
            effect <- apply( abs( regulation.flow[ gene, ] ), 2, mean )
        else
            effect <- regulation.flow[ gene, ]
        
        regulation.flow.table <- data.frame( downstream.gene, 
            regulation = effect )
    }
    else
    {
        upstream.gene <- rownames( regulation.flow )
        
        if ( gene.n > 1 )
            effect <- apply( abs( regulation.flow[ , gene ] ), 1, mean )
        else
            effect <- regulation.flow[ , gene ]
        
        regulation.flow.table <- data.frame( upstream.gene, 
            regulation = effect )
    }
    
    if ( ! dir.exists( write.dir ) )
        dir.create( write.dir, recursive = TRUE )
    
    write.regulation.flow.file.path <- file.path( write.dir, 
        sprintf( "regulation_flow_%s_%s.csv", 
            ifelse( is.from, "from", "to" ), 
            ifelse( gene.n > 1, "gene_set", correct.file.name(gene) ) ) )
    
    write.csv( regulation.flow.table, write.regulation.flow.file.path, 
        row.names = FALSE )
    
    invisible()
}

