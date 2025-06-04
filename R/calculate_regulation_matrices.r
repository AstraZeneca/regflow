
#' Calculate regulation matrices
#' 
#' It calculates one-step and all-step regulation matrices from a regulation 
#' graph. 
#' 
#' @param rg An [igraph] graph object with the regulation graph, usually 
#'     obtained with [build.regulation.graph()]. 
#' @param cache.file.path Character string with the directory and name of an R 
#'     data file to save and reuse the calculated matrices. 
#' @param cache.ignore Boolean indicating whether to force execution and 
#'     ignore a possible previous cache. 
#' 
#' @return A list with two elements: 
#' * `step`: Matrix with one-step regulation. 
#' * `flow`: Matrix with all-step regulation. 
#' 
#' @export

calculate.regulation.matrices <- function( rg, cache.file.path, cache.ignore, alpha=1 )
{
    if ( ! file.exists( cache.file.path ) || cache.ignore )
    {
        rg.edge.degree <- degree( rg, ends( rg, E( rg ) )[ , 1 ], mode = "out" )
        rg.edge.sign <- E( rg )$sign
        
        E( rg )$weight <- rg.edge.sign / rg.edge.degree
        
        rg.matrix <- as_adjacency_matrix( rg, attr = "weight" )
        
        rg.matrix <- rg.matrix[ order( rownames( rg.matrix ) ), 
            order( colnames( rg.matrix ) ) ]
        
        rg.matrix.lm.eigenval <- eigs( rg.matrix, 1 )$values
        
        rg.matrix.scaled <- alpha*0.99 / Mod( rg.matrix.lm.eigenval ) * rg.matrix
        
        regulation.step <- rg.matrix.scaled
        
        regulation.flow <- solve( Diagonal( nrow( regulation.step ) ) - 
            regulation.step )
        rownames( regulation.flow ) <- rownames( regulation.step )
        colnames( regulation.flow ) <- colnames( regulation.step )
        
        regulation.matrices <- list( step = regulation.step, 
            flow = regulation.flow )
        
        save( regulation.matrices, file = cache.file.path )
    }
    else
    {
        object.loaded <- load( cache.file.path )
        stopifnot( object.loaded == "regulation.matrices" )
    }

    regulation.matrices
}

