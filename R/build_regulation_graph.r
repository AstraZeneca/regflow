
#' Build a regulation graph
#' 
#' It builds a regulation graph from interactions imported from a database and 
#' processed to standard format. It also accepts a set of extra interactions to 
#' add, or supersede those of the database in case of collision. A special 
#' value of `sign == 0` in an extra interaction is interpreted as removal of 
#' the interaction if it exists. 
#' 
#' @param interactions.data Data frame with interactions, usually obtained with 
#'     [get.interactions()]. 
#' @param extra.interactions Data frame with extra interactions, usually 
#'     obtained with [get.extra.interactions()]. Ignored if `NULL`. 
#' @param cache.file.path Character string with the directory and name of an R 
#'     data file to save and reuse the graph built. 
#' @param cache.ignore Boolean indicating whether to force execution and 
#'     ignore a possible previous cache. 
#' 
#' @return An [igraph] graph object with the largest connected component (in a 
#'     weakly sense) of the regulation graph. 
#' 
#' @export

build.regulation.graph <- function( interactions.data, extra.interactions, 
    cache.file.path, cache.ignore )
{
    if ( ! file.exists( cache.file.path ) || cache.ignore )
    {
        # add, modify, or remove extra interactions
        
        if ( ! is.null( extra.interactions ) )
        {
            interactions.edge.main <- sprintf( "%s|%s", 
                interactions.data$source, interactions.data$target )
            interactions.edge.extra <- sprintf( "%s|%s", 
                extra.interactions$source, extra.interactions$target )
            
            interactions.data.remove.idx <- numeric()
            interactions.extra.add.idx <- numeric()
            
            for ( ie.extra.idx in 1 : length( interactions.edge.extra ) )
            {
                ie.main.idx <- match( interactions.edge.extra[ ie.extra.idx ], 
                    interactions.edge.main, nomatch = 0 )
                
                if ( ie.main.idx != 0 )
                {
                    ie.extra.sign <- extra.interactions$sign[ ie.extra.idx ]
                    
                    if ( ie.extra.sign != 0 )
                    {
                        interactions.data$sign[ ie.main.idx ] <- ie.extra.sign
                        interactions.data$type[ ie.main.idx ] <- 
                            extra.interactions$type[ ie.extra.idx ]
                        interactions.data$evidence[ ie.main.idx ] <- 
                            extra.interactions$evidence[ ie.extra.idx ]
                    }
                    else
                        interactions.data.remove.idx <- 
                            c( interactions.data.remove.idx, ie.main.idx )
                }
                else
                    interactions.extra.add.idx <- c( interactions.extra.add.idx, 
                        ie.extra.idx )
            }
            
            if ( length( interactions.data.remove.idx ) > 0 )
                interactions.data <- 
                    interactions.data[ - interactions.data.remove.idx, ]
            
            if ( length( interactions.extra.add.idx ) > 0 )
                interactions.data <- rbind( interactions.data, 
                    extra.interactions[ interactions.extra.add.idx, ] )
        }
        
        # build graph from interactions
        
        rg <- graph_from_data_frame( interactions.data )
        
        # remove components different from the main one
        
        rg.component <- components( rg, mode = "weak" )
        
        rg <- delete_vertices( rg, which( 
            rg.component$membership != which.max( rg.component$csize ) 
        ) )
        
        # save graph cache
        
        save( rg, file = cache.file.path )
    }
    else
    {
        object.loaded <- load( cache.file.path )
        stopifnot( object.loaded == "rg" )
    }

    rg
}

