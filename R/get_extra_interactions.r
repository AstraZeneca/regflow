
#' Get extra interactions
#' 
#' It gets extra interactions from a CSV file. Field names need to be equal to 
#' the standard interaction names as given by 
#' [get.standard.names.in.interactions()]. 
#' 
#' @param extra.interactions.file.path Character string with the name of the 
#'     CSV file with the extra interactions. It may be equal to `NULL`. 
#' 
#' @return A data frame with the extra interactions in standard format. Or 
#'     `NULL` if `extra.interactions.file.path == NULL`. 
#' 
#' @export

get.extra.interactions <- function( extra.interactions.file.path )
{
    if ( is.null( extra.interactions.file.path ) )
        return( NULL )
    
    extra.interactions <- read.csv( extra.interactions.file.path )
    
    stopifnot( 
        all( names( extra.interactions ) == 
            get.standard.names.in.interactions() ) && 
        ! anyNA( extra.interactions ) && 
        anyDuplicated( sprintf( "%s|%s", extra.interactions$source, 
            extra.interactions$target ) ) == 0 && 
        all( extra.interactions$sign %in% c( -1, 0, +1 ) )
    )
    
    extra.interactions
}

