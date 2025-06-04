
check.gene.complex.omnipath <- function( omnipath.gene.complex, hugo.gene )
{
    complex.split <- strsplit( omnipath.gene.complex, "_" )
    
    sapply( complex.split, function( cs ) 
        all( cs %in% hugo.gene )
    )
}


#' Import database from OmniPath
#' 
#' It downloads the interaction database from the OmniPath server, using the 
#' library [OmnipathR]. It saves to a data file the subset of interactions 
#' involving genes with HUGO names, including complexes. This data file will 
#' usually be read when calling [get.interactions()]. 
#' 
#' @param database.file.path Character string with the directory and name of 
#'     the file to save the selected interactions. 
#' @param only.protein.coding Boolean indicating whether to only select 
#'     interactions involving protein-coding genes.  
#' 
#' @return `NULL`, invisibly. 
#' 
#' @export

import.database.omnipath <- function( database.file.path, only.protein.coding )
{
    database.omnipath.all <- OmnipathR::import_all_interactions()
    
    hugo.gene <- get.hugo.gene( only.protein.coding, 
        dirname( database.file.path ) )
    
    database.omnipath <- database.omnipath.all[ 
        check.gene.complex.omnipath( 
            database.omnipath.all$source_genesymbol, hugo.gene ) & 
        check.gene.complex.omnipath( 
            database.omnipath.all$target_genesymbol, hugo.gene ), 
    ]
    
    save( database.omnipath, file = database.file.path )
    
    invisible()
}

