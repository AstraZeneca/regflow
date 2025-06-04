
standardize.complex.name.omnipath <- function( cno )
{
    gsub( "_", ":", cno )
}


get.interactions.omnipath <- function( database.file.path )
{
    database.omnipath <- get( load( database.file.path ) )
    
    # correct bug in n_sources and n_references
    
    database.omnipath$n_resources <- NULL
    
    database.omnipath$n_sources <- sapply( 
        strsplit( database.omnipath$sources, ";", TRUE ), 
        function( dos ) {
            dos.len <- length( dos )
            if ( dos.len == 1 && is.na( dos ) )
                0
            else
                dos.len
        }
    )
    
    database.omnipath$n_references <- sapply( 
        strsplit( database.omnipath$references, ";", TRUE ), 
        function( dor ) {
            dor.len <- length( dor )
            if ( dor.len == 1 && is.na( dor ) )
                0
            else
                dor.len
        }
    )
    
    # select interactions:
    #     directed and signed interactions with consensus
    #     transcriptional, post_transcriptional, or post_translational type
    #     at least one reference
    
    interactions.omnipath <- database.omnipath[ 
        database.omnipath$consensus_direction == 1 & 
        xor( database.omnipath$consensus_stimulation == 1, 
             database.omnipath$consensus_inhibition == 1 ) & 
        database.omnipath$type %in% 
            c( "transcriptional", "post_transcriptional", 
               "post_translational" ) & 
        database.omnipath$n_references >= 1, 
    ]
    
    # build interactions data frame with standard names
    
    interactions.omnipath.std <- data.frame( 
        source = standardize.complex.name.omnipath( 
            interactions.omnipath$source_genesymbol ), 
        target = standardize.complex.name.omnipath( 
            interactions.omnipath$target_genesymbol ), 
        sign = +1 * interactions.omnipath$consensus_stimulation + 
               -1 * interactions.omnipath$consensus_inhibition, 
        type = interactions.omnipath$type, 
        evidence = interactions.omnipath$n_references 
    )
    
    # change colliding names in original data frame
    
    names( interactions.omnipath )[ 
        match( c( "source", "target", "type" ), 
            names( interactions.omnipath ) ) ] <- 
        c( "source_uniprot", "target_uniprot", "type_omnipath" )
    
    # return interactions with both standard and original names
    
    cbind( interactions.omnipath.std, interactions.omnipath )
}

