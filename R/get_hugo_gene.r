
get.hugo.gene <- function( only.protein.coding, cache.dir )
{
    cache.file.name <- sprintf( "hugo_gene%s.rda", 
        ifelse( only.protein.coding, "_only_protein_coding", "" ) )
    
    cache.file.path <- file.path( cache.dir, cache.file.name )
    
    if ( ! file.exists( cache.file.path ) )
    {
        ensembl.mart <- biomaRt::useEnsembl( biomart = "genes", 
            dataset = "hsapiens_gene_ensembl" )
        
        hugo.gene <- biomaRt::getBM( 
            attributes = c( "external_gene_name", "gene_biotype" ), 
            filters = ifelse( only.protein.coding, "biotype", "" ), 
            values = ifelse( only.protein.coding, "protein_coding", "" ), 
            mart = ensembl.mart 
        )$external_gene_name
        
        stopifnot( anyDuplicated( hugo.gene ) == 0 )
        
        save( hugo.gene, file = cache.file.path )
    }
    else
    {
        object.loaded <- load( cache.file.path )
        stopifnot( object.loaded == "hugo.gene" )
    }

    hugo.gene
}

