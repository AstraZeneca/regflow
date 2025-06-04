#' Plot distributions of regulation flow
#' 
#' Plots distributions densities and histograms of regulation flow, separated 
#' in positive and negative values. It overlays gene labels to highlight 
#' particular set of genes. 
#' 
#' The call is ignored if `gene.from` or `gene.to` (depending on `is.from` is 
#' `TRUE` or `FALSE`) is `NULL` or has no elements. 
#' 
#' @param regulation.flow A matrix with regulation flow values, obtained with 
#'     [calculate.regulation.matrices()]. 
#' @param gene.from Character vector with gene names as source of regulation. 
#' @param gene.to Character vector with gene names as target of regulation. 
#' @param is.from Boolean indicating: 
#' * `TRUE`: Downstream regulation from `gene.from`. 
#' * `FALSE`: Upstream regulation to `gene.to`. 
#' @param plot.dir Character string with the directory to save the plot. 
#' @param plot.seed Integer with a seed to randomize the positions of overlaid 
#'     gene labels. 
#' 
#' @return `NULL`, invisibly. 
#' 
#' @export

plot_regulation.flow <- function( regulation.flow, gene.from, gene.to, is.from, 
    plot.dir, plot.seed )
{
    if ( ( is.from && length( gene.from ) == 0 ) || 
         ( ! is.from && length( gene.to ) == 0 ) )
        return( invisible() )
  
    gene.from <- correct.genes(gene.from, regulation.flow)
    gene.to <- correct.genes(gene.to, regulation.flow)
    
    set.seed( plot.seed )
    
    if ( is.from )
    {
        gene.x <- gene.from
        gene.reg <- gene.to
        dir.label <- "from"
        
        gene.n <- length( gene.x )
        if ( gene.n > 1 )
            reg.flow <- apply( abs( regulation.flow[ gene.x, ] ), 2, mean )
        else
            reg.flow <- regulation.flow[ gene.x, ]
    }
    else
    {
        gene.x <- gene.to
        gene.reg <- gene.from
        dir.label <- "to"
        
        gene.n <- length( gene.x )
        if ( gene.n > 1 )
            reg.flow <- apply( abs( regulation.flow[ , gene.x ] ), 1, mean )
        else
            reg.flow <- regulation.flow[ , gene.x ]
    }
    
    if ( gene.n > 1 )
    {
        gene.label <- "Gene Set"
        gene.file.name <- "gene_set"
        
        reg.flow.log <- log10( reg.flow[ reg.flow > 0 ] )
        reg.flow.pos.log <- 0
        reg.flow.neg.log <- 0
    }
    else
    {
        gene.label <- gene.x
        gene.file.name <- correct.file.name(gene.x)
        
        reg.flow.log <- 0
        reg.flow.pos.log <- log10( reg.flow[ reg.flow > 0 ] )
        reg.flow.neg.log <- log10( - reg.flow[ reg.flow < 0 ] )
    }
    
    if ( length( reg.flow.log ) >= 10 || length( reg.flow.pos.log ) >= 10 || 
        length( reg.flow.neg.log ) >= 10 )
    {
        if ( ! dir.exists( plot.dir ) )
            dir.create( plot.dir, recursive = TRUE )
        
        plot.regulation.flow.file.path <- file.path( plot.dir, 
            sprintf( "regulation_flow_%s_%s.pdf", dir.label, gene.file.name ) )

        if ( gene.n > 1 ) {
            pdf( file = plot.regulation.flow.file.path, width = 8, height = 5 )
            par( mar = c(3,3,3,1), mgp = c(1.5,0,0), oma = c(0.5,0,3.5,0), 
                srt = 90  )
        }
        else {
            pdf( file = plot.regulation.flow.file.path, width = 8, height = 7 )
            par( mfrow = c(2,1), mar = c(3,3,3,1), mgp = c(1.5,0,0), 
                oma = c(0,0,3.5,0), srt = 90  )
        }
    }
    
    if ( length( reg.flow.log ) >= 10 )
    {
        plot.title <- sprintf( 
            "average absolute regulation %s %s (%d of %d genes)", 
            dir.label, gene.label, length( reg.flow.log ), length( reg.flow ) 
        )
        plot( density( reg.flow.log ), xlim = c( -8, 0 ), ylim = c( 0 , 1 ),
            xlab = sprintf( "log10( average absolute regulation %s %s )", 
                dir.label, gene.label ), 
            main = plot.title, lwd = 2, cex.lab = 1, cex.axis = 1, 
            cex.main = 1.2 )
        hist( reg.flow.log, breaks = 30, freq = FALSE, col = NULL, add = TRUE )
        for ( gener in gene.reg ) {
            gener.rf <- reg.flow[ gener ]
            if ( gener.rf > 0 )
                text( log10( gener.rf ), runif( 1, 0.05, 0.9 ), gener, 
                    adj = c( 0.5, 0.5 ), col = "red3" )
        }
    }

    if ( length( reg.flow.pos.log ) >= 10 )
    {
        plot.title.pos <- sprintf( 
            "positive regulation %s %s (%d of %d genes)", 
            dir.label, gene.label, length( reg.flow.pos.log ), 
            length( reg.flow ) 
        )
        plot( density( reg.flow.pos.log ), xlim = c( -8, 0 ), 
            ylim = c( 0 , 0.8 ),
            xlab = sprintf( "log10( positive regulation %s %s )", dir.label, 
                gene.label ), 
            main = plot.title.pos, lwd = 2, cex.lab = 1, cex.axis = 1, 
            cex.main = 1.2 )
        hist( reg.flow.pos.log, breaks = 30, freq = FALSE, col = NULL, 
            add = TRUE )
        for ( gener in gene.reg ) {
            gener.rf <- reg.flow[ gener ]
            if ( gener.rf > 0 )
                text( log10( gener.rf ), runif( 1, 0.05, 0.7 ), gener, 
                    adj = c( 0.5, 0.5 ), col = "red3" )
        }
    }

    if ( length( reg.flow.neg.log ) >= 10 )
    {
        plot.title.neg <- sprintf( 
            "negative regulation %s %s (%d of %d genes)", 
            dir.label, gene.label, length( reg.flow.neg.log), 
            length( reg.flow ) )
        plot( density( reg.flow.neg.log ), xlim = c( -8, 0 ), 
            ylim = c( 0, 0.8 ), 
            xlab = sprintf( "log10( - negative regulation %s %s )", dir.label, 
                gene.label ), 
            main = plot.title.neg, lwd = 2, cex.lab = 1, cex.axis = 1, 
            cex.main = 1.2 )
        hist( reg.flow.neg.log, breaks = 30, freq = FALSE, col = NULL, 
            add = TRUE )
        for ( gener in gene.reg ) {
            gener.rf <- reg.flow[ gener ]
            if ( gener.rf < 0 )
                text( log10( - gener.rf ), runif( 1, 0.05, 0.7 ), gener, 
                    adj = c( 0.5, 0.5 ), col = "red3" )
        }
    }
    
    if ( length( reg.flow.log ) >= 10 || length( reg.flow.pos.log ) >= 10 || 
        length( reg.flow.neg.log ) >= 10 )
    {
        title( sprintf( "Regulation flow %s %s", dir.label, gene.label ), 
            outer = TRUE, cex.main = 1.5, cex.sub = 1.2 )
        
        dev.off()
    }
    
    invisible()
}

