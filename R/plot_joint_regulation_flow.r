#' Plot joint values of regulation flow for two genes
#' 
#' It obtains scatter plots of regulation flow values, separated by synergistic 
#' and antagonistic effects, positive or negative. It overlays gene labels to 
#' highlight particular set of genes. 
#' 
#' The call is ignored if `gene.from` or `gene.to` (depending on `is.from` is 
#' `TRUE` or `FALSE`) is `NULL` or with length different from 2. 
#' 
#' @param regulation.flow A matrix with regulation flow values, obtained with 
#'     [calculate.regulation.matrices()]. 
#' @param gene.from Character vector with gene names as source of regulation. 
#' @param gene.to Character vector with gene names as target of regulation. 
#' @param is.from Boolean indicating: 
#' * `TRUE`: Downstream regulation from the two genes in `gene.from`. 
#' * `FALSE`: Upstream regulation to the two genes in `gene.to`. 
#' @param plot.dir Character string with the directory to save the plot. 
#' 
#' @return `NULL`. 
#' 
#' @export

plot_joint.regulation.flow <- function( regulation.flow, gene.from, gene.to, 
    is.from, plot.dir )
{
    if ( ( is.from && 
            ( is.null( gene.from ) || length( gene.from ) != 2 ) ) || 
         ( ! is.from && 
            ( is.null( gene.to ) || length( gene.to ) != 2 ) ) )
        return()
        
    if ( is.from ) {
        gene.a <- gene.from[ 1 ]
        gene.b <- gene.from[ 2 ]
        gene.reg <- gene.to
        reg.flow.a <- regulation.flow[ gene.a, ]
        reg.flow.b <- regulation.flow[ gene.b, ]
        reg.flow.gene <- colnames( regulation.flow )
        dir.label <- "from"
    }
    else {
        gene.a <- gene.to[ 1 ]
        gene.b <- gene.to[ 2 ]
        gene.reg <- gene.from
        reg.flow.a <- regulation.flow[ , gene.a ]
        reg.flow.b <- regulation.flow[ , gene.b ]
        reg.flow.gene <- rownames( regulation.flow )
        dir.label <- "to"
    }
    
    reg.flow.q1.bool <- reg.flow.a < 0 & reg.flow.b > 0
    reg.flow.q1.a <- log10( - reg.flow.a[ reg.flow.q1.bool ] )
    reg.flow.q1.b <- log10( reg.flow.b[ reg.flow.q1.bool ] )
    reg.flow.q1.gene <- intersect( gene.reg, reg.flow.gene[ reg.flow.q1.bool ] )
    
    reg.flow.q2.bool <- reg.flow.a > 0 & reg.flow.b > 0
    reg.flow.q2.a <- log10( reg.flow.a[ reg.flow.q2.bool ] )
    reg.flow.q2.b <- log10( reg.flow.b[ reg.flow.q2.bool ] )
    reg.flow.q2.gene <- intersect( gene.reg, reg.flow.gene[ reg.flow.q2.bool ] )
    
    reg.flow.q3.bool <- reg.flow.a < 0 & reg.flow.b < 0
    reg.flow.q3.a <- log10( - reg.flow.a[ reg.flow.q3.bool ] )
    reg.flow.q3.b <- log10( - reg.flow.b[ reg.flow.q3.bool ] )
    reg.flow.q3.gene <- intersect( gene.reg, reg.flow.gene[ reg.flow.q3.bool ] )
    
    reg.flow.q4.bool <- reg.flow.a > 0 & reg.flow.b < 0
    reg.flow.q4.a <- log10( reg.flow.a[ reg.flow.q4.bool ] )
    reg.flow.q4.b <- log10( - reg.flow.b[ reg.flow.q4.bool ] )
    reg.flow.q4.gene <- intersect( gene.reg, reg.flow.gene[ reg.flow.q4.bool ] )
    
    if ( ! dir.exists( plot.dir ) )
        dir.create( plot.dir, recursive = TRUE )
    
    plot.joint.regulation.flow.file.path <- file.path( plot.dir, 
        sprintf( "joint_regulation_flow_%s_%s_and_%s.pdf", 
            dir.label, gene.a, gene.b ) )
    
    pdf( file = plot.joint.regulation.flow.file.path, width = 8, height = 7 )
    
    par( mfrow = c(2,2), mar = c(4,5.5,2,3), mgp = c(2,0.5,0), 
        oma = c(0,0,3,0) )
    
    plot.title.q1 <- sprintf( "%s pos. & %s neg. reg. (%d of %d gen.)", 
        gene.b, gene.a, sum( reg.flow.q1.bool ), length( reg.flow.gene ) )
    plot( reg.flow.q1.a, reg.flow.q1.b, pch = 20, col = "gray60", 
        xlim = c(-9,0.2), ylim = c(-9,0.2), 
        xlab = sprintf( "log10( - negative regulation %s %s )", dir.label, 
            gene.a ), 
        ylab = sprintf( "log10( positive regulation %s %s )", dir.label, 
            gene.b ), 
        main = plot.title.q1 )
    abline( a = 0, b = 1, lty = 2 )
    if ( length( reg.flow.q1.gene > 0 ) )
        text( reg.flow.q1.a[ reg.flow.q1.gene ], 
            reg.flow.q1.b[ reg.flow.q1.gene ], reg.flow.q1.gene, col = "red3", 
            adj = 0.5 )
    mtext( "antagonistic", 2, line = 4, adj = 0.5 )
    
    plot.title.q2 <- sprintf( "%s pos. & %s pos. reg. (%d of %d gen.)", 
        gene.b, gene.a, sum( reg.flow.q2.bool ), length( reg.flow.gene ) )
    plot( reg.flow.q2.a, reg.flow.q2.b, pch = 20, col = "gray60", 
        xlim = c(-9,0.2), ylim = c(-9,0.2), 
        xlab = sprintf( "log10( positive regulation %s %s )", dir.label, 
            gene.a ), 
        ylab = sprintf( "log10( positive regulation %s %s )", dir.label, 
            gene.b ), 
        main = plot.title.q2 )
    abline( a = 0, b = 1, lty = 2 )
    if ( length( reg.flow.q2.gene > 0 ) )
        text( reg.flow.q2.a[ reg.flow.q2.gene ], 
            reg.flow.q2.b[ reg.flow.q2.gene ], reg.flow.q2.gene, col = "red3", 
            adj = 0.5 )
    mtext( "synergistic", 4, line = 1, adj = 0.5 )

    plot.title.q3 <- sprintf( "%s neg. & %s neg. reg. (%d of %d gen.)", 
        gene.b, gene.a, sum( reg.flow.q3.bool ), length( reg.flow.gene ) )
    plot( reg.flow.q3.a, reg.flow.q3.b, pch = 20, col = "gray60", 
        xlim = c(-9,0.2), ylim = c(-9,0.2), 
        xlab = sprintf( "log10( - negative regulation %s %s )", dir.label, 
            gene.a ), 
        ylab = sprintf( "log10( - negative regulation %s %s )", dir.label, 
            gene.b ), 
        main = plot.title.q3 )
    abline( a = 0, b = 1, lty = 2 )
    if ( length( reg.flow.q3.gene > 0 ) )
        text( reg.flow.q3.a[ reg.flow.q3.gene ], 
            reg.flow.q3.b[ reg.flow.q3.gene ], reg.flow.q3.gene, col = "red3", 
            adj = 0.5 )
    mtext( "synergistic", 2, line = 4, adj = 0.5 )
    
    plot.title.q4 <- sprintf( "%s neg. & %s pos. reg. (%d of %d gen.)", 
        gene.b, gene.a, sum( reg.flow.q4.bool ), length( reg.flow.gene ) )
    plot( reg.flow.q4.a, reg.flow.q4.b, pch = 20, col = "gray60", 
        xlim = c(-9,0.2), ylim = c(-9,0.2), 
        xlab = sprintf( "log10( positive regulation %s %s )", dir.label, 
            gene.a ), 
        ylab = sprintf( "log10( - negative regulation %s %s )", dir.label, 
            gene.b ), 
        main =  plot.title.q4 )
    abline( a = 0, b = 1, lty = 2 )
    if ( length( reg.flow.q4.gene > 0 ) )
        text( reg.flow.q4.a[ reg.flow.q4.gene ], 
            reg.flow.q4.b[ reg.flow.q4.gene ], reg.flow.q4.gene, col = "red3", 
            adj = 0.5 )
    mtext( "antagonistic", 4, line = 1, adj = 0.5 )
    
    title( sprintf( "Joint regulation flow %s %s and %s", dir.label, gene.a, 
        gene.b ), outer = TRUE, cex.main = 1.5, cex.sub = 1.2 )
    
    dev.off()
    
    return()
}

