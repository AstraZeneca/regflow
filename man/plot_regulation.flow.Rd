% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_regulation_flow.r
\name{plot_regulation.flow}
\alias{plot_regulation.flow}
\title{Plot distributions of regulation flow}
\usage{
plot_regulation.flow(
  regulation.flow,
  gene.from,
  gene.to,
  is.from,
  plot.dir,
  plot.seed
)
}
\arguments{
\item{regulation.flow}{A matrix with regulation flow values, obtained with
\code{\link[=calculate.regulation.matrices]{calculate.regulation.matrices()}}.}

\item{gene.from}{Character vector with gene names as source of regulation.}

\item{gene.to}{Character vector with gene names as target of regulation.}

\item{is.from}{Boolean indicating:
\itemize{
\item \code{TRUE}: Downstream regulation from \code{gene.from}.
\item \code{FALSE}: Upstream regulation to \code{gene.to}.
}}

\item{plot.dir}{Character string with the directory to save the plot.}

\item{plot.seed}{Integer with a seed to randomize the positions of overlaid
gene labels.}
}
\value{
\code{NULL}, invisibly.
}
\description{
Plots distributions densities and histograms of regulation flow, separated
in positive and negative values. It overlays gene labels to highlight
particular set of genes.
}
\details{
The call is ignored if \code{gene.from} or \code{gene.to} (depending on \code{is.from} is
\code{TRUE} or \code{FALSE}) is \code{NULL} or has no elements.
}
