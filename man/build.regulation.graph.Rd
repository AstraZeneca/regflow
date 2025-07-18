% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_regulation_graph.r
\name{build.regulation.graph}
\alias{build.regulation.graph}
\title{Build a regulation graph}
\usage{
build.regulation.graph(
  interactions.data,
  extra.interactions,
  cache.file.path,
  cache.ignore
)
}
\arguments{
\item{interactions.data}{Data frame with interactions, usually obtained with
\code{\link[=get.interactions]{get.interactions()}}.}

\item{extra.interactions}{Data frame with extra interactions, usually
obtained with \code{\link[=get.extra.interactions]{get.extra.interactions()}}. Ignored if \code{NULL}.}

\item{cache.file.path}{Character string with the directory and name of an R
data file to save and reuse the graph built.}

\item{cache.ignore}{Boolean indicating whether to force execution and
ignore a possible previous cache.}
}
\value{
An \link{igraph} graph object with the largest connected component (in a
weakly sense) of the regulation graph.
}
\description{
It builds a regulation graph from interactions imported from a database and
processed to standard format. It also accepts a set of extra interactions to
add, or supersede those of the database in case of collision. A special
value of \code{sign == 0} in an extra interaction is interpreted as removal of
the interaction if it exists.
}
