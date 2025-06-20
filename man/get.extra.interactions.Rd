% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_extra_interactions.r
\name{get.extra.interactions}
\alias{get.extra.interactions}
\title{Get extra interactions}
\usage{
get.extra.interactions(extra.interactions.file.path)
}
\arguments{
\item{extra.interactions.file.path}{Character string with the name of the
CSV file with the extra interactions. It may be equal to \code{NULL}.}
}
\value{
A data frame with the extra interactions in standard format. Or
\code{NULL} if \code{extra.interactions.file.path == NULL}.
}
\description{
It gets extra interactions from a CSV file. Field names need to be equal to
the standard interaction names as given by
\code{\link[=get.standard.names.in.interactions]{get.standard.names.in.interactions()}}.
}
