% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\docType{methods}
\name{plot-SlingshotDataSet}
\alias{plot-SlingshotDataSet}
\alias{plot,SlingshotDataSet,ANY-method}
\alias{lines,SlingshotDataSet-method}
\title{Plot Slingshot output}
\usage{
\S4method{plot}{SlingshotDataSet,ANY}(x, type = NULL,
  show.constraints = FALSE, add = FALSE, dims = seq_len(2), asp = 1,
  cex = 2, lwd = 2, col = 1, ...)

\S4method{lines}{SlingshotDataSet}(x, type = NULL, dims = seq_len(2), ...)
}
\arguments{
\item{x}{a \code{SlingshotDataSet} with results to be plotted.}

\item{type}{character, the type of output to be plotted, can be one of 
\code{"lineages"}, \code{"curves"}, or \code{"both"} (by partial matching),
see Details for more.}

\item{show.constraints}{logical, whether or not the user-specified initial 
and terminal clusters should be specially denoted by green and red dots,
respectively.}

\item{add}{logical, indicates whether the output should be added to an
existing plot.}

\item{dims}{numeric, which dimensions to plot (default is \code{1:2}).}

\item{asp}{numeric, the y/x aspect ratio, see \code{\link{plot.window}}.}

\item{cex}{numeric, amount by which points should be magnified, see
\code{\link{par}}.}

\item{lwd}{numeric, the line width, see \code{\link{par}}.}

\item{col}{character or numeric, color(s) for lines, see \code{\link{par}}.}

\item{...}{additional parameters to be passed to \code{\link{lines}}.}
}
\value{
returns \code{NULL}.
}
\description{
Tools for visualizing lineages inferred by \code{slingshot}.
}
\details{
If \code{type == 'lineages'}, straight line connectors between
  cluster centers will be plotted. If \code{type == 'curves'}, simultaneous
  principal curves will be plotted.

When \code{type} is not specified, the function will first check the
  \code{curves} slot and plot the curves, if present. Otherwise,
  \code{lineages} will be plotted, if present.
}
\examples{
data("slingshotExample")
sds <- slingshot(rd, cl, start.clus = "1")
plot(sds, type = 'b')

# add to existing plot
plot(rd, col = 'grey50')
lines(sds, lwd = 3)

}
