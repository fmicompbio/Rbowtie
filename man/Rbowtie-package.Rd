\name{Rbowtie-package}
\alias{Rbowtie-package}
\alias{Rbowtie}
\docType{package}

\title{R bowtie wrapper}

\description{
This package provides an R wrapper around the popular bowtie short read aligner and around SpliceMap, a de novo splice junction discovery and alignment tool.

The package is used by the \code{\link[QuasR]{QuasR}} bioconductor package. We recommend to use the \code{\link[QuasR]{QuasR}} package instead of using Rbowtie directly. The \code{\link[QuasR]{QuasR}} package provides a simpler interface then Rbowtie and covers the whole analysis workflow of typical ultra-high throughput sequencing experiments, starting from the raw sequence reads, over pre-processing and alignment, up to quantification.
}

\details{
See packageDescription('Rbowtie') for package details.
}

\author{
Florian Hahne, Anita Lerch, Michael B Stadler
}

%\references{TODO}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
  \code{\link{bowtie}}, \code{\link{SpliceMap}}
}

\examples{
\dontrun{
    example(bowtie)
    example(SpliceMap)
}
}
