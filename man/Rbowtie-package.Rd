\docType{package}
\name{Rbowtie-package}
\alias{Rbowtie-package}
\alias{Rbowtie}

\title{R bowtie wrapper}

\description{
This package provides an R wrapper around the popular bowtie short read aligner
and around SpliceMap, a de novo splice junction discovery and alignment tool.

The package is used by the \pkg{QuasR} bioconductor package. We recommend to use
\pkg{QuasR} instead of using \code{Rbowtie} directly. \pkg{QuasR} provides a simpler
interface than \code{Rbowtie} and covers the whole analysis workflow of typical
high throughput sequencing experiments, starting from the raw sequence reads,
over pre-processing and alignment, up to quantification.

We also recommend to use \pkg{Rhisat2} instead of SpliceMap for spliced alignment,
which provides a more up-to-date and improved spliced short read aligner.
}

\references{
Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment
of short DNA sequences to the human genome. Genome Biology 10(3):R25 (2009).

Au KF, Jiang H, Lin L, Xing Y, Wong WH. Detection of splice junctions from paired-end
RNA-seq data by SpliceMap. Nucleic Acids Research, 38(14):4570-8 (2010).
}

\seealso{
  \code{\link{bowtie}}, \code{\link{SpliceMap}}
}

\author{
Florian Hahne, Anita Lerch, Michael B Stadler
}
