## Rbowtie
[![Build Status](https://travis-ci.com/fmicompbio/Rbowtie.svg?branch=master)](https://travis-ci.com/fmicompbio/Rbowtie)

The `Rbowtie` R package provides an R interface to the
[`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) short-read aligner by
[Langmead et al. (2009)](http://genomebiology.com/2009/10/3/R25), as well as to
[`SpliceMap`](https://web.stanford.edu/group/wonglab/SpliceMap/) by [Au et al. (2010)](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkq211).
The package contains wrapper functions to create a genome index and to perform
the read alignment to the generated index. These are used by the
[`QuasR`](https://bioconductor.org/packages/QuasR/) bioconductor package.
We recommend to use the `QuasR` package instead of using `Rbowtie` directly.

### Source code
The source code for bowtie v1.1.1 was obtained from [https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.1](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.1) on February 6, 2014.

