## Rbowtie

The `Rbowtie` R package provides an R interface to the
[`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) short-read aligner by
[Langmead et al. (2009)](http://genomebiology.com/2009/10/3/R25), as well as to
[`SpliceMap`](https://web.stanford.edu/group/wonglab/SpliceMap/) by [Au et al. (2010)](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkq211).
The package contains wrapper functions to create a genome index and to perform
the read alignment to the generated index. These are used by the
[`QuasR`](https://bioconductor.org/packages/QuasR/) bioconductor package.
We recommend to use the `QuasR` package instead of using `Rbowtie` directly.

### Source code

The source code for bowtie v1.3.0 was obtained from [https://github.com/BenLangmead/bowtie/archive/v1.3.0.tar.gz](https://github.com/BenLangmead/bowtie/archive/v1.3.0.tar.gz) on November 2, 2020. The folders genomes, reads, doc, indexes and scripts were not included into the package to reduce its size.

### Software status

| Platforms        |  OS              | R CMD check      | Coverage         |
|:----------------:|:----------------:|:----------------:|:----------------:|
| GitHub Actions | Linux/Windows/macOS | [![R build status](https://github.com/fmicompbio/Rbowtie/workflows/R-CMD-check/badge.svg)](https://github.com/fmicompbio/Rbowtie/actions) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/Rbowtie/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/Rbowtie) |
| Bioc ([_devel_](http://bioconductor.org/packages/devel/bioc/html/Rbowtie.html)) | Multiple | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/Rbowtie.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/Rbowtie) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/Rbowtie/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/Rbowtie) |
| Bioc ([_release_](http://bioconductor.org/packages/release/bioc/html/Rbowtie.html)) | Multiple | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/Rbowtie.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/Rbowtie) | [![Codecov.io coverage status](https://codecov.io/github/fmicompbio/Rbowtie/coverage.svg?branch=master)](https://codecov.io/github/fmicompbio/Rbowtie) |
