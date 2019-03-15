## Rbowtie
[![Build Status (Linux)](https://travis-ci.com/fmicompbio/Rbowtie.svg?branch=master)](https://travis-ci.com/fmicompbio/Rbowtie)
[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/github/fmicompbio/Rbowtie?branch=master&svg=true)](https://ci.appveyor.com/project/fmicompbio/Rbowtie)

The `Rbowtie` R package provides an R interface to the
[`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) short-read aligner by
[Langmead et al. (2009)](http://genomebiology.com/2009/10/3/R25), as well as to
[`SpliceMap`](https://web.stanford.edu/group/wonglab/SpliceMap/) by [Au et al. (2010)](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkq211).
The package contains wrapper functions to create a genome index and to perform
the read alignment to the generated index. These are used by the
[`QuasR`](https://bioconductor.org/packages/QuasR/) bioconductor package.
We recommend to use the `QuasR` package instead of using `Rbowtie` directly.

### Source code
The source code for bowtie v1.2.2_p1 was obtained from [https://github.com/BenLangmead/bowtie/archive/v1.2.2_p1.tar.gz](https://github.com/BenLangmead/bowtie/archive/v1.2.2_p1.tar.gz) on March 13, 2019, and patched to include additional changes available on github since the release (up to master branch commit 58c6ac9 from Sep 5, 2018, and additional unmerged pull request #74, commit c7004c790a73cd34a0e5dffcb31d7c184d3e7250). The folders genomes, reads, doc, indexes and scripts were not included into the package to reduce its size.

