\name{SpliceMap}
\alias{SpliceMap}

\title{
Interface to SpliceMap
}

\description{
The following function can be used to call the SpliceMap binaries.
}

\usage{
SpliceMap(cfg)
}

\arguments{
  
  \item{cfg}{A list containing named elements with SpliceMap parameters
    (see sQuote{Details} or SpliceMap documentation).}
 
}


\details{

  The \code{SpliceMap} function performes the same steps as the
  \code{runSpliceMap} binary from the \emph{SpliceMap} software package,
  but using \code{R} functions to improve compatibility on Windows.

  While the original \emph{SpliceMap} software package is able to call
  different tools to find sub-read alignments, the \code{SpliceMap}
  function from the \code{Rbowtie} package works only with \code{bowtie}
  contained in the package itself. Further modifications from the
  original version include the reporting of unmapped reads at the end of
  the output sam file, and the restriction to a single (pair) of input
  sequence file(s).
  
  The \code{cfg} argument is a list with \emph{SpliceMap} configuration
  parameters that would normally be specified using the \emph{SpliceMap}
  configuration file. Here is a list of supported parameters extracted
  from the sample config file that is distributed with \emph{SpliceMap}:

  \itemize{
    \item \code{genome_dir} (single character value)
    Directory of the chromosome files in FASTA format, or path to a
    single FASTA file containing all chromosomes. If a directory, each chromosome
    can be in a separate file or chromosomes can be concatenated, ie. chr1.fa, chr2.fa, ...

    \item \code{reads_list1} and \code{reads_list2} (single character values)
    These are the two input sequence files. \code{reads_list2} can be
    missing if reads are not paired-end. Note: \code{reads_list1} must
    be the first pair, and pair-reads should be in the
    \dQuote{forward-reverse} format. 

    \item \code{read_format} (single character value)
    Format of the sequencer reads, also make sure reads are  not split
    over multiple lines. Choices are: \code{FASTA}, \code{FASTQ}, \code{RAW}

    \item \code{quality_format} (single character value)
    Format of the quality string if FASTQ is used. Choices are:
    \itemize{
      \item \code{phred-33} Phred base 33 (same as Sanger format)
      \item \code{phred-64} Phred base 64 (same as Illumina 1.3+)
      \item \code{solexa} Format used by solexa machines
    }

    \item \code{outfile} (single character value)
    Name of the output file (spliced alignments in SAM format). If a
    file with this name already exists, \code{SpliceMap} will stop with
    an exception. Note that unmapped reads will be appended to the end
    of the output file (for paired-end experiments, only pairs without
    alignments for any read are considered unmapped).
    
    \item \code{temp_path} (single character value)
    Directory name of the directory that stores temporary files. All
    temporary files will be created in a subfolder of
    \code{temp_path} and will be removed when \code{SpliceMap} finishes
    successfully or failed at an intermediate step.

    \item \code{max_intron} (single integer value)
    Maximum intron size, this is absolute 99th-percentile
    maximum. Introns beyond this size will be ignored. If you don't set
    this, we will assume a mamalian genome (400,000)

    \item \code{min_intron} (single integer value)
    25-th intron size, this is the lower 25th-percentile intron
    size. This is not the smallest size that SpliceMap will search. That
    is about ~25bp. If you don't set this, we will assume a mammalian
    genome (20,000)

    \item \code{max_multi_hit} (single integer value)
    Maximum number of multi-hits. If a 25-mer seed has more than this
    many multi-hits, it will be discarded.  Default is 10.

    \item \code{full_read_length} (single integer value)
    Full read length. SpliceMap will only use the first
    \dQuote{full_read_length} bp for mapping. If the read is shorter
    than \dQuote{full_read_length}, the full read will be used before
    head clip.  If you don't set this parameter, SpliceMap will use as
    many as possible. This is for the case where the reads might have
    N's at the end.  It is always desireable to cut off the N's

    \item \code{head_clip_length} (single integer value)
    Number of bases to clip off the head of the read. This clipping is
    applied after \dQuote{full_read_length}

    \item \code{seed_mismatch} (single integer value)
    Number of mismatches allowed in half-seeding. Choices are
    0,1(default) or 2

    \item \code{read_mismatch} (single integer value)
    Maximum number of mismatches allowed in entire read. No limit on
    value, however SpliceMap can only identify reads with a maximum of 2
    mismatches per 25bp.  Default is 2.

    \item \code{max_clip_allowed} (single integer value)
    Maximum number of bases allowed to be soft clipped from the ends of
    reads during alignment. This is required as mismatches near
    junctions could cause parts of a a read to not map.  Default is 40.

    \item \code{num_chromosome_together} (single integer value)
    Number of chromosomes to process at once, to take advantage of
    multi-core systems. The child processes are created using the
    \code{\link[parallel]{makeCluster}} function from the \pkg{parallel}
    package. This is not threading, so it will take extra
    memory. However, running 2 at a time should be fine on current
    hardware. Default = 2

    \item \code{bowtie_base_dir} (single character value)
    Base of bowtie index, this should be the same genome as the
    chromosome files eg. if you bowtie files are
    \dQuote{genome/hg18/genome.1.ewbt}, ... then your base dir is
    \dQuote{genome/hg18/genome}

    \item \code{num_threads} (single integer values)
    Number of threads to use for \code{bowtie} mapping. Default value is 2

    \item \code{try_hard} (single character value)
    Try hard? Choices are \dQuote{yes} or \dQuote{no}.  Default value is
    \dQuote{yes} (about 15\% slower).

    \item \code{selectSingleHit} (single logical value)
    If \code{TRUE} and multiple alignments are found for a read (pair),
    only a single alignment is selected randomly and reported. This is a
    new paramater only available in the \code{R} version of
    \code{SpliceMap}, but not in the original implementation that could
    report several alignments per read (pair).
  }

  The following parameters are mandatory:
  \itemize{
    \item \code{genome_dir}
    \item \code{reads_list1} (\code{reads_list2} for paired-end experiments)
    \item \code{read_format}
    \item \code{quality_format}
    \item \code{bowtie_base_dir}
    \item \code{outfile}
    \item \code{temp_path}
    \item \code{num_threads}
    \item \code{quality_format}
    \item \code{selectSingleHit}
  }
  
}


\value{
  An invisible character vector of length one with the file name of the
  generated output SAM file.

  An exception is thrown if a failure is detected at one of the many steps.
  
}

\references{
    Au KF, Jiang H, Lin L, Xing Y, Wong WH. \emph{Detection of splice
    junctions from paired-end RNA-seq data by SpliceMap}. Nucleic Acids
    Research, 38(14):4570-8 (2010).
}

\author{
Michael Stadler
}

\seealso{
  \code{\link[parallel]{makeCluster}} from package \pkg{parallel}
}

\examples{
## Building a bowtie index
refDir <- system.file(package="Rbowtie", "samples", "refs")
indexDir <- file.path(tempdir(), "refsIndex")

tmp <- bowtie_build(references=dir(refDir, full=TRUE), outdir=indexDir, prefix="index", force=TRUE)

## Alignments
readsFiles <- system.file(package="Rbowtie", "samples", "reads", "reads.fastq")
samFiles <- file.path(tempdir(), "splicedAlignments.sam")

cfg <- list(genome_dir=refDir,
            reads_list1=readsFiles,
            read_format="FASTQ",
            quality_format="phred-33",
            outfile=samFiles,
            temp_path=tempdir(),
            max_intron=400000,
            min_intron=20000,
            max_multi_hit=10,
            seed_mismatch=1,
            read_mismatch=2,
            num_chromosome_together=2,
            bowtie_base_dir=file.path(indexDir, "index"),
            num_threads=4,
            try_hard="yes",
            selectSingleHit=TRUE)
            
res <- SpliceMap(cfg)
}

\keyword{programming}

