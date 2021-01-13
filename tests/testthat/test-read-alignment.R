context("read-alignment")

test_that("malformed input gives error", {
    tmp <- tempdir()
    outdir <- file.path(tmp, "bowtieindexdir")
    bowtie_build(
        references=system.file("extdata/refs/chr1.fa", package="Rbowtie"),
        outdir=outdir, force=TRUE, execute=TRUE
    )

    ## malformed sequences
    expect_error(bowtie(sequences=TRUE, outfile="file1",
                        index=file.path(outdir, "index")))

    ## malformed outfile
    expect_error(bowtie(sequences=system.file("extdata/reads/reads1.fastq",
                                              package="Rbowtie"),
                        outfile=c("file1", "file2"),
                        index=file.path(outdir, "index")))

    ## malformed sequences, c=TRUE
    expect_error(bowtie(sequences=list(TRUE, TRUE), outfile="file1", c=TRUE,
                        type="paired", index=file.path(outdir, "index")))

    ## malformed index
    expect_error(bowtie(sequences=system.file("extdata/reads/reads1.fastq",
                                              package="Rbowtie"),
                        outfile=file.path(tmp, "alignments.sam"),
                        index=1, force=TRUE,
                        execute=FALSE, type="single"))

    ## outfile is not overwritten if it exists and force=FALSE
    write.table(1, file=file.path(tmp, "alignments2.sam"))
    expect_error(bowtie(sequences=system.file("extdata/reads/reads1.fastq",
                                              package="Rbowtie"),
                        outfile=file.path(tmp, "alignments2.sam"),
                        index=file.path(outdir, "index"), force=FALSE,
                        execute=FALSE, type="single", strict=TRUE))
})

test_that("correctly formatted input works", {
    tmp <- tempdir()
    outdir <- file.path(tmp, "bowtieindexdir")
    bowtie_build(
        references=system.file("extdata/refs/chr1.fa", package="Rbowtie"),
        outdir=outdir, force=TRUE, execute=TRUE
    )
    ## single-end reads
    expect_equal(bowtie(sequences=system.file("extdata/reads/reads1.fastq",
                                              package="Rbowtie"),
                        outfile=file.path(tmp, "alignments.sam"),
                        index=file.path(outdir, "index"), force=TRUE,
                        execute=FALSE, type="single"),
                 paste0(shQuote(file.path(system.file(package="Rbowtie"), "bowtie")),
                        " -x ", shQuote(file.path(outdir, "index")),
                        " ", shQuote(system.file("extdata/reads/reads1.fastq",
                                                 package="Rbowtie")),
                        "  ", shQuote(file.path(tmp, "alignments.sam"))))

    ## paired-end reads
    expect_equal(bowtie(sequences=list(system.file("extdata/reads/reads1.fastq",
                                                   package="Rbowtie"),
                                       system.file("extdata/reads/reads2.fastq",
                                                   package="Rbowtie")),
                        outfile=file.path(tmp, "alignments.sam"),
                        index=file.path(outdir, "index"), force=TRUE,
                        execute=FALSE, type="paired"),
                 paste0(shQuote(file.path(system.file(package="Rbowtie"), "bowtie")),
                        " -x ", shQuote(file.path(outdir, "index")),
                        " -1 ", shQuote(system.file("extdata/reads/reads1.fastq",
                                                    package="Rbowtie")),
                        " -2 ", shQuote(system.file("extdata/reads/reads2.fastq",
                                                    package="Rbowtie")),
                        "   ", shQuote(file.path(tmp, "alignments.sam"))))

    ## single-end reads but with type="paired" should not work
    expect_error(bowtie(sequences=system.file("extdata/reads/reads1.fastq",
                                              package="Rbowtie"),
                        outfile=file.path(tmp, "alignments.sam"),
                        index=file.path(outdir, "index"), force=TRUE,
                        execute=FALSE, type="paired"))

    ## paired-end reads but with type="single" should not work
    expect_error(bowtie(sequences=list(system.file("extdata/reads/reads1.fastq",
                                                   package="Rbowtie"),
                                       system.file("extdata/reads/reads2.fastq",
                                                   package="Rbowtie")),
                        outfile=file.path(tmp, "alignments.sam"),
                        index=file.path(outdir, "index"), force=TRUE,
                        execute=FALSE, type="single"))

    ## missing outfile (return alignments as vector)
    v <- bowtie(sequences=system.file("extdata/reads/reads1.fastq",
                                      package="Rbowtie"),
                index=file.path(outdir, "index"))
    expect_is(v, "character")
    expect_length(v, 2L)
    
    ## bowtie alignments
    r1 <- system.file("extdata", "reads", "reads1.fastq", package = "Rbowtie")
    r2 <- system.file("extdata", "reads", "reads2.fastq", package = "Rbowtie")
    idx <- paste0(outdir, "/index")
    # ... quality mode
    expect_true(is(res1q <- bowtie(sequences = r1, index = idx), "character"))
    expect_length(res1q, 2L)
    expect_true(is(res2q <- bowtie(sequences = list(r1, r2), index = idx, type = "paired"), "character"))
    expect_length(res2q, 0L)
    # ... mismatch mode 
    expect_true(is(res1v <- bowtie(sequences = r1, index = idx, list(v=2)), "character"))
    expect_length(res1v, 2L)
    expect_true(is(res2v <- bowtie(sequences = list(r1, r2), index = idx, type = "paired", list(v=2)), "character"))
    expect_length(res2v, 0L)
    
    
    ## SpliceMap alignments
    fout1 <- tempfile(fileext = ".sam")
    fout2 <- tempfile(fileext = ".sam")
    fout3 <- tempfile(fileext = ".sam")
    lst <- list(genome_dir = system.file("extdata", "refs", package="Rbowtie"),
                reads_list1 = system.file("extdata", "reads", "reads1.fastq", package = "Rbowtie"),
                read_format = "FASTQ",
                quality_format = "phred-64",
                bowtie_base_dir = file.path(outdir, "index"),
                temp_path = tempdir(),
                num_threads = 2L,
                outfile = fout1,
                selectSingleHit = FALSE)

    expect_error(SpliceMap("error"))
    expect_error(SpliceMap(list()))
    expect_error(SpliceMap(list(outfile = "")))
    expect_error(SpliceMap(list(outfile = "", temp_path = "error")))
    
    expect_is(res1 <- SpliceMap(lst), "character")
    expect_error(SpliceMap(lst), regexp = "output file .+ already exists")
    expect_identical(res1, fout1)

    lst2 <- lst
    lst2$outfile <- fout2
    lst2$num_chromosome_together <- 2L
    
    expect_is(res2 <- SpliceMap(lst2), "character")
    expect_identical(res2, fout2)
    
    rl1 <- readLines(fout1)
    rl2 <- readLines(fout2)
    expect_length(rl1, 9L)
    expect_identical(rl1, rl2)
    
    lst3 <- lst
    lst3$outfile <- fout3
    lst3$reads_list2 <- system.file("extdata", "reads", "reads2.fastq", package = "Rbowtie")
    expect_is(res3 <- SpliceMap(lst3), "character")
    expect_identical(res3, fout3)
    
    rl3 <- readLines(fout3)
    expect_length(rl3, 13L)
    
    unlink(c(fout1, fout2, fout3))
})
