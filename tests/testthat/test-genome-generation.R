context("genome_generation")

test_that("malformed input gives error", {
    expect_error(
        bowtie_build(references=TRUE, outdir="./", force=TRUE, execute=FALSE)
    )
    expect_error(
        bowtie_build(references="nonexistent.fa",
                     outdir="./", force=TRUE, execute=FALSE)
    )
    expect_error(
        bowtie_build(
            references=system.file("extdata/refs/chr1.fa", package="Rbowtie"),
            outdir=1, force=TRUE, execute=FALSE)
    )
    expect_error(
        bowtie_build(
            references=system.file("extdata/refs/chr1.fa", package="Rbowtie"),
            outdir=c("dir1", "dir2"), force=TRUE, execute=FALSE)
    )

    tmp <- tempdir()
    dir.create(file.path(tmp, "mytestdir"), recursive=TRUE, showWarnings=FALSE)
    expect_error(
        bowtie_build(
            references=system.file("extdata/refs/chr1.fa", package="Rbowtie"),
            outdir=file.path(tmp, "mytestdir"), force=FALSE, execute=FALSE)
    )
})

test_that("correctly formatted input works", {
    tmp <- tempdir()
    expect_equal(
        bowtie_build(
            references=system.file("extdata/refs/chr1.fa", package="Rbowtie"),
            outdir=file.path(tmp, "bowtieindexdir"), force=TRUE, execute=FALSE),
        paste(shQuote(file.path(system.file(package="Rbowtie"), "bowtie-build")),
              shQuote(system.file("extdata/refs/chr1.fa", package="Rbowtie")),
              shQuote(paste0(file.path(tmp, "bowtieindexdir"), "/index")))
    )
})
