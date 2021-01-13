context("test_helpers")

test_that(".createFlags works", {
    expect_equal(
        .createFlags(
            flagList=list(p=3, q=TRUE, phred33=TRUE, `no-softclip`=TRUE,
                          secondary=FALSE, vectorArg=c("file1", "file2"))),
        "-q --phred33 --no-softclip -p 3 --vectorArg file1,file2"
    )

    expect_error(
        .createFlags(flagList=list(3, x=4))
    )
})

test_that(".bowtieBin works", {
    expect_error(
        .bowtieBin(bin="bowtie-build", args="")
    )
    expect_error(
        .bowtieBin(bin="bowtie", args="")
    )
    expect_equal(
        .bowtieBin(bin="bowtie", args="-1 file1 -2 file2", execute=FALSE),
        paste(shQuote(file.path(system.file(package="Rbowtie"), "bowtie")),
              "-1 file1 -2 file2")
    )
})

test_that("print usage methods work", {
    expect_is(bowtie_build_usage(), "character")
    expect_is(bowtie_usage(), "character")
    expect_is(bowtie_version(), "character")
})

test_that(".write_cfg works", {
    # prepare fake files and folders
    tf1 <- tempfile(fileext = ".txt") # config file
    tf2 <- tempfile(fileext = ".sam") # alignment file
    td1 <- tempfile() # (fake) bowtie index
    td2 <- tempfile() # fake genome dir
    dir.create(td1)
    dir.create(td2)
    for (i in 1:6) {
        writeLines("fake", file.path(td1, paste0("fake_index_", i, ".ebwt")))
        writeLines("fake", file.path(td2, paste0("chr", i, ".fa")))
    }
    
    # parameter list
    lst <- list(genome_dir = td2,
                reads_list1 = system.file("extdata", "reads", "reads1.fastq", package = "Rbowtie"),
                reads_list2 = system.file("extdata", "reads", "reads2.fastq", package = "Rbowtie"),
                read_format = "FASTQ",
                quality_format = "phred-33",
                bowtie_base_dir = td1,
                temp_path = tempdir(),
                num_threads = 2L,
                outfile = tf2,
                selectSingleHit = FALSE)
    
    # error checks
    expect_error(.write_cfg(list(), tf1))
    lst2 <- lst; lst2$quality_format <- NULL
    expect_error(.write_cfg(lst2, tf1))
    lst2$quality_format <- "error"
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$bowtie_base_dir <- "error"
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$temp_path <- "error"
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$reads_list1 <- "error"
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$reads_list2 <- "error"
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$genome_dir <- "error"
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$genome_dir <- td1
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$genome_dir <- file.path(td2, c("chr1.fa", "chr2.fa")) 
    expect_error(expect_warning(.write_cfg(lst2, tf1)))
    lst2 <- lst; lst2$reads_list1 <- c(lst2$reads_list1, lst2$reads_list1)
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$reads_list2 <- c(lst2$reads_list1, lst2$reads_list1)
    expect_error(.write_cfg(lst2, tf1))
    lst2 <- lst; lst2$selectSingleHit <- "error"
    expect_error(.write_cfg(lst2, tf1))
    
    
    # correct results
    expect_is(res <- .write_cfg(lst, tf1), "list")
    resLines <- readLines(tf1)
    expect_length(resLines, 32L)
    expect_identical(resLines[c(3,6)],
                     system.file("extdata", "reads", paste0("reads",1:2,".fastq"),
                                 package = "Rbowtie"))
    expect_identical(resLines[9:14],
                     normalizePath(file.path(td2, paste0("chr", 1:6, ".fa")), winslash = "/"))
    expect_identical(resLines[16], paste("genome_dir","=",td2))
    expect_identical(resLines[19], paste("bowtie_base_dir","=",td1))
    expect_identical(resLines[22], paste("outfile","=",tf2))
    expect_identical(resLines[23], "selectSingleHit = FALSE")
})