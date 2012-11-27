### main wrapper around SpliceMap

# run splicemap (complete workflow)
# inputs: list with arguments
# output: SAM file path (invisibly) or exeption (on any failure)
SpliceMap <- function(cfg) {
    # minimal checks of input argument (mostly done by .write_cfg)
    if(!is.list(cfg))
        stop("'cfg' must be a list with SpliceMap settings")
    if(!('outfile' %in% names(cfg)))
        stop("no 'outfile' element in SpliceMap config list 'cfg'")
    if(file.exists(cfg[['outfile']]))
        stop(sprintf("output file %s already exists",cfg[['outfile']]))
    if(!('temp_path' %in% names(cfg)))
        stop("no 'temp_path' element in SpliceMap config list 'cfg'")
    cfg[['temp_path']] <- sub("(\\\\|/)$","",cfg[['temp_path']]) # remove trailing slash for file.exists() on windows
    if(!file.exists(cfg[['temp_path']]))
        stop(sprintf("'temp_path' %s does not exist",cfg[['temp_path']]))
    
    # prepare subfolder in cfg[['temp_path']] for all temporary files
    outdir <- cfg[['temp_path']]
    cfg[['temp_path']] <- tempfile(pattern="SpliceMapTemp_", tmpdir=cfg[['temp_path']])
    if(!dir.create(cfg[['temp_path']]))
        stop(sprintf("could not create directory %s",cfg[['temp_path']]))
    on.exit({
        if(file.exists(file.path(cfg[['temp_path']], "junction2.sam")))
            file.rename(file.path(cfg[['temp_path']], "junction2.sam"), cfg[['outfile']])
        unlink(cfg[['temp_path']], recursive=TRUE, force=TRUE)
    })
    
    # perform the steps normally performed by runSpliceMap using R functionality that is compatible with windows
    # 1. prepare files (config file, temp and output directories)
    cfgFname <- file.path(cfg[['temp_path']], "run.cfg")
    cfg <- .write_cfg(cfg, cfgFname)

    # 2. generate SpliceMap sequence, identifier, quality and 25mer files
    message("[SpliceMap] splitting reads into 25mers...", appendLF=FALSE)
    callstr <- paste(shQuote(cfgFname), "generate25mers")
    ret <- system2(file.path(system.file(package="Rbowtie"), "runSpliceMap_QuasR"), callstr, stdout=TRUE, stderr=FALSE)
    if(length(ret)>0 && grepl("^Total 25-mer extraction section time", ret[length(ret)]))
        message("done")
    else
        stop("failed while generating 25mers")

    # 3. map 25-mers using bowtie
    message("[SpliceMap] aligning 25mers...", appendLF=FALSE)
    callstr <- paste(ifelse(!is.null(cfg$try_hard) && cfg$try_hard=="yes", "-y", ""),
                     sprintf("-S -k %d -m %d -v 2 -r -p %d --best --strata",
                             cfg[['max_multi_hit']], cfg[['max_multi_hit']], cfg[['num_threads']]),
                     shQuote(cfg[['bowtie_base_dir']]),
                     shQuote(file.path(cfg[['temp_path']], "25mers.map")),
                     shQuote(file.path(cfg[['temp_path']], "25mers.map_unsorted")))
    ret <- system2(file.path(system.file(package="Rbowtie"),"bowtie"), callstr, stdout=TRUE, stderr=TRUE)
    if(length(ret)>0 && grepl("^(Reported [0-9]+ alignments|No alignments)", ret[length(ret)]))
        message("done")
    else
        stop("failed while aligning 25mers")

    # 4. sort bowtie output
    message("[SpliceMap] sorting 25mer-alignments...", appendLF=FALSE)
    callstr <- paste("-idx",
                     shQuote(file.path(cfg[['temp_path']], "25mers.map_unsorted")),
                     shQuote(file.path(cfg[['temp_path']], "25mers.map.out")))
    ret <- system2(file.path(system.file(package="Rbowtie"), "sortsam"), callstr, stdout=TRUE, stderr=TRUE)
    if(length(ret)>0 && grepl("^Finished", ret[length(ret)]))
        message("done")
    else
        stop("failed while sorting 25mer-alignments")
    unlink(file.path(cfg[['temp_path']], "25mers.map_unsorted"), force=TRUE)

    # 5. index mapped 25-mers
    message("[SpliceMap] indexing 25mer-alignments...", appendLF=FALSE)
    callstr <- paste(shQuote(cfgFname), "index25merAlignments")
    ret <- system2(file.path(system.file(package="Rbowtie"), "runSpliceMap_QuasR"), callstr, stdout=TRUE, stderr=FALSE)
    if(length(ret)>0 && grepl("^Total mapping index creation section execution time", ret[length(ret)]))
        message("done")
    else
        stop("failed while indexing 25mer-alignments")
    unlink(file.path(cfg[['temp_path']], "25mers.map.out"))
    unlink(file.path(cfg[['temp_path']], "25mers.map"))

    # 6. call SpliceMap to create spliced alignments for individual chromosomes
    refL <- scan(file.path(cfg[['temp_path']], "ref_list"),
                 what=list(chrName="", chrFile="", chrDir="", chrS="", chrE=""), sep="\t", quiet=TRUE)
    callstr <- sample(paste(shQuote(cfgFname), refL[['chrName']]))
    if(require(parallel, quietly=TRUE)) {
        message("[SpliceMap] finding spliced alignments for ",
                length(refL[['chrName']]), " chromosomes using ",
                cfg[['num_chromosome_together']], " parallel processes...", appendLF=FALSE)
        cl <- makeCluster(cfg[['num_chromosome_together']])
        ret <- parLapplyLB(cl, callstr, system2, command=file.path(system.file(package="Rbowtie"), "SpliceMap"), stdout=TRUE, stderr=FALSE)
        stopCluster(cl)
    } else {
        message("[SpliceMap] finding spliced alignments for ",
                length(refL[['chrName']]), " chromosomes using serial processes...", appendLF=FALSE)
        ret <- lapply(callstr, system2, command=file.path(system.file(package="Rbowtie"), "SpliceMap"), stdout=TRUE, stderr=FALSE)
    }
    if(length(ret)==length(callstr) && all(sapply(ret, length) > 0) && all(grepl("^Total (.+) execution time.*$",unlist(ret))))
        message("done")
    else
        stop("failed while finding spliced alignments")

    # 6a. identify unmapped reads by comparing input read identifiers with identifiers in alignments --> unmapped.sam
    #     (Remark: semi-mapped read pairs will appear as mapped)
    if("reads_list2" %in% names(cfg)) {
        message("[SpliceMap] extracting unmapped reads (paired-end mode)...", appendLF=FALSE)
        callstr <- paste("P", paste(shQuote(file.path(cfg[['temp_path']],
                                                      c('read_1_1','read_1_1.names','read_1_1.quals',
                                                        'read_1_2','read_1_2.names','read_1_2.quals','unmapped.sam',
                                                        sprintf("%s_%s.sam",refL[['chrFile']],refL[['chrS']])))), collapse=" "),
                         sep=" ")
        ret <- system2(file.path(system.file(package="Rbowtie"), "getSpliceMapUnmapped"), callstr, stdout=TRUE, stderr=FALSE)
    } else {
        message("[SpliceMap] extracting unmapped reads (single read mode)...", appendLF=FALSE)
        callstr <- paste("S", paste(shQuote(file.path(cfg[['temp_path']],
                                                      c('read_1_1','read_1_1.names','read_1_1.quals','unmapped.sam',
                                                        sprintf("%s_%s.sam",refL[['chrFile']],refL[['chrS']])))), collapse=" "),
                         sep=" ")
        ret <- system2(file.path(system.file(package="Rbowtie"), "getSpliceMapUnmapped"), callstr, stdout=TRUE, stderr=FALSE)
    }
    if(length(ret)>0 && grepl("^Finished", ret[length(ret)]))
        message("done")
    else
        stop("failed extracting unmapped reads")

    # 7. combine splice alignment sam files of different chromosomes and create junction.bed/junction.sam files
    #    (Remark: at the same time:
    #              - select one random alignment from multiple ones (if selectSingleHit==TRUE)
    #              - remove semi-mapped pair alignments (only one aligned read in the pair)
    #                and add these sequence pairs to unmapped2.sam. This will not remove alignment pairs with exactly
    #                one alignment per read but with set BAM_FUNMAP flags, e.g. pairs that aligned to different chromosomes.)
    message("[SpliceMap] combining spliced alignments...", appendLF=FALSE)
    callstr <- paste(shQuote(file.path(cfg[['temp_path']], '')),
                     shQuote(file.path(cfg[['temp_path']], 'junction')),
                     as.character(as.integer(cfg[['selectSingleHit']])))
    ret <- system2(file.path(system.file(package="Rbowtie"), "amalgamateSAM"), callstr, stdout=TRUE, stderr=FALSE)
    if(length(ret)>0 && grepl("^SAM File Amalgamation Time", ret[length(ret)]))
        message("done")
    else
        stop("failed combining spliced alignments")

    # 7a. combine junction.sam with unmapped reads (unmapped.sam, unmapped2.sam) in read input order --> junction2.sam
    #     (the on.exit() will rename junction2.sam to cfg[['outfile']])
    message("[SpliceMap] reordering final output...", appendLF=FALSE)
    callstr <- shQuote(file.path(cfg[['temp_path']]))
    ret <- system2(file.path(system.file(package="Rbowtie"), "fuseReorder"), callstr, stdout=TRUE, stderr=FALSE)
    if(length(ret)>0 && grepl("^Finished", ret[length(ret)]))
        message("done")
    else
        stop("failed reordering final output")

    return(invisible(cfg[['outfile']]))
}



### helper functions
# create run.cfg file
.write_cfg <- function(lst, file) {
    # write out SpliceMap config file and perform parameter validity tests normally done by runSpliceMap
    # return list with (potentially completed) config options corresponding to the written config file

    # check if required parameters are given (NULL will select default)
    req <- c("genome_dir","reads_list1","read_format","bowtie_base_dir","temp_path","num_threads","outfile","selectSingleHit")
    if(length(f <- req[!req %in% names(lst)]))
        stop(sprintf("required settings are missing for SpliceMap config file: %s", paste(f,collapse=", ")))

    if(lst[['read_format']] == 'FASTQ') {
        # check quality string
        if(!("quality_format" %in% names(lst)))
            stop("required settings are missing for SpliceMap config file: quality_format")
        if(!(lst[["quality_format"]] %in% c("phred-33","phred-64","solexa")))
            stop("'quality_format' must be one of: phred-33, phred-64, solexa")
    }

    # make sure that only a single (pair of) read file(s) is given (necessary for getSpliceMapUnmapped)
    if(length(lst[["reads_list1"]]) != 1)
        stop(sprintf("'reads_list1' has not exactly one element; all reads need to be contained in a single file"))

    # consistent behavior of file.exists on windows systems. file.exists("dir/") != file.exists("dir").
    for(d in c("genome_dir","bowtie_base_dir","temp_path"))
        lst[[d]] <- sub("(\\\\|/)$","",lst[[d]])

    # check consistency of paired end input
    if("reads_list2" %in% names(lst) && length(lst[["reads_list1"]]) != length(lst[["reads_list2"]]))
        stop("'reads_list1' and 'reads_list2' arguments have not the same length")

    # check existance of bowtie index
    if(!file.exists(lst[["bowtie_base_dir"]])) {
        tmp <- sub(paste(basename(lst[["bowtie_base_dir"]]),"$",sep=""),"",lst[["bowtie_base_dir"]])
        if(length(tmp2 <- list.files(tmp, pattern="*.ebwt")) < 6)
            stop("Invalid bowtie index in 'bowtie_base_dir'.")
    }
  
    # check existance of temp_path
    if(!file.exists(lst[["temp_path"]]))
        stop(sprintf("'temp_path' (%s) does not exists", lst[["temp_path"]]))

    if(!file.exists(tmp <- file.path(lst[["temp_path"]], "debug_logs")))
        if(!dir.create(tmp))
            stop(sprintf("could not create directory %s",tmp))
    on.exit(unlink(tmp, recursive=TRUE, force=TRUE))

    # collapse input file list(s)
    if(any(f <- !file.exists(lst[["reads_list1"]])))
        stop(sprintf("non-existing sequence files in 'reads_list1': %s", paste(lst[["reads_list1"]][f], collapse=", ")))
    lst[["reads_list1"]] <- paste(lst[["reads_list1"]], collapse="\n")
    if("reads_list2" %in% names(lst)) {
        lst[["reads_list2"]] <- paste(lst[["reads_list2"]], collapse="\n")
        if(any(f <- !file.exists(lst[["reads_list2"]])))
            stop(sprintf("non-existing sequence files in 'reads_list2': %s", paste(lst[["reads_list2"]][f], collapse=", ")))
    }

    # add chromosome file list
    if(!file.exists(lst[["genome_dir"]]))
        stop(sprintf("non-existing 'genome_dir': %s",lst[["genome_dir"]]))
       
    if(!(file.info(lst[["genome_dir"]])$isdir)) {
        # 'genome_dir' is a file
        lst[["genome_files"]] <- normalizePath(lst[["genome_dir"]], winslash="/")
    } else {
        # 'genome_dir' is a directory
        chrs <- normalizePath(list.files(lst[["genome_dir"]], pattern="\\.fa$|\\.fna$|\\.fasta$", full.names = TRUE))
        if(length(chrs) == 0)
            stop(sprintf("No chromosome files found in '%s' using pattern '\\.fa$|\\.fna$|\\.fasta$'", lst[["genome_dir"]]), winslash="/")
        lst[["genome_files"]] <- paste(chrs, collapse="\n")
    }

    # check selectSingleHit
    if(!is.logical(lst[["selectSingleHit"]]) || length(lst[["selectSingleHit"]])!=1)
        stop(sprintf("'selectSingleHit' must be TRUE or FALSE"))
    
    # set defaults (replace NULL)
    if(!"max_intron" %in% names(lst)  ||  is.null(lst[["max_intron"]]))
        lst[["max_intron"]] <- 400000
    if(!"min_intron" %in% names(lst)  ||  is.null(lst[["min_intron"]]))
        lst[["min_intron"]] <- 20000
    if(!"max_multi_hit" %in% names(lst)  ||  is.null(lst[["max_multi_hit"]]))
        lst[["max_multi_hit"]] <- 10
    if(!"seed_mismatch" %in% names(lst)  ||  is.null(lst[["seed_mismatch"]]))
        lst[["seed_mismatch"]] <- 1
    if(!"read_mismatch" %in% names(lst)  ||  is.null(lst[["read_mismatch"]]))
        lst[["read_mismatch"]] <- 2
    if(!"num_chromosome_together" %in% names(lst)  ||  is.null(lst[["num_chromosome_together"]]))
        lst[["num_chromosome_together"]] <- 2
    if(!"try_hard" %in% names(lst)  ||  is.null(lst[["try_hard"]]))
        lst[["try_hard"]] <- "yes"

    # set (overrule) expected settings
    lst[["sam_file"]] <- "sam"
    lst[["ud_coverage"]] <- "no"

    # copy and transform to character
    lst2 <- lst
    if(is.numeric(lst2[["min_intron"]]))
        lst2[["min_intron"]] <- sprintf("%d",lst2[["min_intron"]])
    if(is.numeric(lst2[["max_intron"]]))
        lst2[["max_intron"]] <- sprintf("%d",lst2[["max_intron"]])

    # make sure only single values are given
    if(any(f <- sapply(lst2, length) != 1))
        stop(sprintf("only a single value can be given to settings: %s", paste(names(lst2)[f], collapse=", ")))

    # create output string
    i <- grep("reads_list|genome_files",names(lst2))
    out <- c("# SpliceMap config file, created by QuasR",
             "> reads_list1", lst2[["reads_list1"]], "<",
             if("reads_list2" %in% names(lst2)) 
             c("> reads_list2", lst2[["reads_list2"]], "<")
             else
             character(0),
             "> genome_files", lst2[["genome_files"]], "<",
             sprintf("%s = %s", names(lst2)[-i], sapply(lst2[-i], as.character)))
    write.table(out, file, quote=FALSE, col.names=FALSE, row.names=FALSE)

    on.exit()
    return(lst)
}
