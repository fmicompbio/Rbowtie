## The main wrapper around bowtie

## The main wrapper around bowtie-build
bowtie_build <- function(references, outdir, ..., prefix="index", force=FALSE, strict=TRUE)
{
    if(strict && (!is.character(references) || !all(file.exists(references))))
        stop("Argument 'references' has to be a character vector of filenames ",
             "for building the sequence index.")
    if(strict && (!is.character(outdir) || length(outdir)!=1))
        stop("Argument 'outdir' must be a character scalar giving the output ",
             "directory to store the bowtie indices in.")
    if(strict && (file.exists(outdir) && !force))
        stop("Directory '", outdir, "' exists. Use 'force=TRUE' to overwrite.")
    if(is.null(list(...)[["usage"]]) || !list(...)[["usage"]])
        dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
    args <- sprintf("%s %s %s/%s", .createFlags(list(...)), paste(references, collapse=","),
                    outdir, prefix)
    return(invisible(.bowtieBin("bowtie-build", args)))
}


## The main wrapper around bowtie
bowtie <- function(sequences, index, ..., type=c("single", "paired", "crossbow"), outfile,
                   force=FALSE, strict=TRUE)
{
    type <- match.arg(type)
    args <- list(...)
    args <- args[setdiff(names(args), c("1", "2", "12"))]
    seqIn <- !is.null(args[["c"]]) && args[["c"]]
    seqArg <- ""
    if(strict)
    {
        seqArg <- switch(type,
                         single={
                             if(!is.character(sequences) || (!seqIn && !all(file.exists(sequences))))
                                 stop("Argument 'sequences' has to be a character vector of filenames ",
                                      "to align against the bowtie index or a character of read ",
                                      "sequences if the additional argument c==TRUE.")
                             paste(sequences, collapse=",")
                         },
                         paired={
                         if(!is.list(sequences) || length(sequences)!=2)
                             stop("Argument 'sequences' must be a list of length 2.")
                         tmp <- NULL
                         for(i in 1:2)
                         {
                             if(!is.character(sequences[[i]]) || (!seqIn && !all(file.exists(sequences[[i]]))))
                                 stop("Argument 'sequences[[", i, "]]' has to be a character vector of filenames ",
                                      "to align against the bowtie index or a character of read ",
                                      "sequences if the additional argument c==TRUE.")
                             tmp <- paste(tmp,  "-", i, " ", paste(sequences[[i]], collapse=","), " ", sep="")
                         }
                         tmp
                     },
                     crossbow={
                         if(!is.character(sequences) || (!seqIn && !all(file.exists(sequences))))
                                 stop("Argument 'sequences' has to be a character vector of filenames ",
                                      "to align against the bowtie index or a character of read ",
                                      "sequences if the additional argument c==TRUE.")
                         paste("-12 ", paste(sequences, collapse=","))
           })
    
        if(!is.character(index) || !file.exists(dirname(index)))
            stop("Argument 'index' has to be a character scalar giving the path to the index directory.")
    }
    outfile <- if(!missing(outfile))
    {
        if(strict && (!is.character(outfile) || length(outfile)!=1))
            stop("Argument 'outfile' must be a character scalar giving the output ",
                 "file name to store the bowtie alignments in.")
        if(strict && (file.exists(outfile) && !force))
            stop("File '", outfile, "' exists. Use 'force=TRUE' to overwrite.")
        sprintf(" %s", outfile)
    } else ""
   
    
    args <- sprintf("%s %s %s %s", .createFlags(args), index, seqArg, outfile)
    return(invisible(.bowtieBin("bowtie", args)))
}

## Little helpers that return a description of the intended usage for bowtie and bowtie-build
bowtie_build_usage <- function()
    print(bowtie_build("dummy", "dummy", force=TRUE, usage=TRUE, strict=FALSE))

bowtie_usage <- function()
    print(bowtie("dummy", "dummy", force=TRUE, usage=TRUE, strict=FALSE))


## A helper function to create a scalar of command line arguments from a named list.
## Logical list entries are being interpreted as flags, all other entries are being
## collapsed into the form '<entryName>=<entryValue>'. Vectors of non-logical entry
## values will be collapsed into a single comma-separated scalar.
.createFlags <- function(flagList)
{
    if(!length(flagList))
        return("")
    if(is.null(names(flagList)) || any(names(flagList)==""))
        stop("Unable to create command line arguments from input.")
    logFlags <- sapply(flagList, is.logical)
    flags <- NULL
    if(any(logFlags))
    {
        fnames <- names(flagList)[logFlags][sapply(flagList[logFlags], function(x) x[1])]
        flags <- paste(sapply(fnames, function(x) ifelse(nchar(x)==1, sprintf("-%s", x), sprintf("--%s", x))),
                       collapse=" ")
    }
    fnames <- sapply(names(flagList)[!logFlags], function(x) ifelse(nchar(x)==1, sprintf("-%s", x),
                                                                    sprintf("--%s", x)))
    flags <- paste(flags, paste(fnames, sapply(flagList[!logFlags], paste, collapse=","),
                                collapse=" ", sep=" "), collapse=" ")
    return(gsub("^ *| *$", "", flags))
}


## A helper function to call one of the two bowtie binaries with additional arguments.
.bowtieBin <- function(bin=c("bowtie", "bowtie-build"), args="")
{
    if(is.null(args) || args=="")
        stop("The bowtie binaries need to be called with additional arguments")
    args <- gsub("^ *| *$", "", args)
    bin <- match.arg(bin)
    call <- paste(file.path(system.file(package="Rbowtie"), bin), args)
    #return(call)
    output <- system(call, intern=TRUE)
    return(output)
}

## The direct binary call function
.execute <- function(callstr){
  call <- file.path(system.file(package="Rbowtie"), callstr)
  return(system(shQuote(call), intern=TRUE))
}

