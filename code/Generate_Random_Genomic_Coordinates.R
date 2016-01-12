########################################################################
##		To generate random genomic coordinates
##		http://www.niravmalani.org/generate-random-genomic-coordinates/
########################################################################

## load up necessary libraries ##
library(GenomicRanges)

## set UCSC freeze to obtain the genome size/attributes ##
freeze <- "hg19"

## set the number of sites to generate ##
n <- 4000

## keep chrM, unknown, random or haploid chromosomes? ##
keep_abnormals <- FALSE

z <- gzcon(url(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/", freeze, 
    "/database/chromInfo.txt.gz")))
zlines <- try(readLines(z))
close(z)
if (class(zlines) == "try-error") stop("Could not get thru to UCSC server - try later!")
raw.data <- textConnection(zlines)
chromSizes <- read.delim(raw.data, header = FALSE, stringsAsFactors = FALSE)[, 
    1:2]
chromSizes <- with(chromSizes, structure(V2, names = V1))
close(raw.data)
head(chromSizes)
tail(chromSizes)

# order & filter chromosomes #
if (!keep_abnormals) {
    chromSizes <- chromSizes[!grepl("M|un|random|hap", names(chromSizes), ignore.case = T)]
}

chromSizes <- chromSizes[order(suppressWarnings(as.numeric(sub("chr(\\d+)\\_?.*", 
    "\\1", names(chromSizes)))))]
chromSizes

cumulativeSize <- cumsum(as.numeric(chromSizes))
names(cumulativeSize) <- names(chromSizes)
cumulativeSize

genomeLength <- sum(as.numeric(chromSizes))
genomeLength

z <- gzcon(url(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/", freeze, 
    "/database/gap.txt.gz")))
zlines <- try(readLines(z))
close(z)
if (class(zlines) == "try-error") stop("Could not get thru to UCSC server - try later!")
raw.data <- textConnection(zlines)
gaps <- read.delim(raw.data, header = FALSE, stringsAsFactors = FALSE)
close(raw.data)

gaps <- with(gaps, GRanges(seqnames = V2, IRanges(V3, V4)))
gaps <- keepSeqlevels(gaps, names(chromSizes))
gaps

test <- TRUE
sites.good <- GRanges()
while (any(test)) {
    cat(".")
    ## sample n sites randomly from a linearized genome
    sites <- sample(genomeLength, n, replace = TRUE)

    ## detect which part of the genome is the site sampled from
    cuts <- structure(names(cumulativeSize), names = as.character(cut(cumulativeSize, 
        c(0, cumulativeSize))))
    sites.chr <- cut(sites, c(0, cumulativeSize))
    names(sites) <- cuts[as.character(sites.chr)]

    ## get the actual region of the chromosome by subtracting cumulative size
    ## ##
    sites <- cumulativeSize[names(sites)] - sites

    sites <- GRanges(seqnames = names(sites), IRanges(start = sites, width = 1))

    ## check if the sites are in a gap by any chance ##
    test <- overlapsAny(sites, gaps, type = "any", ignore.strand = TRUE)
    n <- n - table(test)["FALSE"]
    sites.good <- c(sites.good, sites[!test, ])
}


## sort the random good sites and add strand info if needed ##
sites.good <- sort(sites.good)
strand(sites.good) <- sample(c("+", "-"), length(sites.good), replace = TRUE)

## result ##
sites.good

#To store sites as a vector
Output_sites.good=as.data.frame(sites.good)
write.table(Output_sites.good, file="Random_Genome_4000_SV.txt")
