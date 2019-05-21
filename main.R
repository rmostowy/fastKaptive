##################################################
############ Written by Rafal Mostowy ############
################# Version: 0.2.2 #################
######## Date: Tue Jun 26 16:23:10 2018 ##########
##################################################
##################################################

suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(optparse)))

cmdArgs <- commandArgs(trailingOnly = FALSE)
thisFile <- function() {
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

############################## GLOBAL PARAMETERS #####################################
## -- 1. cds.coverage.threshold --
## to consider a gene present in the assembly, I assume 
## a minimum coverage of this gene in the assembly; otherwise ignoring it
## -- 2. klocus.overlap.threshold --
## a gene is considered located outside the K-locus 
## if less than 'klocus.overlap.threshold' of its 
## length overlaps with the best-match alignment
## -- 3. cds.overlap.threshold --
## a gene is considered "extra and within locus" if 
## (a) fulfils the above assumption and 
## (b) has less than 'cds.overlap.threshold' overlap 
## -- 4. locus.length.comparison.margin --
## when comparing best hit with bm, allow this error 
## of margin (quivalent to 'start_end_margin' in Kaptive)
## -- 5. hit.overlap.proximity --
## when comparing locus.regions, merge them if they're close;
## closeness is defined as a maximum distance of 'hit.overlap.proximity'
## (quivalent to 'gap_fill_size' in Kaptive)
## -- 6. contamination.jaccard.threshold --
## when comparing locus.regions, merge them if they're close;
## closeness is defined as a maximum distance of 'hit.overlap.proximity'
## (quivalent to 'gap_fill_size' in Kaptive)
cds.coverage.threshold <- 95
klocus.overlap.threshold <- 0.5
cds.overlap.threshold <- 0.7
locus.length.comparison.margin <- 100
hit.overlap.proximity <- 120
contamination.jaccard.threshold <- 0.7
#######################################################################################
# core.genes.set1 <- c("galF", "gnd", "ugd", "hisI")
# core.genes.set2 <- core.genes.set1[1:3]

if(length(cmdArgs)==2){
  this.dir <- getwd()
} else{
  this.dir <- thisFile()
  filename <- strsplit(this.dir,"/")[[1]]
  filename <- filename[length(filename)]
  this.dir <- gsub(sprintf("/%s", filename), "", this.dir)
}

source(sprintf("%s/functions/general-functions.R", this.dir))
source(sprintf("%s/functions/specialist-functions.R", this.dir))

# setwd("/Users/rmostowy/Dropbox/Projects/CpsEvolution/CapsuleEntero/data/2_getGenera/genus15_Klebsiella/zz_test")

option_list <- list(
  make_option(c("-l", "--locus.name"), default="cps",
              help="Which locus to search: cps [default] or kps"),
  make_option(c("-e", "--extended.db"), action="store_true", default=FALSE,
              help="Use an extended locus search database  [default %default]"),
  make_option(c("-D", "--distant.search"), action="store_true", default=FALSE,
              help="Search in distant bacteria using blastn instead of megablast [default %default]"),
  make_option(c("-p", "--path.file"), default="paths.input.txt",
              help="File with path names to all assemblies [default %default]"),
  make_option(c("-o", "--output.folder"), default="extract.output",
              help="Name of the output folder [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

locus.name <- opt$locus.name
extended.db <- opt$extended.db
search.in.distant.species <- opt$distant.search
genome.paths.filename <- opt$path.file
output.folder.name <- opt$output.folder

# error handling
locus.test <- locus.name %in% c("cps", "kps")
if(!locus.test) stop("Wrong locus name, needs to be cps or kps!")
paths.test <- file.exists(genome.paths.filename) & file.size(genome.paths.filename)>100
if(!paths.test) stop("Assembly paths need to be provided!\nSee '-h' option for further help.")

# locus.name="kps"
# extended.db=TRUE
# search.in.distant.species=FALSE
# genome.paths.filename="assembly-paths.txt"

if(locus.name=="cps"){
  database.folder.path <- sprintf("%s/reference-db/cps/0_combined/database", this.dir)
  if(extended.db) database.folder.path <- gsub("0_combined", "00_extended", database.folder.path)
  core.genes.set1 <- c("galF", "gnd", "ugd", "hisI")
  core.genes.set2 <- core.genes.set1[1:3]
  if(search.in.distant.species) cds.coverage.threshold <- 90
} else if(locus.name=="kps"){
  database.folder.path <- sprintf("%s/reference-db/kps/0_combined/database", this.dir)
  if(extended.db) database.folder.path <- gsub("0_combined", "00_extended", database.folder.path)
  cds.coverage.threshold <- 70
  core.genes.set1 <- c("FEDUCSTM")
  core.genes.set1 <- strsplit(core.genes.set1, "")[[1]]
  core.genes.set1 <- sprintf("kps%s", core.genes.set1)
  core.genes.set2 <- c("DMTECS")
  core.genes.set2 <- strsplit(core.genes.set2, "")[[1]]
  core.genes.set2 <- sprintf("kps%s", core.genes.set2)
} else{
  cat("Wrong locus name! Searching for cps locus.\n")
}

reference.database.filename <- sprintf("%s/references.fa", database.folder.path)
cds.database.filename <- sprintf("%s/genes.fa", database.folder.path)
# jaccard.dist.filename <- sprintf("%s/ref-jaccard.txt", database.folder.path)

# load reference database files
get.names.cmd <- sprintf("cat %s | grep '>'", reference.database.filename)
reference.names <- system(get.names.cmd, intern = T)
reference.names <- substr(reference.names, 2, nchar(reference.names))
cds.seq <- read.fasta(cds.database.filename)
cds.seq.n <- names(cds.seq)
cds.seq.n.ref <- sapply(1:length(cds.seq.n), function(k) strsplit(cds.seq.n[k], "__")[[1]][1])
cds.seq.n.gn <- sapply(1:length(cds.seq.n), function(k) strsplit(cds.seq.n[k], "__")[[1]][2])
# ref.jaccard <- read.table(jaccard.dist.filename, header = T, stringsAsFactors = F)

assembly.path.file.exists <- file.exists(genome.paths.filename)
assembly.path.file.empty <- file.size(genome.paths.filename)==0
if(!assembly.path.file.exists | assembly.path.file.empty){
  stop("Assembly path file does't exist or is empty! Exiting...\n")
} else{
  all.assembly.filesnames <- readLines(genome.paths.filename)
  assembly.filesnames.exist <- file.exists(all.assembly.filesnames)
  assembly.filesnames.nonempty <- file.size(all.assembly.filesnames)>100
  all.assembly.filesnames.exist <- all(assembly.filesnames.exist & assembly.filesnames.nonempty)
  if(!any(assembly.filesnames.exist))
    stop("None of the assembly files exists! Exiting...\n")
  if(!all.assembly.filesnames.exist){
    filenames.which.do.not.exist <- all.assembly.filesnames[!assembly.filesnames.exist]
    cat("Warning! The following assembly filenames do not exist or are empty:", filenames.which.do.not.exist, sep="\n")
  }
}

# all.assembly.filesnames <- Sys.glob("assemblies/GC*")
no.assemblies <- length(all.assembly.filesnames)
system(sprintf("mkdir -p %s", output.folder.name))
temp.dir <- tempdir()
blast.tmp.dir <- sprintf("%s/%s", temp.dir, get.random.dir.name(12))

summary.table <- NULL
cds.table <- NULL
contig.hits.table <- NULL
pb <- txtProgressBar(min = 0, max = no.assemblies, style = 3, initial = 0)
for(assembly.index in 1:no.assemblies){
  # write(assembly.index, "latest.txt", ncolumns = 1)
  # cat(assembly.index,"\n")
  # risk.of.contamination.val <- FALSE
  original.assembly.filename <- all.assembly.filesnames[assembly.index]
  is.compressed <- grepl("^.*(.gz|.bz2|.tar|.zip|.tgz|.gzip|.7z)[[:space:]]*$", original.assembly.filename)
  if(is.compressed){
    uncompress.cmd <- sprintf("gzip -dfk %s", original.assembly.filename)
    system(uncompress.cmd)
    assembly.filename <- strsplit(original.assembly.filename,"\\.gz$")[[1]]
  } else{
    assembly.filename <- original.assembly.filename
  }
  assembly.filename.raw <- strsplit(assembly.filename, "/")[[1]]
  # assembly.filename.raw <- assembly.filename.raw[-length(assembly.filename.raw)]
  assembly.filename.raw <- assembly.filename.raw[length(assembly.filename.raw)]
  assembly.seq <- read.fasta(assembly.filename)
  assembly.seq.n <- names(assembly.seq)
  system(sprintf("mkdir -p %s", blast.tmp.dir))
  
  # blast references against the assembly
  blast.output.file <- sprintf("%s/blast.out", blast.tmp.dir)
  my.blast(assembly.filename, reference.database.filename, blast.output.file, 1e-50, megablast = !search.in.distant.species)
  reference.blast.output <- read.blast.output(blast.output.file)
  
  output.table <- NULL
  locus.regions.hits <- NULL
  summary.df <- data.frame(assembly = assembly.filename.raw,
                           best.match = NA,
                           best.match.cov = NA,
                           no.ref.genes.missing = NA,
                           no.extra.genes.within = NA,
                           no.extra.genes.outside = NA,
                           no.regions.with.bm = NA,
                           regions.with.core.set1 = NA,
                           regions.with.core.set2 = NA,
                           regions.with.all.core = NA,
                           stringsAsFactors = F)
  
  if(!is.null(reference.blast.output)){
    
    # calculate coverage of all references in assembly and pick the best match
    ref.coverage.table <- calculate.reference.coverage(reference.names, reference.blast.output)
    bm.name <- ref.coverage.table$reference.name[1]
    bm.cov <- ref.coverage.table$coverage[1]
    bm.data <- reference.blast.output[which(reference.blast.output$seqid == bm.name),]
    contigs.of.interest <- unique(bm.data$geneid)
    bm.length <- bm.data$seqlen[1]
    
    # remove likely transposons, ie identical hits elsewhere in the assembly
    bm.data <- remove.likely.transposons(bm.data)
    # define locus regions
    locus.regions <- bm.data.to.locus.region(bm.data)
    
    # if one alignment framgent overlaps perfectly with the reference, it's a perfect match
    locus.regions <- perfect.match.test(locus.regions)
    # if start and end are in a single contig and overlap with reference, merge them
    # locus.regions.pre <- locus.regions ###### <- to remove
    locus.regions <- start.end.match.test(bm.data, locus.regions)
    locus.regions <- fill.gaps.novel(bm.data, locus.regions)
    # locus.regions.post <- locus.regions ###### <- to remove
    # if(nrow(locus.regions.pre)!=nrow(locus.regions.post)) cat(assembly.filename.raw,"\n") ###### <- to remove
    # merge all fragments which 
    locus.regions <- merge.overlaps(locus.regions)
    
    # prepare ref cds sequences and contigs of interest
    bm.cds <- cds.seq[which(cds.seq.n.ref == bm.name)]
    bm.cds.n <- names(bm.cds)
    contigs.of.interest.seq.filename <- sprintf("%s/coi.fa", blast.tmp.dir)
    contigs.of.interest.seq <- assembly.seq[which(names(assembly.seq) %in% contigs.of.interest)]
    write.fasta(bm.cds, names=names(bm.cds), sprintf("%s/best.match.cds.fa", blast.tmp.dir))
    write.fasta(contigs.of.interest.seq, names=names(contigs.of.interest.seq), contigs.of.interest.seq.filename)
    
    # blast all CDS against contigs of interest
    gene.blast.output.file <- sprintf("%s/blast.cds.out", blast.tmp.dir)
    my.blast(contigs.of.interest.seq.filename, cds.database.filename, gene.blast.output.file, 1e-50)
    gene.blast.output <- read.blast.output(gene.blast.output.file)
    gene.blast.output <- gene.blast.output[gene.blast.output$qcovs>=cds.coverage.threshold,]

    summary.df$best.match <- bm.name
    summary.df$best.match.cov <- bm.cov
    summary.df$no.ref.genes.missing <- length(bm.cds.n)
    summary.df$no.extra.genes.within <- summary.df$no.extra.genes.outside <- summary.df$no.regions.with.bm <- 0
    
    gene.blast.output.is.ok <- FALSE
    gene.blast.output.is.null <- is.null(gene.blast.output)
    if(!gene.blast.output.is.null)
      if(nrow(gene.blast.output)>0) gene.blast.output.is.ok <- TRUE
    
    if(gene.blast.output.is.ok){
      # extract hits which match the BM
      hits.to.bm <- get.bm.hits(gene.blast.output, bm.cds.n, locus.regions)
      hits.to.bm <- remove.duplicates.from.hits(hits.to.bm)
      hits.to.bm <- hits.to.bm[hits.to.bm$klocus.overlap>=klocus.overlap.threshold,]
      locus.regions <- locus.regions[locus.regions$geneid %in% hits.to.bm$geneid,]
      if(nrow(locus.regions)>0){
        # regions where no bm cds were found with blastn should not be locus regions
        locus.regions <- locus.regions[locus.regions$geneid %in% hits.to.bm$geneid,]
        # extract all other hits 
        hits.to.others <- get.other.hits(gene.blast.output, bm.cds.n, locus.regions, hits.to.bm)
        
        hits.outside.klocus <- hits.to.others[hits.to.others$klocus.overlap<klocus.overlap.threshold,]
        hits.outside.klocus <- remove.duplicates.from.hits(hits.outside.klocus)
        
        hits.inside.extra <- hits.to.others[(hits.to.others$klocus.overlap>=klocus.overlap.threshold) & (hits.to.others$bm.cds.overlap<cds.overlap.threshold),]
        hits.inside.extra <- remove.duplicates.from.hits(hits.inside.extra)
        
        hits.to.bm.duplicated <- hits.to.bm[duplicated(hits.to.bm$seqid),]
        hits.to.bm.duplicated <- hits.to.bm.duplicated[hits.to.bm.duplicated$klocus.overlap>=klocus.overlap.threshold,]
        if(nrow(hits.to.bm.duplicated)>0) hits.inside.extra <- rbind(hits.to.bm.duplicated, hits.inside.extra)
        hits.to.bm <- hits.to.bm[!duplicated(hits.to.bm$seqid),]
        # hits.to.bm <- hits.to.bm[hits.to.bm$klocus.overlap>=klocus.overlap.threshold,]
        
        ref.output.table <- get.ref.output(hits.to.bm, bm.cds.n)
        extra.output.table <- get.extra.output(hits.inside.extra)
        output.table <- rbind(ref.output.table, extra.output.table)
        output.table <- cbind(rep(assembly.filename.raw, nrow(output.table)), output.table)
        colnames(output.table)[1] <- "assembly"
        
        locus.regions.output <- sprintf("%s/locus_%s.fa", output.folder.name, assembly.filename.raw)
        locus.regions.names <- save.locus.region.sequences(locus.regions, assembly.seq, locus.regions.output)
        all.hits.locus <- rbind(hits.to.bm, hits.to.others[hits.to.others$klocus.overlap>=klocus.overlap.threshold,])
        all.hits.locus$gene.name <- sapply(1:nrow(all.hits.locus), function(k) strsplit(all.hits.locus$seqid[k],"__")[[1]][2])
        all.hits.locus$gene.name.raw <- sapply(1:nrow(all.hits.locus), function(k) strsplit(all.hits.locus$gene.name[k],"_")[[1]][1])
        all.hits.locus$overlapping.region <- sapply(1:nrow(all.hits.locus), function(hit.index){
          df <- all.hits.locus[hit.index,]
          region.overlap <- sapply(1:nrow(locus.regions), function(k) overlap(df$sstart, df$send, locus.regions$sstart[k], locus.regions$send[k])>0)
          contig.overlap <- df$geneid == locus.regions$geneid
          which(region.overlap & contig.overlap)[1]
        })
        locus.regions.all.gn <- lapply(1:nrow(locus.regions), function(region.index) unique(all.hits.locus$gene.name.raw[which(all.hits.locus$overlapping.region==region.index)]))
        locus.regions.hits <- get.locus.region.hits(all.hits.locus, locus.regions, locus.regions.names, assembly.filename.raw)

        regions.contain.set1 <- sapply(1:length(locus.regions.all.gn), function(k) all(core.genes.set1 %in% locus.regions.all.gn[[k]]))
        regions.contain.set2 <- sapply(1:length(locus.regions.all.gn), function(k) all(core.genes.set2 %in% locus.regions.all.gn[[k]]))
        regions.names.contain.set1 <- paste(locus.regions.names[regions.contain.set1], collapse=";")
        regions.names.contain.set2 <- paste(locus.regions.names[regions.contain.set2], collapse=";")
        regions.names.contain.all.core <- ""
        overlapping.regions.of.all.core.genes <- all.hits.locus$overlapping.region[all.hits.locus$gene.name.raw %in% core.genes.set1]
        if(length(overlapping.regions.of.all.core.genes)>0){
          if(all(overlapping.regions.of.all.core.genes[1]==overlapping.regions.of.all.core.genes)){
            regions.names.contain.all.core <- locus.regions.names[overlapping.regions.of.all.core.genes[1]]
          }  
        }
        summary.df$no.ref.genes.missing <- length(which(is.na(ref.output.table$pident)))
        summary.df$no.extra.genes.within <- no.elements(hits.inside.extra)
        summary.df$no.extra.genes.outside <- no.elements(hits.outside.klocus)
        summary.df$no.regions.with.bm <- length(locus.regions.names)
        summary.df$regions.with.core.set1 <- regions.names.contain.set1
        summary.df$regions.with.core.set2 <- regions.names.contain.set2
        summary.df$regions.with.all.core <- regions.names.contain.all.core
      }
    }
  }
  cds.table <- rbind(cds.table, output.table)
  summary.table <- rbind(summary.table, summary.df)
  contig.hits.table <- rbind(contig.hits.table, locus.regions.hits)
  if(is.compressed){
    delete.cmd <- sprintf("rm -rf %s", assembly.filename)
    system(delete.cmd)
  }
  delete.file(blast.tmp.dir)
  setTxtProgressBar(pb, assembly.index)
}
cat("\n\nDone!\n")

write.table(summary.table, sprintf("%s/summary.txt", output.folder.name), quote = F, row.names = F)
write.csv(summary.table, sprintf("%s/summary.csv", output.folder.name), row.names = F)
write.table(cds.table, sprintf("%s/cds.txt", output.folder.name), quote = F, row.names = F)
write.table(contig.hits.table, sprintf("%s/contig-hits.txt", output.folder.name), quote = F, row.names = F)
