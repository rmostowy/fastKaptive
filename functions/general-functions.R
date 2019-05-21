suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
col.red <- brewer.pal(9, "Set1")[1]
col.blue <- brewer.pal(9, "Set1")[2]
col.green <- brewer.pal(9, "Set1")[3]
col.purple <- brewer.pal(9, "Set1")[4]
col.orange <- brewer.pal(9, "Set1")[5]
col.yellow <- brewer.pal(9, "Set1")[6]
col.brown <- brewer.pal(9, "Set1")[7]
col.pink <- brewer.pal(9, "Set1")[8]
col.grey <- brewer.pal(9, "Set1")[9]

make.dir <- function(dir.name) system(sprintf("mkdir -p %s", dir.name))

delete.file <- function(filename) system(sprintf("rm -rf %s", filename))

get.random.dir.name <- function(no.characters = 9){
  folder.vector <- sample(c(letters, LETTERS, 0:9), no.characters, replace = TRUE)
  return(paste(folder.vector, collapse=""))
}

write.bash <- function(content, file.path, overwrite = TRUE){
  contains.right.extension <- grepl(".sh", file.path)
  content.is.a.string <- is.character(content)
  file.already.exists <- file.exists(file.path)
  proceed <- TRUE
  if(file.already.exists & !overwrite) proceed <- FALSE
  if(!contains.right.extension) stop("Needs the right extension! Exiting...\n")
  if(!content.is.a.string) stop("Content needed in the form of a string! Exiting...\n")
  if(!proceed) stop("File already exists! Use 'overwrite = TRUE' option. Exiting...\n")
  header <- "#!/bin/bash\n"
  full.content <- c(header, content)
  write(full.content, file.path, ncolumns = 1)
  system(sprintf("chmod +x %s", file.path))
}

equiv.vec <- function(vec1, vec2){
  both.are.vectors <- (is.vector(vec1) | is.null(vec1)) & (is.vector(vec2) | is.null(vec2))
  both.are.nonempty <- length(vec1)>0 & length(vec2)>0
  if(!both.are.vectors) stop("One of the elements is not a vector!")
  if(!both.are.nonempty) stop("One of the vectors is empty")
  all(vec1 %in% vec2) & all(vec2 %in% vec1)
}

read.prokka.output <- function(prokka.folder.path){
  gff.filename <- sprintf("%s/fast-kaptive.gff", prokka.folder.path)
  ffn.filename <- sprintf("%s/fast-kaptive.ffn", prokka.folder.path)
  cds.table <- NULL
  gff.exists <- file.exists(gff.filename)
  ffn.exists <- file.exists(ffn.filename)
  both.files.exist <- gff.exists & ffn.exists
  if(both.files.exist){
    cds.positions <- readLines(gff.filename)
    this.cds.seq.nt <- read.fasta(ffn.filename)
    start.line <- grep("##sequence-region", cds.positions)
    end.line <- grep("##FASTA", cds.positions)
    cds.positions <- cds.positions[(max(start.line)+1):(end.line-1)]
    cds.table <- sapply(1:length(cds.positions), function(k){
      s <- strsplit(cds.positions[k], "\t")[[1]]
      s.info <- strsplit(s[9],";")[[1]]
      contig.name <- s[1]
      start.pos <- as.numeric(s[4])
      end.pos <- as.numeric(s[5])
      strand <- s[7]
      locus.tag.index <- grep("locus_tag", s.info)
      locus.tag.exists <- as.logical(length(locus.tag.index))
      out <- rep(NA, 6)
      if(locus.tag.exists){
        product.string <- gene.string <- name.string <- ""
        product.index <- grep("product", s.info)
        product.string.exists <- as.logical(length(product.index))
        if(product.string.exists) product.string <- strsplit(s.info[product.index], "=")[[1]][2]
        name.index <- grep("Name", s.info)
        name.string.exists <- as.logical(length(name.index))
        if(name.string.exists) name.string <- strsplit(s.info[name.index], "=")[[1]][2]
        gene.index <- grep("gene", s.info)
        gene.string.exists <- as.logical(length(gene.index))
        is.hypothetical <- as.logical(length(grep("hypothetical", product.string)))
        is.putative <- as.logical(length(grep("putative", product.string, ignore.case = T)))
        if(gene.string.exists) gene.string <- strsplit(s.info[gene.index], "=")[[1]][2]
        locus.tag <- strsplit(s.info[locus.tag.index], "=")[[1]][2]
        if(name.string.exists){
          product.name <- name.string
        } else if(gene.string.exists){
          product.name <- gene.string
        } else {
          product.name <- "HG"
        }
        out <- c(contig.name, start.pos, end.pos, strand, locus.tag, product.name)
      }
      out
    })
    cds.table <- as.data.frame(t(cds.table), stringsAsFactors = F)
    colnames(cds.table) <- c("contig.name", "start.pos", "end.pos", "strand", "locus.tag", "product.name")
    cds.table <- cds.table[!is.na(cds.table$contig.name),]
  } else{
    cat("No input files or wrong Prokka path!\n")
  }
  return(cds.table)
}

my.table <- function(vector, freq = FALSE, no.digits = 4){
  this.table <- table(vector)
  this.table.sorted <- sort(this.table, decreasing = T)
  out <- this.table.sorted
  if(freq) out <- round(out/sum(out), d=no.digits)
  return(out)
}

get.hc.clusters <- function(cluster.output){
  no.of.clusters <- max(cluster.output)
  cluster.members <- lapply(1:no.of.clusters, function(cluster.index){
    names(which(cluster.output==cluster.index))
  })
  cluster.size <- sapply(cluster.members, length)
  cluster.numbers <- lapply(1:no.of.clusters, function(k){
    rep(k, cluster.size[k])
  })  
  cluster.members.df <- data.frame(cluster = unlist(cluster.numbers), isolate = unlist(cluster.members))
  return(cluster.members.df)  
}

s2v <- function(s) strsplit(s, "")[[1]]
v2s <- function(v) paste(v, collapse="")

bootstrap.resample <- function(vec){
  vector.length <- length(vec)
  output <- vec
  if(vector.length>1) output <- sample(vec, vector.length, replace = T)
  return(output)
}

get.med <- function(v) as.numeric(quantile(v, probs = 0.5, na.rm = T))
get.lower <- function(v) as.numeric(quantile(v, probs = 0.025, na.rm = T))
get.upper <- function(v) as.numeric(quantile(v, probs = 0.975, na.rm = T))

lw <- function(x) length(which(x))

equiv <- function(x, y) all(x %in% y) & all(y %in% x)
