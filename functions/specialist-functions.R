file.test <- function(filename){
  test.result <- FALSE
  file.exists <- file.exists(filename)
  if(file.exists){
    file.nonempty <- file.size(filename)>0
    if(file.nonempty)
      test.result <- TRUE
  }
  return(test.result)
}
data.frame.nonempty <- function(input.df){
  input.ok <- FALSE
  input.nonnull <- !is.null(input.df)
  if(input.nonnull)
    if(nrow(input.df)>0)
      input.ok <- TRUE  
  return(input.ok)
}
no.elements <- function(input.df){
  out <- 0
  if(!is.null(input.df)) out <- nrow(input.df)
  out
}
fasta.file.test <- function(filename){
  test.result <- FALSE
  file.passes.test <- file.test(filename)
  if(file.passes.test){
    no.seqs <- length(system(sprintf("cat %s | grep '>'", filename), intern = T))
    if(no.seqs>0) test.result <- TRUE
  }
  return(test.result)
}
overlap <- function(min1, max1, min2, max2){
  range1 <- range(min1, max1)
  range2 <- range(min2, max2)
  true.min1 <- range1[1]
  true.max1 <- range1[2]
  true.min2 <- range2[1]
  true.max2 <- range2[2]
  max(0, min(true.max1, true.max2) - max(true.min1, true.min2))
}
does.it.overlap <- function(subject.start, subject.end, query.start, query.end) !((query.end < subject.start) | (query.start > subject.end))
does.it.lie.in.close.proximity <- function(subject.start, subject.end, query.start, query.end, distance.threshold){
  out <- FALSE
  there.is.overlap <- does.it.overlap(subject.start, subject.end, query.start, query.end)
  if(there.is.overlap){
    out <- TRUE
  } else{
    subject.first <- subject.start<query.start
    if(subject.first){
      dist <- query.start-subject.end
    } else{
      dist <- subject.start-query.end
    }
    if(dist<=distance.threshold) out <- TRUE    
  }
  out
} 
is.contained <- function(subject.start, subject.end, query.start, query.end){
  containment1 <- (subject.start<=query.start)&(subject.end>=query.end)
  containment2 <- (subject.start>=query.start)&(subject.end<=query.end)
  containment1 | containment2
}
are.genes.same <- function(subject.start, subject.end, query.start, query.end, overlap.threshold){
  genes.are.same <- FALSE
  genes.overlap <- does.it.overlap(subject.start, subject.end, query.start, query.end)
  genes.are.contained <- is.contained(subject.start, subject.end, query.start, query.end)
  if(genes.overlap){
    if(genes.are.contained){
      genes.are.same <- TRUE
    } else{
      subject.length <- subject.end-subject.start+1
      query.length <- query.end-query.start+1
      if(query.start<subject.start) al.length <- query.end-subject.start+1
      if(query.start>subject.start) al.length <- subject.end-query.start+1
      max.overlap <- max(al.length/subject.length, al.length/query.length)
      if(max.overlap>=overlap.threshold) genes.are.same <- TRUE
    }
  }
  return(genes.are.same)  
}
remove.duplicate.rows <- function(blast.table){
  blast.table.dupl <- blast.table[,colnames(blast.table) %in% c("geneid", "sstart", "send")]
  return(blast.table[!duplicated(blast.table.dupl),])
}
get.strand <- function(start,end){
  strand <- "+"
  if(start>end) strand <- "-"
  return(strand)
}
perfect.match.test <- function(locus.regions){
  locus.regions.sorted <- locus.regions[order(-locus.regions$length),]
  longest.hit.length <- locus.regions.sorted$length[1]
  is.bm.match <- abs(longest.hit.length-bm.length)<=locus.length.comparison.margin
  if(is.bm.match) locus.regions <- locus.regions.sorted[1,]
  return(locus.regions)
}

fill.gaps.novel <- function(bm.data, locus.regions){
  bm.data.sorted <- bm.data[order(-bm.data$length),]
  bm.data.sorted$strand <- rep("+", nrow(bm.data.sorted))
  bm.data.sorted$strand[bm.data.sorted$sstart>bm.data.sorted$send] <- "-"
  test.result <- FALSE
  first.piece.index <- which(min(bm.data.sorted$seqstart)==bm.data.sorted$seqstart)[1]
  last.piece.index <- which(max(bm.data.sorted$seqend)==bm.data.sorted$seqend)[1]
  if(first.piece.index != last.piece.index){
    # first.piece.range <- c(bm.data.sorted$seqstart[first.piece.index],bm.data.sorted$seqend[first.piece.index])
    first.piece.range <- c(bm.data.sorted$sstart[first.piece.index],bm.data.sorted$send[first.piece.index])
    first.piece.strand <- bm.data.sorted$strand[first.piece.index] #get.strand(bm.data.sorted$sstart[first.piece.index],bm.data.sorted$send[first.piece.index])
    # last.piece.range <- c(bm.data.sorted$seqstart[last.piece.index],bm.data.sorted$seqend[last.piece.index])
    last.piece.range <- c(bm.data.sorted$sstart[last.piece.index],bm.data.sorted$send[last.piece.index])
    last.piece.strand <- bm.data.sorted$strand[last.piece.index]  #get.strand(bm.data.sorted$sstart[last.piece.index],bm.data.sorted$send[last.piece.index])
    found.in.same.contig <- bm.data.sorted$geneid[first.piece.index]==bm.data.sorted$geneid[last.piece.index]
    same.strand <- first.piece.strand==last.piece.strand
    span.range <- range(first.piece.range, last.piece.range)
    span.length <- span.range[2]-span.range[1]+1
    if(found.in.same.contig & same.strand & span.length < 50000) test.result <- TRUE
  }
  if(test.result){
    contig.range <- range(c(bm.data.sorted$sstart[first.piece.index],bm.data.sorted$sstart[first.piece.index],
                            bm.data.sorted$send[last.piece.index],bm.data.sorted$send[last.piece.index]))
    contig.pos <- contig.range[1]:contig.range[2]
    contig.name <- bm.data.sorted$geneid[first.piece.index]
    which.strand <- bm.data.sorted$strand[first.piece.index]
    locus.regions.rows.to.merge <- sapply(1:nrow(locus.regions), function(locus.regions.index){
      pos.overlap <- any(locus.regions$sstart[locus.regions.index]:locus.regions$send[locus.regions.index] %in% contig.pos)
      contig.overlap <- bm.data.sorted$geneid[locus.regions.index]==contig.name
      strand.agreement <- bm.data.sorted$strand[locus.regions.index]==which.strand
      pos.overlap & contig.overlap & strand.agreement
    })
    if(which.strand=="+"){
      sstart <- min(locus.regions$sstart[locus.regions.rows.to.merge])
      send <- max(locus.regions$send[locus.regions.rows.to.merge])
      length <- send-sstart+1
    } else{
      sstart <- max(locus.regions$sstart[locus.regions.rows.to.merge])
      send <- min(locus.regions$send[locus.regions.rows.to.merge])
      length <- sstart-send+1
    }
    locus.regions <- rbind(data.frame(geneid=contig.name, sstart, send, length, strand=which.strand, stringsAsFactors = F), locus.regions[!locus.regions.rows.to.merge,])
  }
  return(locus.regions)
}

start.end.match.test <- function(bm.data, locus.regions){
  bm.data.sorted <- bm.data[order(-bm.data$length),]
  bm.data.sorted$strand <- rep("+", nrow(bm.data.sorted))
  bm.data.sorted$strand[bm.data.sorted$sstart>bm.data.sorted$send] <- "-"
  test.result <- FALSE
  first.piece.index <- which(min(bm.data.sorted$seqstart)==bm.data.sorted$seqstart)[1]
  last.piece.index <- which(max(bm.data.sorted$seqend)==bm.data.sorted$seqend)[1]
  if(first.piece.index != last.piece.index){
    # first.piece.range <- c(bm.data.sorted$seqstart[first.piece.index],bm.data.sorted$seqend[first.piece.index])
    first.piece.range <- c(bm.data.sorted$sstart[first.piece.index],bm.data.sorted$send[first.piece.index])
    first.piece.strand <- bm.data.sorted$strand[first.piece.index] #get.strand(bm.data.sorted$sstart[first.piece.index],bm.data.sorted$send[first.piece.index])
    # last.piece.range <- c(bm.data.sorted$seqstart[last.piece.index],bm.data.sorted$seqend[last.piece.index])
    last.piece.range <- c(bm.data.sorted$sstart[last.piece.index],bm.data.sorted$send[last.piece.index])
    last.piece.strand <- bm.data.sorted$strand[last.piece.index]  #get.strand(bm.data.sorted$sstart[last.piece.index],bm.data.sorted$send[last.piece.index])
    found.in.same.contig <- bm.data.sorted$geneid[first.piece.index]==bm.data.sorted$geneid[last.piece.index]
    same.strand <- first.piece.strand==last.piece.strand
    if(found.in.same.contig & same.strand){
      span.range <- range(first.piece.range, last.piece.range)
      span.length <- span.range[2]-span.range[1]+1
      if(abs(bm.length-span.length)<locus.length.comparison.margin) test.result <- TRUE
    }  
  }
  if(test.result){
    contig.range <- range(c(bm.data.sorted$sstart[first.piece.index],bm.data.sorted$sstart[first.piece.index],
                                  bm.data.sorted$send[last.piece.index],bm.data.sorted$send[last.piece.index]))
    contig.pos <- contig.range[1]:contig.range[2]
    contig.name <- bm.data.sorted$geneid[first.piece.index]
    which.strand <- bm.data.sorted$strand[first.piece.index]
    
    locus.regions.rows.to.merge <- sapply(1:nrow(locus.regions), function(locus.regions.index){
      pos.overlap <- any(locus.regions$sstart[locus.regions.index]:locus.regions$send[locus.regions.index] %in% contig.pos)
      contig.overlap <- bm.data.sorted$geneid[locus.regions.index]==contig.name
      strand.agreement <- bm.data.sorted$strand[locus.regions.index]==which.strand
      pos.overlap & contig.overlap & strand.agreement
    })
    if(which.strand=="+"){
      sstart <- min(locus.regions$sstart[locus.regions.rows.to.merge])
      send <- max(locus.regions$send[locus.regions.rows.to.merge])
      length <- send-sstart+1
    } else{
      sstart <- max(locus.regions$sstart[locus.regions.rows.to.merge])
      send <- min(locus.regions$send[locus.regions.rows.to.merge])
      length <- sstart-send+1
    }
    locus.regions <- rbind(data.frame(geneid=contig.name, sstart, send, length, strand=which.strand, stringsAsFactors = F), locus.regions[!locus.regions.rows.to.merge,])
  }
  return(locus.regions)
}
merge.overlaps <- function(locus.regions){
  if(nrow(locus.regions)>1){
    overlapping.matrix <- matrix(0, nrow = nrow(locus.regions), ncol = nrow(locus.regions))  
    for(this.piece.index in 1:nrow(locus.regions)){
      for(other.piece.index in 1:nrow(locus.regions)){
        this.piece.strand <- locus.regions$strand[this.piece.index]
        this.piece.range <- sort(c(locus.regions$sstart[this.piece.index],locus.regions$send[this.piece.index]))
        this.piece.contig <- locus.regions$geneid[this.piece.index]
        other.piece.strand <- locus.regions$strand[other.piece.index]
        other.piece.range <- sort(c(locus.regions$sstart[other.piece.index],locus.regions$send[other.piece.index]))
        other.piece.contig <- locus.regions$geneid[other.piece.index]
        # positions.overlap <- does.it.overlap(this.piece.range[1], this.piece.range[2], other.piece.range[1], other.piece.range[2])
        strands.are.same <- this.piece.strand==other.piece.strand
        contigs.are.same <- this.piece.contig==other.piece.contig
        positions.close.or.overlap <- does.it.lie.in.close.proximity(this.piece.range[1], this.piece.range[2], other.piece.range[1], other.piece.range[2], hit.overlap.proximity)
        if(positions.close.or.overlap & strands.are.same & contigs.are.same) overlapping.matrix[this.piece.index, other.piece.index] <- 1  
      }
    }
    g <- graph_from_adjacency_matrix(overlapping.matrix)
    g.clusters <- clusters(g)
    
    no.pieces <- g.clusters$no
    locus.regions.new <- NULL
    for(piece.index in 1:no.pieces){
      row.indices <- which(g.clusters$membership==piece.index)
      geneid <- locus.regions$geneid[row.indices][1]
      strand <- locus.regions$strand[row.indices][1]
      if(strand=="+"){
        sstart <- min(locus.regions$sstart[row.indices])
        send <- max(locus.regions$send[row.indices])
        length <- send-sstart+1
      } else{
        sstart <- max(locus.regions$sstart[row.indices])
        send <- min(locus.regions$send[row.indices])
        length <- sstart-send+1
      }
      locus.regions.new <- rbind(locus.regions.new, data.frame(geneid,sstart,send,length,strand, stringsAsFactors = F))
    }
    locus.regions <- locus.regions.new
  }
  return(locus.regions)
}
remove.duplicates.from.hits <- function(hit.table){
  if(!is.null(hit.table)){
    duplicate.names <- NULL
    if(nrow(hit.table)>1){
      duplicate.names <- lapply(1:(nrow(hit.table)-1), function(subject.row.index){
        # cat(subject.row.index,"\n")
        # hit.table[subject.row.index,]
        subject.contig <- hit.table$geneid[subject.row.index]
        subject.range <- c(hit.table$sstart[subject.row.index], hit.table$send[subject.row.index])
        subject.start <- min(subject.range)
        subject.end <- max(subject.range)
        remaining.hits <- hit.table[(subject.row.index+1):nrow(hit.table),]
        no.remaining.rows <- nrow(remaining.hits)
        remaining.hits.are.duplicates <- sapply(1:no.remaining.rows, function(query.row.index){
          query.contig <- remaining.hits$geneid[query.row.index]
          query.range <- c(remaining.hits$sstart[query.row.index], remaining.hits$send[query.row.index])
          query.start <- min(query.range)
          query.end <- max(query.range)
          query.is.duplicate <- FALSE
          if(subject.contig==query.contig) query.is.duplicate <- are.genes.same(subject.start, subject.end, query.start, query.end, cds.overlap.threshold)
          query.is.duplicate  
        })
        remaining.hits$seqid[remaining.hits.are.duplicates]
      })
      duplicate.names <- unique(unlist(duplicate.names))
      hit.table <- hit.table[!hit.table$seqid %in% duplicate.names,]  
    } 
  }
  return(hit.table)
}
get.locus.region.hits <- function(all.hits.locus, locus.regions, locus.regions.names, assembly.filename.raw){
  locus.region.hits <- NULL
  for(region.index in 1:nrow(locus.regions)){
    this.region.table <- all.hits.locus[which(all.hits.locus$overlapping.region==region.index),]
    out <- NULL
    if(nrow(this.region.table)>0){
      this.region.strand <- locus.regions$strand[region.index]
      this.region.table.unique <- remove.duplicates.from.hits(this.region.table)
      this.region.table.unique <- this.region.table.unique[order(this.region.table.unique$sstart),]
      if(this.region.strand=="-") this.region.table.unique <- this.region.table.unique[nrow(this.region.table.unique):1,]
      out <- data.frame(locus.region.name=rep(locus.regions.names[region.index], nrow(this.region.table.unique)),
                        gene.name=edit.gene.names(this.region.table.unique$gene.name.raw), 
                        gene.start=this.region.table.unique$sstart, 
                        gene.end=this.region.table.unique$send, 
                        gene.query=this.region.table.unique$seqid,
                        gene.qcovs=this.region.table.unique$qcovs,
                        gene.eval=this.region.table.unique$evalue,
                        gene.pident=this.region.table.unique$pident,
                        stringsAsFactors = F)
    }
    locus.region.hits <- rbind(locus.region.hits, out)
  }
  locus.region.hits <- data.frame(assembly = rep(assembly.filename.raw, nrow(locus.region.hits)), locus.region.hits, stringsAsFactors = F)
  return(locus.region.hits)
}

my.blast <- function(subject.filename, query.filename, output.filename, e.value, megablast = FALSE, max.target = 40, sseq = FALSE, qseq = FALSE){
  options <- "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
  if(sseq) options <- sprintf("%s sseq", options)
  if(qseq) options <- sprintf("%s qseq", options)
  subject.file.test <- fasta.file.test(subject.filename)
  query.file.test <- fasta.file.test(query.filename)
  task <- "blastn"
  if(megablast) task <- "megablast"
  if(subject.file.test & query.file.test){
    blast.cmd <- sprintf("blastn -subject %s -query %s -evalue %.0e -task %s -max_target_seqs %d -out %s -outfmt '%s'", subject.filename, query.filename, e.value, task, max.target, output.filename, options)
    system(blast.cmd)
  } else{
    cat("Query and subject files do not exist or are empty!\n")
  }
}
my.blat <- function(subject.filename, query.filename, output.filename){
  subject.file.test <- fasta.file.test(subject.filename)
  query.file.test <- fasta.file.test(query.filename)
  if(subject.file.test & query.file.test){
    blat.cmd <- sprintf("blat -t=dna -q=dna -out=blast8 %s %s %s", subject.filename, query.filename, output.filename)
    system(blat.cmd, ignore.stdout = TRUE)
  } else{
    cat("Query and subject files do not exist or are empty!\n")
  }
}

read.blast.output <- function(blast.output.filename, sseq = FALSE, qseq = FALSE){
  blast.output.passes.test <- file.test(blast.output.filename)
  blast.table <- NULL
  colnames.core <- c("seqid", "seqlen", "geneid", "genelen", "frames", "pident", "nident", "length", "mismatch", "gapopen", "seqstart", "seqend", "sstart", "send", "evalue", "bitscore", "qcovs")
  if(sseq & qseq){
    colnames.final <- c(colnames.core, "sseq", "qseq")
  } else if(sseq & !qseq){
    colnames.final <- c(colnames.core, "sseq")
  } else if(!sseq & qseq){
    colnames.final <- c(colnames.core,"qseq")
  } else{
    colnames.final <- colnames.core
  }
  if(blast.output.passes.test){
    blast.table <- read.table(blast.output.filename, stringsAsFactors=F, comment.char="")
    colnames(blast.table) <- colnames.final
  }
  return(blast.table)
}

read.blat.output <- function(blat.output.filename, e.value){
  blat.output.filename <- file.test(blat.output.filename)
  blat.table <- NULL
  if(blat.output.filename){
    blat.table <- read.table(blat.output.filename, stringsAsFactors=F)
    colnames(blat.table) <- c("seqid", "geneid", "pident", "length", "mismatch", "gapopen", "seqstart", "seqend", "sstart", "send", "evalue", "bitscore")
    blat.table <- blat.table[blat.table$evalue<=e.value,]
  } else{
    cat("Blat output does not exist or is empty!\n")
  }
  return(blat.table)
}

calculate.reference.coverage <- function(all.reference.names, blast.output){
  blast.output.nonempty <- nrow(blast.output)>0
  match.between.files <- any(blast.output$seqid %in% all.reference.names)
  matching.results <- NULL
  N <- length(all.reference.names)
  if(blast.output.nonempty & match.between.files){
    ref.covs <- sapply(1:N, function(reference.index){
      ref.name <- all.reference.names[reference.index]
      any.match <- ref.name %in% blast.output$seqid
      ref.coverage <- no.contigs <- seqlen <- 0
      if(any.match){
        df <- blast.output[which(blast.output$seqid == ref.name),]
        all.pos <- unlist(lapply(1:nrow(df), function(k) df$seqstart[k]:df$seqend[k]))
        ref.coverage <- length(unique(all.pos))/df$seqlen[1]
        no.contigs <- length(unique(df$geneid))
        seqlen <- df$seqlen[1]
      }
      c(ref.coverage,no.contigs,seqlen)
    })
    matching.results <- data.frame(reference.name=all.reference.names, coverage=round(ref.covs[1,]*100), no.contigs=ref.covs[2,], ref.length = ref.covs[3,], stringsAsFactors = F)
    matching.results <- matching.results[order(-matching.results$coverage, -matching.results$ref.length),]
    rownames(matching.results) <- NULL
  } else{
    cat("No blast output or wrong references used!\n")
  }
  return(matching.results)
}

remove.likely.transposons <- function(best.match.data){
  # if there are two overlapping matches to a reference but one is longer, then ignore the shorter one
  # this is to ignore short matches to duplicated genes (likely transposons)
  if(nrow(best.match.data)>1){
    best.match.data <- best.match.data[order(-best.match.data$length), ]
    best.match.data.rows.are.redundant <- sapply(2:nrow(best.match.data), function(hit.index){
      previous.rows <- 1:(hit.index-1)
      this.hit.range <- best.match.data$seqstart[hit.index]:best.match.data$seqend[hit.index]
      is.redundant <- FALSE
      overlaps.with.previous <- sapply(previous.rows,function(previous.row.index){
        this.previous.row <- previous.rows[previous.row.index]
        this.previous.hit.range <- best.match.data$seqstart[this.previous.row]:best.match.data$seqend[this.previous.row]
        prop.overlap <- length(intersect(this.previous.hit.range, this.hit.range))/length(this.hit.range)
        prop.overlap
      })
      if(any(overlaps.with.previous>cds.overlap.threshold)) is.redundant <- TRUE
      is.redundant
    })
    best.match.data.rows.are.redundant <- c(FALSE,best.match.data.rows.are.redundant)
    best.match.data <- best.match.data[!best.match.data.rows.are.redundant,]
  }
  return(best.match.data)
}

bm.data.to.locus.region <- function(best.match.data){
  locus.regions <- NULL
  if(nrow(best.match.data)>0){
  locus.regions <- best.match.data[,colnames(best.match.data) %in% c("geneid", "sstart", "send")]
  locus.regions$length <- abs(locus.regions$send-locus.regions$sstart+1)
  locus.regions$strand <- rep("+", nrow(locus.regions))
  locus.regions$strand[locus.regions$sstart>locus.regions$send] <- "-"
  } else{
    cat("bm.data is an empty object!\n")
  }
  return(locus.regions)
}

get.bm.hits <- function(gene.blast.output, best.match.cds.n, locus.regions){
  hits.to.bm <- NULL
  hits.to.bm.exist <- any(best.match.cds.n %in% gene.blast.output$seqid)
  if(hits.to.bm.exist){
    hits.to.bm <- gene.blast.output[gene.blast.output$seqid %in% best.match.cds.n,]
    hits.to.bm$klocus.overlap <- sapply(1:nrow(hits.to.bm), function(al.index){
      # cat(al.index,"\n")
      out <- 0
      # hits.to.bm[al.index,]
      this.contig <- hits.to.bm$geneid[al.index]
      this.range <- hits.to.bm$sstart[al.index]:hits.to.bm$send[al.index]
      relevant.bm.data <- locus.regions[locus.regions$geneid==this.contig,]
      if(nrow(relevant.bm.data)>0){
        relevant.bm.data.pos <- lapply(1:nrow(relevant.bm.data), function(k) relevant.bm.data$sstart[k]:relevant.bm.data$send[k])
        relevant.bm.data.pos <- unique(unlist(relevant.bm.data.pos))
        out <- length(intersect(this.range, relevant.bm.data.pos))/length(this.range)
      }
      out
    })
    hits.to.bm$bm.cds.overlap <- rep(NA, nrow(hits.to.bm))  
  }
  return(hits.to.bm)
}

get.other.hits <- function(gene.blast.output, best.match.cds.n, locus.regions, hits.to.bm){
  hits.to.others <- gene.blast.output[!gene.blast.output$seqid %in% best.match.cds.n,]
  hits.to.others.dupl.test <- hits.to.others[,colnames(hits.to.others) %in% c("geneid", "sstart", "send")]
  hits.to.others <- hits.to.others[!duplicated(hits.to.others.dupl.test),]
  if(nrow(hits.to.others)>0){
    hits.to.others$klocus.overlap <- sapply(1:nrow(hits.to.others), function(al.index){
      out <- 0
      this.contig <- hits.to.others$geneid[al.index]
      this.range <- hits.to.others$sstart[al.index]:hits.to.others$send[al.index]
      relevant.bm.data <- locus.regions[locus.regions$geneid==this.contig,]
      if(nrow(relevant.bm.data)>0){
        relevant.bm.data.pos <- lapply(1:nrow(relevant.bm.data), function(k) relevant.bm.data$sstart[k]:relevant.bm.data$send[k])
        relevant.bm.data.pos <- unique(unlist(relevant.bm.data.pos))
        out <- length(intersect(this.range, relevant.bm.data.pos))/length(this.range)
      }
      out
    })
    
    hits.to.others$bm.cds.overlap <- sapply(1:nrow(hits.to.others), function(al.index){
      out <- 0
      this.contig <- hits.to.others$geneid[al.index]
      this.range <- hits.to.others$sstart[al.index]:hits.to.others$send[al.index]
      relevant.bm.data <- hits.to.bm[hits.to.bm$geneid==this.contig,]
      if(nrow(relevant.bm.data)>0){
        relevant.bm.data.pos <- lapply(1:nrow(relevant.bm.data), function(k) relevant.bm.data$sstart[k]:relevant.bm.data$send[k])
        relevant.bm.data.pos <- unique(unlist(relevant.bm.data.pos))
        out <- length(intersect(this.range, relevant.bm.data.pos))/length(this.range)  
      }
      out
    })
  } else{
    hits.to.others <- NULL
  }
  return(hits.to.others)
}
risk.of.contamination <- function(ref.coverage.table, min.coverage, contamination.thr){
  # Of all references which span a min.coverage of the assembly in question (ie all close matches),
  # test for contamination by checking if other matches are similar; 
  # if not, then there's a risk of contamination.
  # Similarity is defined by jaccard distance less or equal to contamination.thr
  out <- FALSE
  shortlist <- ref.coverage.table[ref.coverage.table$coverage>=min.coverage,]
  if(nrow(shortlist)>1){
    no.pairs <- nrow(shortlist)-1
    pair.distances <- sapply(1:no.pairs, function(pair.index){
      best.match <- shortlist$reference.name[1]
      query <- shortlist$reference.name[1+pair.index]
      match.ab <- ref.jaccard$from==best.match & ref.jaccard$to==query
      match.ba <- ref.jaccard$to==best.match & ref.jaccard$from==query
      ref.jaccard$jaccard.dist[which(match.ab | match.ba)]
    })
    if(any(pair.distances>contamination.thr)) out <- TRUE
  }
  out
}

get.ref.output <- function(hits.to.bm, bm.cds.n){
  n <- length(bm.cds.n)
  # output.reference 
  ref.output.table <- NULL
  for(ref.gene.index in 1:n){
    this.gene.name <- bm.cds.n[ref.gene.index]
    this.gene.found <- this.gene.name %in% hits.to.bm$seqid
    if(this.gene.found){
      this.gene.index <- which(hits.to.bm$seqid == this.gene.name)
      this.gene.pident <- hits.to.bm$pident[this.gene.index]
      this.gene.qcovs <- hits.to.bm$qcovs[this.gene.index]
      this.gene.contig <- hits.to.bm$geneid[this.gene.index]
      this.gene.start <- hits.to.bm$sstart[this.gene.index]
      this.gene.end <- hits.to.bm$send[this.gene.index]
      this.gene.strand <- "+"
      if(this.gene.start>this.gene.end) this.gene.strand <- "-"
      this.gene.assembly.index <- which(assembly.seq.n == this.gene.contig)
      this.gen.seq <- assembly.seq[[this.gene.assembly.index]][this.gene.start:this.gene.end]
      if(this.gene.strand=="-") this.gen.seq <- comp(this.gen.seq)
      df <- data.frame(category="ref", name=this.gene.name, pident=this.gene.pident, qcovs=this.gene.qcovs, 
                       contig.name=this.gene.contig, contig.start=this.gene.start, contig.end=this.gene.end, strand=this.gene.strand, seq=c2s(this.gen.seq), stringsAsFactors = F)
    } else{
      df <- data.frame(category="ref", name=this.gene.name, pident=NA, qcovs=NA, 
                       contig.name=NA, contig.start=NA, contig.end=NA, strand=NA, seq=NA, stringsAsFactors = F) 
    }
    ref.output.table <- rbind(ref.output.table, df)
  } 
  return(ref.output.table)
}

get.extra.output <- function(hits.inside.extra){
  extra.output.table <- NULL
  input.ok <- data.frame.nonempty(hits.inside.extra)
  if(input.ok){
    for(row.index in 1:nrow(hits.inside.extra)){
      df <- data.frame(category="extra", name=hits.inside.extra$seqid[row.index], pident=hits.inside.extra$pident[row.index], qcovs=hits.inside.extra$qcovs[row.index], 
                       contig.name=hits.inside.extra$geneid[row.index], contig.start=hits.inside.extra$sstart[row.index], contig.end=hits.inside.extra$send[row.index], 
                       strand="+", seq=NA, stringsAsFactors = F)
      if(df$contig.start>df$contig.end) df$strand <- "-"
      this.gene.assembly.index <- which(assembly.seq.n == df$contig.name)
      this.gen.seq <- assembly.seq[[this.gene.assembly.index]][df$contig.start:df$contig.end]
      if(df$strand=="-") this.gen.seq <- comp(this.gen.seq)
      df$seq <- paste(this.gen.seq, collapse="")
      extra.output.table <- rbind(extra.output.table, df)
    }  
  } 
  return(extra.output.table)
}

save.locus.region.sequences <- function(locus.regions, assembly.seq, file.output){
  assembly.seq.n <- names(assembly.seq)
  locus.region.sequences <- lapply(1:nrow(locus.regions), function(region.index){
    this.contig.name <- locus.regions$geneid[region.index]
    this.start <- locus.regions$sstart[region.index]
    this.end <- locus.regions$send[region.index]
    this.strand <- locus.regions$strand[region.index]
    this.seq <- assembly.seq[[which(assembly.seq.n == this.contig.name)]][this.start:this.end]
    if(this.strand=="-") this.seq <- comp(this.seq)
    this.seq
  })
  locus.region.names <- sapply(1:nrow(locus.regions), function(region.index){
    this.contig.name <- locus.regions$geneid[region.index]
    this.strand <- locus.regions$strand[region.index]
    this.start <- locus.regions$sstart[region.index]
    this.end <- locus.regions$send[region.index]
    if(this.strand=="-"){
      this.start <- locus.regions$send[region.index]
      this.end <- locus.regions$sstart[region.index]
    }
    sprintf("%s_%d_to_%d_%s_strand", this.contig.name, this.start, this.end, this.strand)
  })
  write.fasta(locus.region.sequences, names=locus.region.names, file.output)
  return(locus.region.names)
}

edit.gene.names <- function(gene.names){
  no.genes <- length(gene.names)
  gene.names.new <- character(0)
  if(no.genes>0){
    gene.names.new <- sapply(1:no.genes, function(k) strsplit(gene.names[k],"_")[[1]][1])
    gene.names.new.dupl <- duplicated(gene.names.new)
    if(any(gene.names.new.dupl)){
      dupl.names <- unique(gene.names.new[which(gene.names.new.dupl)])
      no.dupl <- length(dupl.names)
      for(i in 1:no.dupl){
        this.dupl <- dupl.names[i]
        this.dupl.index <- which(gene.names.new==this.dupl)
        gene.names.new[this.dupl.index] <- sprintf("%s_%d", this.dupl, 1:length(this.dupl.index))
      }
    }
  }
  return(gene.names.new)
}

gene.table.from.blast <- function(this.blast.output){
  blast.table <- NULL
  if(nrow(this.blast.output)>0){
    start.end.matrix <- sapply(1:nrow(this.blast.output), function(k){
      start.end <- c(this.blast.output$sstart[k], this.blast.output$send[k])
      c(min(start.end), max(start.end))
    })
    this.blast.output$sstart.norm <- start.end.matrix[1,]
    this.blast.output$send.norm <- start.end.matrix[2,]
    
    no.blast.hits <- nrow(this.blast.output)
    if(no.blast.hits>1){
      overlapping.matrix <- matrix(0, nrow=no.blast.hits, ncol=no.blast.hits)
      for(index.i in 1:(no.blast.hits-1)){
        for(index.j in (index.i+1):no.blast.hits){
          start.norm.i <- this.blast.output$sstart.norm[index.i]
          end.norm.i <- this.blast.output$send.norm[index.i]
          start.norm.j <- this.blast.output$sstart.norm[index.j]
          end.norm.j <- this.blast.output$send.norm[index.j]
          hits.are.duplicates <- are.genes.same(start.norm.i, end.norm.i, start.norm.j, end.norm.j, cds.overlap.threshold)
          if(hits.are.duplicates) overlapping.matrix[index.i, index.j] <- 1    
        }
      }
      g <- graph_from_adjacency_matrix(overlapping.matrix)
      clustering.g <- clusters(g)
      no.blast.cds <- clustering.g$no
      cluster.members <- lapply(1:no.blast.cds, function(k) which(clustering.g$membership==k))
      
      full.table <- NULL
      for(blast.cds.index in 1:no.blast.cds){
        this.cds.blast.hits <- this.blast.output[cluster.members[[blast.cds.index]],]
        this.cds.blast.hits.sorted <- this.cds.blast.hits[order(-this.cds.blast.hits$bitscore, -this.cds.blast.hits$qcovs),]
        # this.cds.blast.hits.sorted <- this.cds.blast.hits[order(-this.cds.blast.hits$qcovs, -this.cds.blast.hits$pident),]
        this.cds.df <- data.frame(contig.name = this.cds.blast.hits.sorted$geneid[1],
                                  start.pos = this.cds.blast.hits.sorted$sstart[1],
                                  end.pos = this.cds.blast.hits.sorted$send[1],
                                  strand = get.strand(this.cds.blast.hits.sorted$sstart[1], this.cds.blast.hits.sorted$send[1]),
                                  locus.tag = "blast",
                                  product.name = this.cds.blast.hits.sorted$gn[1],
                                  stringsAsFactors = F)
        full.table <- rbind(full.table, this.cds.df)
      }
      full.table <- full.table[order(full.table$start.pos),]
    } else{
      full.table <- data.frame(contig.name = this.blast.output$geneid[1],
                                start.pos = this.blast.output$sstart[1],
                                end.pos = this.blast.output$send[1],
                                strand = get.strand(this.blast.output$sstart[1], this.blast.output$send[1]),
                                locus.tag = "blast",
                                product.name = this.blast.output$gn[1],
                                stringsAsFactors = F)
    }
    
    full.table$product.name <- edit.gene.names(full.table$product.name)
    blast.table <- full.table
  } 
  return(blast.table)
}

combine.prokka.blast <- function(prokka.table, blast.table){
  combined.table <- prokka.table
  if(is.null(blast.table)){
    no.blast.hits <- 0
  } else{
    no.blast.hits <- nrow(blast.table)
  }
  
  if(no.blast.hits>0){
    # blast.table.gn.raw <- sapply(1:nrow(blast.table), function(k) strsplit(blast.table$product.name[k], "_")[[1]][1])
    
    ### 1. find CDS which overlap with multiple Prokka hits:
    ## a) are there genes 
    # duplicated(blast.table.gn.raw)
    
    ### 2. find CDS which Prokka did not find
    blast.table.prokka.overlap <- sapply(1:no.blast.hits, function(blast.gene.index){
      this.gene.start.pos <- blast.table$start.pos[blast.gene.index]
      this.gene.end.pos <- blast.table$end.pos[blast.gene.index]
      overlap.with.prokka <- sapply(1:nrow(prokka.table), function(k) are.genes.same(this.gene.start.pos, this.gene.end.pos, prokka.table$start.pos[k], prokka.table$end.pos[k], cds.overlap.threshold))
      any(overlap.with.prokka)  
    })
    if(any(!blast.table.prokka.overlap)){
      combined.table <- rbind(prokka.table, blast.table[which(!blast.table.prokka.overlap),])
    }
    combined.table <- combined.table[order(combined.table$start.pos),]  
  }
  return(combined.table)
}

