#' @export
print.mappoly2.data <- function(x,
                                type = c("screened", "raw"),
                                detailed = FALSE,
                                ...) {
  #### Raw ####
  msg("Data summary", col = "blue")
  txt <- list(
    paste0("    Ploidy level of ", x$name.p1, ":"),
    paste0("    Ploidy level of ", x$name.p2, ":"),
    "    No. individuals:",
    "    No. unique markers:",
    "    Percentage of redundant:",
    "    Percentage of missing:",
    "    Chromosome info:",
    "    Genome position:",
    "    Recombination farction:",
    "    Marker scope:")
  n <- sapply(txt, nchar)
  for (i in 1:length(txt)) {
    txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
  }
  id <- is.na(x$geno.dose)
  cat(" ", txt[[1]], " ", x$ploidy.p1, sep ="")
  cat("\n", txt[[2]], x$ploidy.p2)
  cat("\n", txt[[3]], x$n.ind)
  cat("\n", txt[[4]], length(x$mrk.names))


  if (is.null(x$redundant))
    cat("\n", txt[[5]], "unavailable")
  else if (is.numeric(x$redundant)) {
    cat("\n", txt[[5]], "0%")
  } else {
    cat("\n ", txt[[5]], " ",  round(100*nrow(x$redundant)/(length(x$mrk.names) + nrow(x$redundant)),1), "%" ,sep = "")
  }
  cat("\n ", txt[[6]], " ",   round(100*sum(id)/length(id),1), "%" ,sep = "")
   chrom.flag <- FALSE
  if (all(is.null(x$chrom)) || all(is.na(x$chrom)))
    cat("\n", txt[[7]], "unavailable")
  else{
    cat("\n", txt[[7]], "available")
    chrom.flag <- TRUE
  }
  if(any(is.numeric(x$genome.pos)))
    cat("\n",  txt[[8]], "available")
  else
    cat("\n",  txt[[8]], "unavailable")

  #### RF ####
  if(has.mappoly2.rf(x)){
    cat("\n",  txt[[9]], "available")
    cat("\n ",  txt[[10]], " ",
        x$pairwise.rf$mrk.scope,
        " (",
        paste(unique(x$chrom[colnames(x$pairwise.rf$rec.mat)]), collapse = ", "),
        ")\n", sep = "")
  } else
    cat("\n",  txt[[9]], "unavailable\n")
  #### Screened ####
  if(has.mappoly2.screened(x)){
    msg("Filtering information", col = "blue")
    txt <- list(
      "    Thresholds",
      "       missing mrk: ",
      "       missing ind: ",
      "       chi-square pval: ",
      "       read depth: ",
      "       non full-sib: ",
      "       LOD phase: ",
      "       LOD rf: ",
      "       rf: ",
      "       prob.lower: ",
      "       prob.upper: ",
      "    Screened mrk: ",
      "    Screened ind: ")
    n <- sapply(txt, nchar)
    for (i in 1:length(txt)) {
      txt[[i]] <- paste(txt[[i]], paste0(rep(" ", max(n) - n[i]), collapse = ""))
    }
    id <- is.na(x$geno.dose)
    cat(txt[[1]])
    cat("\n", txt[[2]], x$screened.data$thresholds$miss.mrk)
    cat("\n", txt[[3]], x$screened.data$thresholds$miss.ind)
    cat("\n", txt[[4]], format(x$screened.data$thresholds$chisq.pval, digits = 3))
    if(!all(is.na(x$QAQC.values$markers$read.depth)))
      cat("\n", txt[[5]], x$screened.data$thresholds$read.depth)
    cat("\n", txt[[6]], ifelse(all(is.na(x$QAQC.values$individuals$full.sib)),
                               "-", sum(!x$QAQC.values$individuals$full.sib)))
    cat("\n", txt[[7]], x$screened.data$thresholds$LOD.ph)
    cat("\n", txt[[8]], x$screened.data$thresholds$LOD.rf)
    cat("\n", txt[[9]], x$screened.data$thresholds$rf)
    cat("\n", txt[[10]], x$screened.data$thresholds$prob.lower)
    cat("\n", txt[[11]], x$screened.data$thresholds$prob.upper)
    cat("\n", txt[[12]], length(x$screened.data$mrk.names))
    cat("\n", txt[[13]], length(x$screened.data$ind.names), "\n")
  }
  #### Detailed ####
  w <- table(x$chrom, useNA = "always")
  w <- w[order(embedded_to_numeric(names(w)))]
  names(w)[is.na(names(w))] <- "NoCrh"
  if(detailed){
    if(chrom.flag){
      msg("No. markers per chromosome", col = "blue")
      print(data.frame(chrom = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
    }
    msg("No. of markers per dosage in both parents", col = "blue")
    freq <- table(paste(x$dosage.p1,
                        x$dosage.p2, sep = "-"))
    d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
    d.temp <- data.frame(paste0("    ", d.temp[, 1]),
                         d.temp[, 2],
                         as.numeric(freq))
    colnames(d.temp) <- c(x$name.p1, x$name.p2, "freq")
    print(d.temp, row.names = FALSE)
  }
  msg("", col = "blue")
}
#' @export
print.mappoly2.group <- function(x,
                                 detailed = TRUE,
                                 ...) {
  msg("Grouping summary", col = "blue")
  ## criteria
  cat("       - Number of linkage groups:  ", length(unique(x$groups.snp)), "\n")
  cat("       - Markers per linkage groups: \n")
  w <- table(x$groups.snp, useNA = "ifany")
  w <- data.frame(group = names(w), n_mrk = as.numeric(w), row.names = NULL)
  print_matrix(mat = w, 8, row.names = FALSE)
  cat("\n")
  ## printing summary
  if(!is.null(x$seq.vs.grouped.snp)){
    print_matrix(mat = x$seq.vs.grouped.snp, 8)
  }
  msg("", col = "blue")
}

pad_strings <- function(vec, n) {
  # Find the maximum length of strings in the vector without spaces
  max_length_no_spaces <- max(nchar(gsub(" ", "", vec)))

  # Calculate the desired total length
  desired_length <- max_length_no_spaces + n

  # Pad each string to the desired length
  z<-sapply(vec, function(x) {
    # Remove trailing spaces and then pad
    trimmed_x <- gsub("\\s+$", "", x)
    sprintf("%-*s", desired_length, trimmed_x)
  })
  names(z) <- NULL
  z
}

print_matrix <- function(mat,
                         spaces = 5,
                         zero.print = ".",
                         row.names = TRUE,
                         header = FALSE,
                         equal.space = TRUE,
                         footer = FALSE,
                         title=NULL){
  mat[mat==0]<-zero.print
  txt1 <- NULL
  for(i in 1:ncol(mat))
    txt1 <- c(txt1, colnames(mat)[i], mat[,i])
  n1 <- nchar(txt1)
  for (i in 1:length(txt1))
    txt1[i] <- paste(txt1[i], paste0(rep(" ", max(n1) - n1[i]), collapse = ""))
  dim(txt1) <- c(nrow(mat)+1, ncol(mat))
  if(!equal.space)
    for(i in 1:ncol(txt1))
      txt1[,i] <- pad_strings(txt1[,i], 1)
  txt2 <- rownames(mat)
  n2 <- nchar(txt2)
  for (i in 1:length(txt2))
    txt2[i] <- paste(txt2[i], paste0(rep(" ", max(n2) - n2[i]), collapse = ""))
  txt2 <- c(paste0(rep(" ", (max(n2)+1)), collapse = ""), txt2)
  for (i in 1:length(txt2))
    txt2[i] <- paste0(paste0(rep(" ", spaces), collapse = ""), txt2[i])
  txt3 <- paste0(paste0(rep(" ", spaces), collapse = ""),
                 paste0(rep("-", sum(nchar(c(txt2[1], txt1[1,])))+ncol(mat)-spaces), collapse = ""))
  if(!is.null(title)){
    cat(paste0(paste0(rep(" ", spaces), collapse = ""),
               paste0(rep("-", 10), collapse = ""),
               paste0(" ", title, " "),
               paste0(rep("-", nchar(txt3) - 12 - nchar(title)), collapse = "")), "\n")
  }
  cat(txt2[1], txt1[1,], "\n")
  cat(txt3, "\n")
  for(i in 2:nrow(txt1)){
    if(header & i == 3)
      cat(txt3, "\n")
    if(row.names) cat(txt2[i], txt1[i,], "\n")
    else cat(paste0(rep(" ", nchar(txt2[i])), collapse = ""), txt1[i,], "\n")
    if(footer & (i == (nrow(txt1)-1))){
      cat(txt3, "\n")
    }
  }

  cat(txt3, "\n")
  invisible(nchar(txt3))
}

get_sequence_mat <- function(x,
                             type,
                             parent){
  dn <- list(c(names(x$maps)),
             c("Ch","n.mrk", "ord", "phase", "map (cM)", "haplo"))
  Chrom <- sapply(x$maps, function(y) paste0(embedded_to_numeric(unique(x$data$chrom[y[[type]]$mkr.names])), collapse = "/"))
  n.mrk <- sapply(x$maps, function(y) length(y[[type]]$mkr.names))
  ord <- sapply(x$maps, function(y) ifelse(is.null(y[[type]]$order), "N", "Y"))
  u1 <- sapply(x$maps, function(z) ifelse(is.null(z[[type]][[parent]][["rf.phase"]]), 0,
                                          length(z[[type]][[parent]][["rf.phase"]])))
  u2 <- sapply(x$maps, function(z) ifelse(is.null(z[[type]][[parent]][["rf.phase"]]), 0,
                                          nrow(z[[type]][[parent]][["rf.phase"]][[1]]$p1)))
  u2 <- round((100*u2/n.mrk), 1)
  phase <- paste0(u1, "/", paste0(u2,"%"), sep = "")
  map.len <- sapply(x$maps, function(z) ifelse(is.null(z[[type]][[parent]][["hmm.phase"]][[1]]$rf), 0,
                                               round(sum(imf_h(z[[type]][[parent]][["hmm.phase"]][[1]]$rf)),1)))
  map.len <- format(map.len, digits = 2)
  haplo <- sapply(x$maps, function(z) ifelse(is.null(z[[type]][[parent]][["hmm.phase"]][[1]]$haploprob), "N", "Y"))
  M <- cbind(Chrom, n.mrk, ord, phase, map.len, haplo)
  dimnames(M) <- dn
  return(M)
}

#' @export
print.mappoly2.sequence <- function(x,
                                    type = c("mds", "genome", "data"), ...){
  type <- match.arg(type)
  if(type == "data")
    print(x$data)
  mat.p1 <- get_sequence_mat(x, type, "p1")
  mat.p2 <- get_sequence_mat(x, type, "p2")[,4:6,drop = FALSE]
  mat.p1p2 <- get_sequence_mat(x, type, "p1p2")[,4:6,drop = FALSE]
  mat<-cbind(mat.p1,mat.p2,mat.p1p2)
  mat <- rbind(colnames(mat), mat)
  colnames(mat) <- c("", "", "", "p1","","", "p2", "", "", "p1p2","", "")
  {
    if(type == "mds")
      print_matrix(mat, spaces = 0, zero.print = ".",
                   row.names = TRUE,  equal.space = FALSE,
                   header = TRUE, title = "MDS")
    else if(type == "genome" )
      print_matrix(mat, spaces = 0, zero.print = ".",
                   row.names = TRUE,  equal.space = FALSE,
                   header = TRUE, title = "Genome")

  }
  invisible(mat)
}

#' Summarize Genetic Mapping Data
#'
#' This function provides a summary of genetic mapping data, including information
#' about marker distribution, map lengths, and gap sizes. It's designed to work with
#' specific types of genetic data structures in R, offering a concise overview of
#' the mapping results.
#'
#' @param x A genetic mapping data object, typically containing information about
#'          maps, markers, and related genetic data.
#' @param type The type of summary required: 'both', 'mds', or 'genome'.
#'             'mds' refers to Multi-Dimensional Scaling, and 'genome' refers to
#'             the entire genomic data. 'both' will provide information for both
#'             mds and genome.
#' @param parent Specifies which parent's data to use in the summary.
#'               Options are 'p1p2' (both parents), 'p1' (first parent), or 'p2'
#'               (second parent).
#'
#' @return An invisible data frame containing the mapping summary,
#'         including linkage groups (LG), chromosome data (Chrom), map lengths,
#'         number of markers per centiMorgan (Markers/cM), information on
#'         simplex markers for each parent, double-simplex, multiplex, total
#'         marker count, and the maximum gap size.
#'
#' @details The function processes the mapping data based on the specified 'type'
#'          and 'parent' parameters. It calculates various statistics including
#'          the maximum gap size in the map and the distribution of different
#'          types of markers. The output is formatted as a table for easy viewing
#'          and interpretation.
#' @export
map_summary <- function(x,
                        type = c("both", "mds", "genome"),
                        parent = c("p1p2", "p1", "p2")) {

  # Ensure that 'x' is a valid mappoly2 sequence
  if (!is.mappoly2.sequence(x)) {
    stop("The input data is not a valid mappoly2 sequence")
  }
  
  # Select the default type and parent if not specified
  type <- match.arg(type)
  parent <- match.arg(parent)

  # If both types are requested, check if either has been completed
  if (type == "both") {
    completed_summaries <- list()
    
    if ("mds" %in% dimnames(x$maps[[1]])[[2]]) {
      x1 <- tryCatch({
        map_summary(x, type = "mds", parent)
      }, error = function(e) {
        message("Skipping MDS map summary: ", e$message)
        NULL
      })
      completed_summaries$mds <- x1
    }
    
    if ("genome" %in% dimnames(x$maps[[1]])[[2]]) {
      x2 <- tryCatch({
        map_summary(x, type = "genome", parent)
      }, error = function(e) {
        message("Skipping Genome map summary: ", e$message)
        NULL
      })
      completed_summaries$genome <- x2
    }
    
    if (length(completed_summaries) == 0) {
      stop("No valid map summary data available for either MDS or Genome.")
    }
    
    return(invisible(completed_summaries))
  }

  # Detect if the map has been estimated for the given type and parent
  v <- detect_hmm_est_map(x)
  
  # Handle cases where only one type (genome or mds) is computed
  if (!type %in% dimnames(v)[[2]]) {
    stop(paste("Map type", type, "has not been computed for linkage group", parent))
  }
  
  # Extract mapping status and ensure NA values are replaced with FALSE
  u <- apply(v[parent, , , drop = FALSE], 1, function(x) {
    all(!is.na(x)) && all(x)  # Ensure no NA values and all are TRUE
  })
  
  # Ensure 'u' doesn't have NA values by replacing them with FALSE
  u[is.na(u)] <- FALSE
  
  # Identify which map computations are missing
  h <- names(u)[1:2][!u[1:2]]

  # Check for missing map estimation and handle error cases accordingly
  if (length(h) == 1) {
    if (!u[type]) {
      stop(paste(h, "order has not been computed for", parent))
    }
  } else if (length(h) == 2) {
    if (!u[type]) {
      stop(paste(h[1], "and", h[2], "orders have not been computed for", parent))
    }
  }

  # Proceed with generating map summary if no errors occur
  w <- lapply(x$maps, function(y) y[[type]])
  mrk.id <- m <- vector("list", length(w))
  names(mrk.id) <- names(m) <- names(w)
  
  for (i in names(w)) {
    mrk.id[[i]] <- rownames(w[[i]][[parent]]$hmm.phase[[1]]$p1)
    if (is.mapped.sequence(x, i, type, parent)) {
      m[[i]] <- round(imf_h(w[[i]][[parent]]$hmm.phase[[1]]$rf), 1)
    } else {
      m[[i]] <- 0
    }
  }
  
  mg <- sapply(m, max)
  ml <- sapply(m, sum)
  mn <- sapply(mrk.id, function(y) length(y))
  md <- sapply(mrk.id, function(y, x) sapply(get_dosage_type(x, mrk.names = y), length), x)
  md <- cbind(md, apply(md, 1, sum))
  y <- c(round(mn / ml, 3), round(sum(mn / sum(ml))))
  y[is.infinite(y)] <- 0

  # Create a summary table
  mat <- data.frame(
    "LG" = c(names(w), "Total"),
    "Chrom" = c(sapply(mrk.id, function(y) paste0(embedded_to_numeric(unique(x$data$chrom[y])), collapse = "/")), ""),
    "Map_length_(cM)" = round(c(ml, sum(ml)), 1),
    "Markers/cM" = y,
    "Simplex_P1" = md[1, ],
    "Simplex_P2" = md[2, ],
    "Double-simplex" = md[3, ],
    "Multiplex" = md[4, ],
    "Total" = c(mn, sum(mn)),
    "Max_gap" = c(mg, max(mg)),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # Create titles for the output depending on the type (mds or genome)
  p <- c(x$data$name.p1, x$data$name.p2, paste(x$data$name.p1, x$data$name.p2, sep = " x "))
  names(p) <- c("p1", "p2", "p1p2")
  
  if (type == "mds") {
    print_matrix(mat, spaces = 0, zero.print = ".", row.names = FALSE, equal.space = FALSE, header = FALSE, footer = TRUE, title = paste0("MDS --- ", p[parent]))
  } else if (type == "genome") {
    print_matrix(mat, spaces = 0, zero.print = ".", row.names = FALSE, equal.space = FALSE, header = FALSE, footer = TRUE, title = paste0("Genome --- ", p[parent]))
  }

  class(mat) <- "mappoly2.map.summary"  # Assign a custom class to the returned object

  print(mat)  # Automatically print the summary

  invisible(mat)  # Return the summary invisibly
}
