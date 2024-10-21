# augment_phased_map_one.R

augment_phased_map_one <- function(map, mrk, mat, geno, max.phases,
                                   ploidy.p1, ploidy.p2, dosage.p1,
                                   dosage.p2, tol, thresh.LOD.ph.to.insert,
                                   thresh.rf.to.insert, verbose, n.ind) {
  # Extract positioned and unpositioned markers
  mrk.pos <- rownames(map$p1)  # Positioned markers
  mrk.id <- setdiff(mrk, mrk.pos)  # Markers to be positioned

  # Check if there are any markers to be positioned
  if(length(mrk.id) == 0) {
    return(map)  # No unpositioned markers found
  }

  # Function to perform two-point phasing for a parent
  perform_phasing <- function(dosage, map, M, parent, mrk.id, mrk.pos) {
    dose.vec <- dosage[mrk.id]
    InitPh <- map
    S <- M[[paste0("Sh.", parent)]][mrk.id, mrk.pos, drop = FALSE]
    phasing_one(mrk.id, dose.vec, S, InitPh, verbose = FALSE)
  }

  # Two-point phasing for both parents
  L1 <- perform_phasing(dosage.p1, map$p1, mat, "p1", mrk.id, mrk.pos)
  L2 <- perform_phasing(dosage.p2, map$p2, mat, "p2", mrk.id, mrk.pos)

  # Selecting phase configurations
  n.conf <- sapply(L1, nrow) * sapply(L2, nrow)
  if(verbose) {
    cat("Distribution of phase configurations:\n")
    temp <- as.data.frame(table(n.conf))
    colnames(temp) <- c("n. phase. conf.", "frequency")
    print(temp)
  }

  # Handle cases with no selected markers
  mrk.sel <- which(n.conf <= max.phases)
  if(length(mrk.sel) == 0) {
    warning("No markers were selected for 'max.phases' = ", max.phases,
            "\n increasing 'max.phases' to ", min(n.conf) + 1)
    max.phases <- min(n.conf) + 1
    mrk.sel <- which(n.conf <= max.phases)
  }

  # Update phase information for selected markers
  L1 <- L1[mrk.sel]
  L2 <- L2[mrk.sel]
  mrk.id <- mrk.id[mrk.sel]

  # Create pedigree matrix
  pedigree <- matrix(rep(c(1, 2, ploidy.p1, ploidy.p2, 1), n.ind),
                     nrow = n.ind, byrow = TRUE)

  # Find flanking markers
  flanking <- find_flanking_markers(mrk, mrk.pos, mrk.id)
  phasing_results <- vector("list", length(flanking))
  names(phasing_results) <- names(flanking)

  # Initialize progress bar if verbose mode is enabled
  if(verbose) pb <- utils::txtProgressBar(min = 0, max = length(L1), style = 3)

  # Iterating over each set of phasing results
  for(j in seq_along(L1)) {
    # Extract genotype data for the current marker
    G <- geno[mrk.id[j], , drop = TRUE]
    u <- match(unlist(flanking[[mrk.id[j]]]), mrk.pos)

    # Determine the position of the marker and set homolog probabilities
    if(is.na(u)[1]) {  # Marker at the beginning of the linkage group
      homolog_prob <- as.matrix(map$haploprob[, c(na.omit(u), na.omit(u) + 1) + 3])
      idx <- c(0, 1, 2)
    } else if(is.na(u)[2]) {  # Marker at the end of the linkage group
      homolog_prob <- as.matrix(map$haploprob[, c(na.omit(u) - 1, na.omit(u)) + 3])
      idx <- c(0, 2, 1)
    } else {  # Marker in the middle of the linkage group
      homolog_prob <- as.matrix(map$haploprob[, u + 3])
      idx <- c(1, 0, 2)
    }

    # Initialize variables for phasing computations
    w2 <- w1 <- NULL
    z <- vector("list", nrow(L1[[j]]) * nrow(L2[[j]]))
    count <- 1

    # Nested loops for computing phasing results
    for(l in seq_len(nrow(L1[[j]]))) {
      for(k in seq_len(nrow(L2[[j]]))) {
        PH <- list(L1[[j]][l, ], L2[[j]][k, ])
        z[[count]] <- est_hmm_map_biallelic_insert_marker(PH, G, pedigree, homolog_prob,
                                                          rf = c(0.01, 0.01), idx, verbose = FALSE,
                                                          detailed_verbose = FALSE, tol = tol, ret_H0 = FALSE)
        w1 <- rbind(w1, L1[[j]][l, ])
        w2 <- rbind(w2, L2[[j]][k, ])
        count <- count + 1
      }
    }

    # Process the phasing results
    v <- sapply(z, function(x) x[[1]])
    v <- max(v) - v
    id <- order(v)
    phasing_results[[mrk.id[j]]] <- list(loglike = v[id],
                                         rf.vec = t(sapply(z[id], function(x) x[[2]])),
                                         phases = list(p1 = w1[id, , drop = FALSE],
                                                       p2 = w2[id, , drop = FALSE]))

    # Update progress bar if verbose mode is enabled
    if(verbose) utils::setTxtProgressBar(pb, j)
  }

  # Close the progress bar if verbose mode is enabled
  if(verbose) close(pb)

  # Selecting the list of phasing results based on loglike criteria
  selected.list <- phasing_results[sapply(phasing_results, function(x) {
    length(x$loglike) == 1 || x$loglike[2] > thresh.LOD.ph.to.insert
  })]

  # Set default value for thresh.rf.to.insert if it is NULL
  if(is.null(thresh.rf.to.insert)) {
    thresh.rf.to.insert <- max(map$rf)
  }

  # Validate thresh.rf.to.insert value
  if(thresh.rf.to.insert < 0 || thresh.rf.to.insert >= 0.5) {
    stop("'thresh.rf.to.insert' parameter must be between 0 and 0.5")
  }

  # Further filtering of phasing results based on recombination frequency
  selected.list <- phasing_results[sapply(phasing_results, function(x) {
    max(x$rf.vec[1,]) <= thresh.rf.to.insert
  })]

  # Iterate over selected list to update map information
  for(j in names(selected.list)) {
    cur.mrk <- rownames(map$p1)
    pos <- find_flanking_markers(mrk, cur.mrk, j)

    # Skip if no flanking markers are found
    if(length(unlist(pos)) == 0) {
      next()
    }

    # Inserting markers at the beginning, end, or middle of the linkage group
    if(is.na(pos[[1]]$preceding)) {  # Beginning
      map$p1 <- rbind(selected.list[[j]]$phases$p1[1,], map$p1)
      map$p2 <- rbind(selected.list[[j]]$phases$p2[1,], map$p2)
      rownames(map$p1) <- rownames(map$p2) <- c(j, cur.mrk)
    } else if(is.na(pos[[1]]$succeeding)) {  # End
      map$p1 <- rbind(map$p1, selected.list[[j]]$phases$p1[1,])
      map$p2 <- rbind(map$p2, selected.list[[j]]$phases$p2[1,])
      rownames(map$p1) <- rownames(map$p2) <- c(cur.mrk, j)
    } else {  # Middle
      preceding <- cur.mrk[1:match(pos[[1]]$preceding, cur.mrk)]
      succeeding <- cur.mrk[(match(pos[[1]]$succeeding, cur.mrk)):length(cur.mrk)]
      map$p1 <- rbind(map$p1[preceding,], selected.list[[j]]$phases$p1[1,], map$p1[succeeding,])
      map$p2 <- rbind(map$p2[preceding,], selected.list[[j]]$phases$p2[1,], map$p2[succeeding,])
      rownames(map$p1) <- rownames(map$p2) <- c(preceding, j, succeeding)
    }
  }

  # Reset potentially redundant fields
  map$loglike <- NULL
  map$rf <- NULL
  map$error <- NULL
  map$haploprob <- NULL

  return(map)
}
