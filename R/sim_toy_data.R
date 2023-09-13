# NAIR: Network Analysis of Immune Repertoire
# Copyright (C) 2023 Li Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

simulateToyData <- function(
    samples = 2,
    chains = 1,
    sample_size = 100,
    prefix_length = 7,
    prefix_chars = c("G", "A", "T", "C"),
    prefix_probs = rbind(
      "sample1" = c(12, 4, 1, 1),
      "sample2" = c(4, 12, 1, 1)),
    affixes = c("AATTGG", "AATCGG", "AATTCG",
                "AATTGC", "AATTG", "AATTC"),
    affix_probs = rbind(
      "sample1" = c(10, 4, 2, 2, 1, 1),
      "sample2" = c(1, 1, 1, 2, 2.5, 2.5)),
    num_edits = 0,
    edit_pos_probs = function(seq_length) {
      stats::dnorm(seq(-4, 4, length.out = seq_length))
    },
    edit_ops = c("insertion", "deletion", "transmutation"),
    edit_probs = c(5, 1, 4),
    new_chars = prefix_chars,
    new_probs = prefix_probs,
    output_dir = NULL,
    no_return = FALSE
) {

  # Initialize storage objects
  seqs <- character(samples * sample_size)
  ranges <- (matrix(seq_len(samples * sample_size), ncol = samples))

  # Define function for performing edit operations, if applicable
  if (num_edits > 0) {
    edit_op <- function(seq, pos, char) {
      op <- sample(edit_ops, size = 1, prob = edit_probs)
      pos2 <- pos + 1; pos3 <- nchar(seq)
      if (op == "deletion") { char <- NULL }
      if (op != "insertion") { pos <- pos - 1 }
      seq <- paste0(substr(seq, 1, pos), char, substr(seq, pos2, pos3))
    }
  }

  for (i in 1:samples) {
    # Randomly generate sequence prefix for each observation,
    # according to specified prefix chars, length and generation probabilities
    prefix <- apply(matrix(sample(prefix_chars,
                                  size = prefix_length * sample_size,
                                  replace = TRUE, prob = prefix_probs[i, ]),
                           ncol = sample_size),
                    MARGIN = 2,
                    FUN = function(x) paste0(x, collapse = ""))
    # Randomly generate sequence affix for each observation
    seqs[ranges[ , i]] <-
      paste0(prefix, sample(affixes, size = sample_size,
                            replace = TRUE, prob = affix_probs[i, ]))
    # Apply randomized edit operations, if specified
    if (num_edits > 0) {
      for (j in 1:num_edits) {
        char_draws <- sample(new_chars, size = sample_size, replace = TRUE,
                             prob = new_probs[i, ])
        edit_pos <- sapply(seqs[ranges[ , i]],
                           function(x) {
                             sample(seq_len(nchar(x)), size = 1,
                                    prob = edit_pos_probs(nchar(x)))
                           })
        tmp <- c(seqs[ranges[ , i]], as.character(edit_pos), char_draws)
        tmp <- matrix(tmp, ncol = 3)
        seqs[ranges[ , i]] <-
          apply(tmp, MARGIN = 1,
                FUN = function(x) { edit_op(x[[1]], as.integer(x[[2]]), x[[3]]) })

      }
    }
  }

  # Generate count data
  counts <- stats::rgamma(samples * sample_size, shape = 100, rate = 1/100)
  counts <- ceiling(max(counts) - counts) + 1
  if (chains == 1) {
    dat <- data.frame(CloneSeq = seqs)
    dat$CloneCount <- dat$CloneFrequency <- counts
    for (i in 1:samples) {
      dat$CloneFrequency[ranges[ , i]] <-
        counts[ranges[ , i]] / sum(counts[ranges[ , i]])
    }
  } else if (chains == 2) {
    alpha_seqs <- beta_seqs <- seqs
    modify_indices <- sample(
      c(TRUE, FALSE), size = sample_size * samples, replace = TRUE)
    beta_seqs[modify_indices] <-
      paste0(beta_seqs[modify_indices], "G")
    dat <- data.frame(AlphaSeq = alpha_seqs, BetaSeq = beta_seqs)
    dat$Count <- counts
    dat$UMIs <- pmax(1, stats::rbinom(sample_size * samples,
                                      size = 100, prob = 0.03))
  } else {
    stop("`chains` must be 1 or 2")
  }
  dat$SampleID <- as.character(rep(paste0("Sample", 1:samples),
                                   each = sample_size))

  # Save to file if specified
  if (!is.null(output_dir)) {
    for (i in 1:samples) {
      tmp <- dat[ranges[ , i], , drop = FALSE]
      rownames(tmp) <- 1:nrow(tmp)
      saveRDS(tmp, file.path(output_dir, paste0("Sample", i, ".rds")))
    }
  }

  if (!no_return) {
    return(dat)
  } else {
    return(TRUE)
  }
}

