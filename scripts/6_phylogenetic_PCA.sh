#!/usr/bin/env bash
# =============================================================================
# 6_phylogenetic_PCA.sh
# Materials and Methods: "Phylogenetic PCA and statistical analyses"
#
# Phylogenetic principal component analysis (pPCA; Revell 2009) of the plant
# cell wall-degrading enzyme (PCWDE) families on the maximum-likelihood
# supermatrix tree, phylogenetic signal of each family and of the pPCA scores,
# and Spearman's rank correlation between BUSCO completeness and the numbers of
# pectin-degrading and total PCWDE genes in the 51 lichenized genomes.
#
# Usage (from the repository root):
#   bash scripts/6_phylogenetic_PCA.sh [OUT_DIR]          # default results/pPCA
#
# Inputs (this repository):
#   data/CAZyme_family_counts.csv           genes per CAZyme family and genome (script 2)
#   data/genome_statistics.csv              BUSCO completeness and PCWDE totals (script 2)
#   data/genomes.tsv                        lifestyle of each genome
#   data/CAZyme_family_groups.tsv           the 49 PCWDE families
#   trees/supermatrix_ML_rooted_94taxa.nwk  rooted maximum-likelihood tree (script 3)
#
# Steps:
#   1. Counts of the 49 PCWDE families for the 94 genomes in the tree; families
#      with non-zero counts in at least four genomes are kept (46 families) and
#      log10(x + 1)-transformed.
#   2. Pagel's lambda and Blomberg's K of every family (phytools::phylosig).
#   3. pPCA on the correlation matrix under a lambda model:
#      phytools::phyl.pca(tree, X, method = "lambda", mode = "corr").
#   4. Pagel's lambda of the scores on PC1 and PC2.
#   5. Spearman's rank correlation, BUSCO completeness vs pectin and PCWDE genes.
#
# Outputs (OUT_DIR): pPCA_scores.csv, pPCA_loadings.csv, pPCA_variance_explained.csv,
#   phylogenetic_signal_per_family.csv, phylogenetic_signal_of_scores.csv,
#   spearman_BUSCO.tsv, pPCA_summary.txt
#   (the signs of the axes are arbitrary; Figure 3 shows PC2 with the opposite sign)
#
# Software: R 4.5.3 with ape 5.8-1 and phytools 2.5-2 (environment: envs/lichen_ppca.yml).
# =============================================================================
set -euo pipefail

OUT_DIR=${1:-results/pPCA}
RSCRIPT=${RSCRIPT:-Rscript}

"$RSCRIPT" --vanilla - data/CAZyme_family_counts.csv data/genome_statistics.csv data/genomes.tsv \
        data/CAZyme_family_groups.tsv trees/supermatrix_ML_rooted_94taxa.nwk "$OUT_DIR" <<'RSCRIPT'
suppressPackageStartupMessages({ library(ape); library(phytools) })

a <- commandArgs(trailingOnly = TRUE)
counts_file <- a[1]; stats_file <- a[2]; genomes_file <- a[3]; groups_file <- a[4]; tree_file <- a[5]; out_dir <- a[6]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_NONZERO <- 4      # families must have non-zero counts in at least four genomes
N_SIM_K     <- 999    # simulations for the test of Blomberg's K

# ----------------------------------------------------------------------------
# 1. Data: PCWDE family counts of the 94 genomes in the tree
# ----------------------------------------------------------------------------
tree <- read.tree(tree_file)
tree$tip.label <- gsub("-", "_", tree$tip.label)
stopifnot(is.rooted(tree))

counts <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
counts$Name <- gsub("-", "_", counts$Name)
pcwde <- read.delim(groups_file, stringsAsFactors = FALSE)$family
families <- intersect(sort(unique(pcwde)), colnames(counts))   # GH61 (the former name of AA9) has no column

X <- as.matrix(counts[, families])
storage.mode(X) <- "numeric"
rownames(X) <- counts$Name

species <- intersect(tree$tip.label, rownames(X))   # the counts table also lists two conspecific strains not in the tree
stopifnot(length(species) == Ntip(tree))
X <- X[species, , drop = FALSE]

# ----------------------------------------------------------------------------
# 2. Families with non-zero counts in at least four genomes; log10(x + 1)
# ----------------------------------------------------------------------------
keep <- colSums(X > 0) >= MIN_NONZERO & apply(X, 2, sd) > 0
X <- log10(X[, keep, drop = FALSE] + 1)
cat(sprintf("pPCA input: %d genomes x %d PCWDE families (of %d listed)\n",
            nrow(X), ncol(X), length(unique(pcwde))))

# ----------------------------------------------------------------------------
# 3. Phylogenetic signal of each family: Pagel's lambda and Blomberg's K
# ----------------------------------------------------------------------------
signal <- do.call(rbind, lapply(colnames(X), function(f) {
  x  <- setNames(X[, f], rownames(X))
  pl <- phylosig(tree, x, method = "lambda", test = TRUE)
  pk <- phylosig(tree, x, method = "K", test = TRUE, nsim = N_SIM_K)
  data.frame(family = f, lambda = pl$lambda, p_lambda = pl$P, K = pk$K, p_K = pk$P)
}))
write.csv(signal, file.path(out_dir, "phylogenetic_signal_per_family.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 4. pPCA on the correlation matrix under a lambda model
# ----------------------------------------------------------------------------
ppca <- phyl.pca(tree, X, method = "lambda", mode = "corr")
eig  <- diag(ppca$Eval)
var_explained <- eig / sum(eig)

scores   <- as.matrix(ppca$S);  colnames(scores)   <- paste0("PC", seq_len(ncol(scores)))
loadings <- as.matrix(ppca$L);  colnames(loadings) <- paste0("PC", seq_len(ncol(loadings)))
write.csv(data.frame(Name = rownames(scores), scores, check.names = FALSE),
          file.path(out_dir, "pPCA_scores.csv"), row.names = FALSE)
write.csv(data.frame(family = rownames(loadings), loadings, check.names = FALSE),
          file.path(out_dir, "pPCA_loadings.csv"), row.names = FALSE)
write.csv(data.frame(PC = seq_along(eig), eigenvalue = eig, var_explained = var_explained,
                     cumulative_var = cumsum(var_explained)),
          file.path(out_dir, "pPCA_variance_explained.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 5. Phylogenetic signal of the scores on the first two axes
# ----------------------------------------------------------------------------
score_signal <- do.call(rbind, lapply(c("PC1", "PC2"), function(pc) {
  pl <- phylosig(tree, scores[, pc], method = "lambda", test = TRUE)
  data.frame(axis = pc, lambda = pl$lambda, p_lambda = pl$P)
}))
write.csv(score_signal, file.path(out_dir, "phylogenetic_signal_of_scores.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 6. Spearman: BUSCO completeness vs pectin-degrading and total PCWDE genes
#    in the 51 lichenized genomes
# ----------------------------------------------------------------------------
st <- read.csv(stats_file, check.names = FALSE, stringsAsFactors = FALSE)
lichens <- with(read.delim(genomes_file, stringsAsFactors = FALSE), Name[Lifestyle == "lichen"])
st <- st[gsub("-", "_", st$Name) %in% intersect(species, gsub("-", "_", lichens)), ]
stopifnot(nrow(st) == 51)
spearman <- do.call(rbind, lapply(list(pectin = st$pec, PCWDE = st$cell + st$pec + st$lign), function(y) {
  ct <- suppressWarnings(cor.test(st[["C(%)"]], y, method = "spearman", exact = FALSE))
  data.frame(n = nrow(st), rho = unname(ct$estimate), p = ct$p.value)
}))
spearman <- cbind(genes = rownames(spearman), spearman)
write.table(spearman, file.path(out_dir, "spearman_BUSCO.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
summary_lines <- c(
  "Phylogenetic PCA of PCWDE families (phytools::phyl.pca, method = 'lambda', mode = 'corr')",
  sprintf("genomes: %d; families: %d (non-zero in >= %d genomes; log10(x + 1))", nrow(X), ncol(X), MIN_NONZERO),
  sprintf("Pagel's lambda of the trait covariance: %.4f", ppca$lambda),
  sprintf("variance explained: PC1 %.2f%%, PC2 %.2f%%", 100 * var_explained[1], 100 * var_explained[2]),
  sprintf("families with significant Pagel's lambda (p < 0.05): %d of %d", sum(signal$p_lambda < 0.05), nrow(signal)),
  sprintf("Pagel's lambda of the scores: PC1 %.3f, PC2 %.3f", score_signal$lambda[1], score_signal$lambda[2]),
  sprintf("Spearman, BUSCO vs %s genes (n = %d): rho = %.3f, p = %.3f",
          spearman$genes, spearman$n, spearman$rho, spearman$p))
writeLines(summary_lines, file.path(out_dir, "pPCA_summary.txt"))
cat(summary_lines, sep = "\n")
RSCRIPT
