#!/usr/bin/env bash
# =============================================================================
# 4_ultrametric_transformation.sh
# Materials and Methods: "Ultrametric transformation of the phylogeny"
#
# Converts the rooted maximum-likelihood supermatrix tree into a chronogram by
# penalized likelihood in treePL (input of the CAFE analyses, script 5).
#
#   Calibrations, both fixed at the published median ages:
#     root (most recent common ancestor of all 94 taxa)  458 Ma    (95% HPD 400-583;   Beimforde et al. 2014)
#     crown Lecanoromycetes (the 37 sampled genomes)     248.6 Ma  (95% HPD 199.7-303.0; Nelsen et al. 2020a)
#   Smoothing chosen by cross-validation among 100,000 ... 0.001 (factor 10).
#   57,900 alignment sites (the supermatrix of script 3).
#
#   Sensitivity chronograms: smoothing 0.1 and 10; both calibrations as minimum
#   and maximum ages spanning their 95% HPDs; both nodes fixed at the lower
#   (400 and 199.7 Ma) or upper (583 and 303.0 Ma) ends of the HPDs.
#
# Usage (from the repository root):
#   bash scripts/4_ultrametric_transformation.sh [TREE] [OUT_DIR]
#     TREE     default trees/supermatrix_ML_rooted_94taxa.nwk (script 3)
#     OUT_DIR  default results/4_chronogram
#
# Outputs (OUT_DIR):
#   chronogram.tre              final chronogram (trees/chronogram_treePL.nwk)
#   replicate_seed7.tre, replicate_seed99.tre   the same analysis with two other seeds (identical)
#   sensitivity/*.tre           sensitivity chronograms (trees/chronogram_sensitivity/)
#   cv_scores.txt, *.cfg, *.log treePL input and output
#
# Software: treePL (https://github.com/blackrim/treePL, commit f41af04, built with
# the bundled ADOL-C 2.6.3 and NLopt 2.4.2; set TREEPL to the executable); R with ape 5.8-1
# (environment: envs/lichen_ppca.yml).
# =============================================================================
set -euo pipefail

TREE=${1:-trees/supermatrix_ML_rooted_94taxa.nwk}
OUT_DIR=${2:-results/4_chronogram}
GENOMES=$(pwd)/${GENOMES:-data/genomes.tsv}
TREEPL=${TREEPL:-treePL}
RSCRIPT=${RSCRIPT:-Rscript}  # an R with ape (envs/lichen_ppca.yml)
command -v "$TREEPL" >/dev/null || { echo "ERROR: treePL not found (set TREEPL=/path/to/treePL)"; exit 1; }

NUMSITES=57900
SEED=20261002
mkdir -p "$OUT_DIR/sensitivity"
TREE=$(cd "$(dirname "$TREE")" && pwd)/$(basename "$TREE")
cd "$OUT_DIR"

# -----------------------------------------------------------------------------
# 1. Input tree and calibration nodes
#    Rooting on Sordariomyceta leaves the root branch on one side only; its
#    length is split equally between the two root edges. Each calibrated node
#    is given as two tips whose most recent common ancestor is that node.
# -----------------------------------------------------------------------------
"$RSCRIPT" --vanilla - "$TREE" "$GENOMES" > calibration_nodes.txt <<'RSCRIPT'
suppressPackageStartupMessages(library(ape))
a <- commandArgs(trailingOnly = TRUE)
tree <- read.tree(a[1]); tree$tip.label <- gsub("-", "_", tree$tip.label)
root <- Ntip(tree) + 1
re <- which(tree$edge[, 1] == root)
tree$edge.length[re] <- sum(tree$edge.length[re]) / 2
tree$node.label <- NULL
write.tree(tree, "ml94_rootsplit.nwk")
meta <- read.delim(a[2], stringsAsFactors = FALSE); meta$Name <- gsub("-", "_", meta$Name)
two_tips <- function(node) {            # one tip from each of the two descendant lineages
  kids <- tree$edge[tree$edge[, 1] == node, 2]
  sapply(kids[1:2], function(k) if (k <= Ntip(tree)) tree$tip.label[k] else extract.clade(tree, k)$tip.label[1])
}
lecano <- getMRCA(tree, intersect(meta$Name[meta$Class == "Lecanoromycetes"], tree$tip.label))
r2 <- two_tips(root); l2 <- two_tips(lecano)
stopifnot(getMRCA(tree, r2) == root, getMRCA(tree, l2) == lecano)
cat(sprintf("ROOT %s %s\nLECANO %s %s\n", r2[1], r2[2], l2[1], l2[2]))
RSCRIPT
read -r _ ROOT_A ROOT_B     < <(grep '^ROOT '   calibration_nodes.txt)
read -r _ LECANO_A LECANO_B < <(grep '^LECANO ' calibration_nodes.txt)

# treePL configuration: config NAME SMOOTH ROOT_MIN ROOT_MAX LECANO_MIN LECANO_MAX [extra lines...]
config () {
  local name=$1 smooth=$2 rmin=$3 rmax=$4 lmin=$5 lmax=$6; shift 6
  {
    echo "treefile = ml94_rootsplit.nwk"
    echo "smooth = $smooth"
    echo "numsites = $NUMSITES"
    echo "mrca = ROOT $ROOT_A $ROOT_B"
    echo "min = ROOT $rmin"
    echo "max = ROOT $rmax"
    echo "mrca = LECANO $LECANO_A $LECANO_B"
    echo "min = LECANO $lmin"
    echo "max = LECANO $lmax"
    echo "outfile = $name.tre"
    echo "thorough"
    printf '%s\n' "$@"
  } > "$name.cfg"
}
MEDIANS="458 458 248.6 248.6"

# -----------------------------------------------------------------------------
# 2. Priming: optimizer settings recommended by treePL
# -----------------------------------------------------------------------------
config prime 100 $MEDIANS prime "nthreads = 8" "seed = $SEED"
"$TREEPL" prime.cfg > prime.log 2>&1 < /dev/null
sed -n '/PLACE THE LINES BELOW/,$p' prime.log | tail -n +2 > prime_settings.txt
[[ -s prime_settings.txt ]] || { echo "ERROR: treePL prime gave no settings"; exit 1; }
OPT=(); while IFS= read -r line; do OPT+=("$line"); done < prime_settings.txt

# -----------------------------------------------------------------------------
# 3. Cross-validation of the smoothing value
# -----------------------------------------------------------------------------
config cv 100 $MEDIANS "${OPT[@]}" cv "cvstart = 100000" "cvstop = 0.001" "cvmultstep = 0.1" \
       "cvoutfile = cv_scores.txt" "nthreads = 16" "seed = $SEED"
"$TREEPL" cv.cfg > cv.log 2>&1 < /dev/null
SMOOTH=$(sed -n 's/^chisq: (\(.*\)) \(.*\)$/\2 \1/p' cv_scores.txt | sort -g | head -n 1 | cut -d' ' -f2)
echo "smoothing chosen by cross-validation: $SMOOTH" | tee smoothing_chosen.txt

# -----------------------------------------------------------------------------
# 4. Final chronogram (seed 20261002) and two replicate runs
# -----------------------------------------------------------------------------
for seed in "$SEED" 7 99; do
  name=chronogram; [[ $seed == "$SEED" ]] || name=replicate_seed$seed
  config "$name" "$SMOOTH" $MEDIANS "${OPT[@]}" "nthreads = 8" "seed = $seed"
  "$TREEPL" "$name.cfg" > "$name.log" 2>&1 < /dev/null
done

# -----------------------------------------------------------------------------
# 5. Sensitivity chronograms
# -----------------------------------------------------------------------------
sensitivity () {   # NAME SMOOTH ROOT_MIN ROOT_MAX LECANO_MIN LECANO_MAX
  config "sensitivity/$1" "$2" "$3" "$4" "$5" "$6" "${OPT[@]}" "nthreads = 4" "seed = $SEED"
  "$TREEPL" "sensitivity/$1.cfg" > "sensitivity/$1.log" 2>&1 < /dev/null
}
sensitivity smoothing_0.1 0.1       $MEDIANS
sensitivity smoothing_10  10        $MEDIANS
sensitivity HPD_bounds    "$SMOOTH" 400 583 199.7 303.0    # calibrations as minimum and maximum ages
sensitivity HPD_lower     "$SMOOTH" 400 400 199.7 199.7    # both nodes at the lower ends of the HPDs
sensitivity HPD_upper     "$SMOOTH" 583 583 303.0 303.0    # both nodes at the upper ends of the HPDs

echo "chronogram: $OUT_DIR/chronogram.tre; sensitivity chronograms: $OUT_DIR/sensitivity/"
