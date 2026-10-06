#!/usr/bin/env bash
# =============================================================================
# 3_phylogenetic_reconstruction.sh
# Materials and Methods: "Phylogenetic reconstruction"
#
# A. Supermatrix tree. Single-copy orthogroups of the 96 predicted proteomes
#    (OrthoFinder) are aligned (MAFFT), trimmed (trimAl -automated1) and
#    concatenated (AMAS); the maximum-likelihood tree is estimated with RAxML-NG
#    under a per-gene-partitioned LG+G+F model with 100 bootstrap replicates.
#    The tree is rooted on Sordariomyceta (all Sordariomycetes and Leotiomycetes)
#    and the second strains of two species are removed (96 -> 94 tips).
# B. UnFATE tree. The 195 single-copy UnFATE markers (Ametrano et al. 2025) are
#    retrieved from each genome assembly with miniprot, aligned, trimmed,
#    concatenated and analysed as in A.
# C. Robinson-Foulds distance between the two trees.
#
# Usage (from the repository root):
#   bash scripts/3_phylogenetic_reconstruction.sh  ASSEMBLY_DIR  UNFATE_MARKERS  [ANNOTATION_DIR]  [OUT_DIR]
#
#   ASSEMBLY_DIR    genome assemblies, <genome>.fasta (the input of script 2)
#   UNFATE_MARKERS  UnFATE_markers_195.fas from https://github.com/claudioametrano/UnFATE
#   ANNOTATION_DIR  output of script 2 (default results/2_annotation); the predicted
#                   proteomes are read from funannotate/<genome>/predict_results/
#   OUT_DIR         default results/3_phylogeny
#   The 96 genomes are those of data/genomes.tsv (the 94 analysed genomes and the
#   second strains of two species).
#
# Outputs (OUT_DIR):
#   supermatrix/supermatrix_ML.raxml.support    96-taxon tree with bootstrap support
#                                               (trees/supermatrix_ML_96taxa.raxml.support)
#   supermatrix/supermatrix_ML_rooted_94taxa.nwk rooted, pruned tree used by scripts 4-6
#                                               (trees/supermatrix_ML_rooted_94taxa.nwk)
#   UnFATE/UnFATE_195markers.raxml.support      UnFATE tree (trees/UnFATE_195markers.raxml.support)
#   tree_comparison.txt                         Robinson-Foulds distance
#
# Software: OrthoFinder v3.1.2, MAFFT v7.526, trimAl v1.5, AMAS v1.0, RAxML-NG v2.0.1,
# miniprot v0.18 (environment: envs/phy.yml); R with ape 5.8-1 and phangorn 2.12.1
# (environment: envs/lichen_ppca.yml; set RSCRIPT if it is not the default Rscript).
#
# Notes on the published trees:
#  - The supermatrix has 109 of the 110 single-copy orthogroups: OG0002496, the
#    last line of Orthogroups_SingleCopyOrthologues.txt, was not read by the
#    original script (the file has no final newline). It is excluded below so that
#    the published supermatrix (96 taxa x 57,900 sites) is reproduced.
#  - RAxML-NG seeds: the published runs drew the seeds 1782879946 (supermatrix)
#    and 1782910460 (UnFATE), as recorded in their logs; they are fixed here.
# =============================================================================
set -euo pipefail

ASSEMBLY_DIR=${1:?usage: bash $0 ASSEMBLY_DIR UNFATE_MARKERS [ANNOTATION_DIR] [OUT_DIR]}
UNFATE_MARKERS=${2:?usage: bash $0 ASSEMBLY_DIR UNFATE_MARKERS [ANNOTATION_DIR] [OUT_DIR]}
ANNOTATION_DIR=${3:-results/2_annotation}
OUT_DIR=${4:-results/3_phylogeny}
GENOMES=${GENOMES:-data/genomes.tsv}
THREADS=${THREADS:-48}
RSCRIPT=${RSCRIPT:-Rscript}  # an R with ape and phangorn (envs/lichen_ppca.yml)

MODEL="LG+G+F"
BOOTSTRAPS=100
OUTGROUP="Sordaria_macrospora_GCA_000182805.2,Podospora_anserina_GCA_000226545.1"   # RAxML-NG output orientation only
SORDARIOMYCETA=data/outgroup_Sordariomyceta.txt                                      # rooting clade
SECOND_STRAINS="Aureobasidium_pullulans_GCA_003336255.1,Trypethelium_eluteriae_GCA_051942245.1"
EXCLUDED_OG=OG0002496
SEED_SUPERMATRIX=1782879946
SEED_UNFATE=1782910460
MIN_IDENTITY=0.40           # UnFATE: minimum miniprot identity of a marker hit
MIN_COVERAGE=0.30           # UnFATE: minimum coverage of the marker protein
MIN_OCCUPANCY=0.60          # UnFATE: markers found in at least 60% of the genomes

for tool in orthofinder mafft trimal raxml-ng miniprot python3 "$RSCRIPT"; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not found (conda env create -f envs/phy.yml)"; exit 1; }
done
if command -v AMAS.py >/dev/null; then AMAS=AMAS.py; else AMAS=amas; fi

SM="$OUT_DIR/supermatrix"; UF="$OUT_DIR/UnFATE"
mkdir -p "$SM"/{proteomes,og,aligned,trimmed} "$UF"/{miniprot,proteins,markers,aligned,trimmed}

# AMAS partition file -> RAxML-NG partition file with one LG+G+F model per gene
raxml_partitions () {   # $1 AMAS partitions, $2 output
  python3 - "$1" "$2" "$MODEL" <<'PY'
import re, sys
src, out, model = sys.argv[1:4]
with open(src) as fh, open(out, "w") as oh:
    for line in fh:
        line = re.sub(r"^charset\s+", "", line.strip().rstrip(";"), flags=re.I)
        if "," in line and re.match(r"^\s*\S+\s*=", line.split(",", 1)[1]):
            line = line.split(",", 1)[1].strip()          # drop a leading model tag, e.g. "WAG, "
        m = re.match(r"^([\w\-.]+)\s*=\s*([\d\-,\s]+)$", line)
        if m:
            oh.write(f"{model}, {m.group(1)} = {re.sub(r'\s+', '', m.group(2))}\n")
PY
}

# =============================================================================
# A. Supermatrix tree
# =============================================================================
# A1. proteomes with the genome name in every header (">genome|protein")
for genome in $(tail -n +2 "$GENOMES" | cut -f1); do
  f=$(ls "$ANNOTATION_DIR/funannotate/$genome"/predict_results/*.proteins.fa | head -n 1)
  awk -v g="$genome" '/^>/ { id = substr($1, 2); print ">" g "|" id; next } { print }' "$f" > "$SM/proteomes/$genome.fa"
done

# A2. orthogroups
orthofinder -f "$SM/proteomes" -t "$THREADS" -a "$THREADS"
OF=$(ls -d "$SM"/proteomes/OrthoFinder/Results_* | tail -n 1)

# A3. single-copy orthogroups: genome names as headers, MAFFT, trimAl
for og in $(grep -v "^$EXCLUDED_OG\$" "$OF/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"); do
  sed -E 's/^>([^|]+)\|.*/>\1/' "$OF/Single_Copy_Orthologue_Sequences/$og.fa" > "$SM/og/$og.fa"
  mafft --auto --thread "$THREADS" "$SM/og/$og.fa" > "$SM/aligned/$og.fa" 2>/dev/null
  trimal -in "$SM/aligned/$og.fa" -out "$SM/trimmed/$og.fa" -automated1 >/dev/null
  [[ -s "$SM/trimmed/$og.fa" ]] || cp "$SM/aligned/$og.fa" "$SM/trimmed/$og.fa"   # keep the alignment if trimAl returns nothing
done
echo "single-copy orthogroups in the supermatrix: $(ls "$SM"/trimmed/*.fa | wc -l)"

# A4. supermatrix and partitions
$AMAS concat -i "$SM"/trimmed/*.fa -f fasta -d aa -u phylip -t "$SM/supermatrix.phy" -p "$SM/partitions.amas.txt"
raxml_partitions "$SM/partitions.amas.txt" "$SM/partitions.txt"

# A5. maximum-likelihood tree with 100 bootstrap replicates
raxml-ng --all --msa "$SM/supermatrix.phy" --model "$SM/partitions.txt" --bs-trees "$BOOTSTRAPS" \
         --outgroup "$OUTGROUP" --seed "$SEED_SUPERMATRIX" --threads "$THREADS" --prefix "$SM/supermatrix_ML"

# A6. root on Sordariomyceta and remove the second strains
"$RSCRIPT" --vanilla - "$SM/supermatrix_ML.raxml.support" "$SORDARIOMYCETA" "$SECOND_STRAINS" \
        "$SM/supermatrix_ML_rooted_94taxa.nwk" <<'RSCRIPT'
suppressPackageStartupMessages(library(ape))
a <- commandArgs(trailingOnly = TRUE)
tr <- read.tree(a[1])                                    # bootstrap support as node labels
og <- intersect(readLines(a[2]), tr$tip.label)
# root first outside Sordariomyceta, then on the most recent common ancestor of Sordariomyceta
tr <- root(tr, outgroup = setdiff(tr$tip.label, og)[1], resolve.root = TRUE)
m  <- getMRCA(tr, og)
stopifnot(length(extract.clade(tr, m)$tip.label) == length(og))   # Sordariomyceta is monophyletic
tr <- ladderize(root(tr, node = m, resolve.root = TRUE, edgelabel = TRUE))
tr <- drop.tip(tr, strsplit(a[3], ",")[[1]])
write.tree(tr, a[4])
cat(sprintf("rooted on Sordariomyceta (%d taxa); %d tips written to %s\n", length(og), Ntip(tr), a[4]))
RSCRIPT

# =============================================================================
# B. UnFATE tree
# =============================================================================
# B1. protein-to-genome alignment of the marker proteins
for genome in $(tail -n +2 "$GENOMES" | cut -f1); do
  fasta="$ASSEMBLY_DIR/$genome.fasta"
  miniprot -t "$THREADS" --gff "$fasta" "$UNFATE_MARKERS" > "$UF/miniprot/$genome.gff" 2> "$UF/miniprot/$genome.log"
done

# B2. best hit per marker and genome (identity >= 0.40, marker coverage >= 0.30),
#     spliced coding sequence translated to protein
python3 - "$ASSEMBLY_DIR" "$UF" "$UNFATE_MARKERS" "$MIN_IDENTITY" "$MIN_COVERAGE" <<'PY'
import glob, os, re, sys

assembly_dir, uf, markers, min_id, min_cov = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
CODON = dict(zip([a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
                 "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def marker_id(name):            # ">GCA_001680685.1-1627at4890" -> "1627at4890"
    return name.rsplit("-", 1)[-1]

def read_fasta(path):
    seqs, name, parts = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None: seqs[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        elif line:
            parts.append(line)
    if name is not None: seqs[name] = "".join(parts)
    return seqs

def translate(s):
    s = s.upper()
    return "".join("X" if "N" in s[i:i + 3] else CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))

# length of each marker = longest of its reference proteins
marker_len = {}
for name, seq in read_fasta(markers).items():
    m = marker_id(name)
    marker_len[m] = max(marker_len.get(m, 0), len(seq))

for gff in sorted(glob.glob(os.path.join(uf, "miniprot", "*.gff"))):
    genome = os.path.basename(gff)[:-4]
    mrna, cds = {}, {}
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            target = re.match(r"^(\S+)\s+(\d+)\s+(\d+)", attrs["Target"])
            mrna[attrs["ID"]] = dict(marker=marker_id(target.group(1)), q_start=int(target.group(2)),
                                     q_end=int(target.group(3)), identity=float(attrs.get("Identity", 0)),
                                     score=float(f[5]) if f[5] not in (".", "") else 0.0,
                                     contig=f[0], strand=f[6], id=attrs["ID"])
        elif f[2] == "CDS":
            cds.setdefault(attrs["Parent"], []).append((int(f[3]), int(f[4])))
    best = {}                                   # highest-scoring alignment per marker
    for m in mrna.values():
        if m["marker"] not in best or m["score"] > best[m["marker"]]["score"]:
            best[m["marker"]] = m
    if not best: continue
    genome_seqs = read_fasta(os.path.join(assembly_dir, genome + ".fasta"))
    os.makedirs(os.path.join(uf, "proteins", genome), exist_ok=True)
    for marker, m in best.items():
        coverage = (m["q_end"] - m["q_start"] + 1) / marker_len[marker]
        if m["identity"] < min_id or coverage < min_cov or m["contig"] not in genome_seqs: continue
        contig = genome_seqs[m["contig"]]
        seq = "".join(contig[s - 1:e] for s, e in sorted(cds.get(m["id"], [])))
        if m["strand"] == "-": seq = seq.translate(COMP)[::-1]
        protein = translate(seq).rstrip("*").replace("*", "X")
        if protein:
            with open(os.path.join(uf, "proteins", genome, marker + ".faa"), "w") as out:
                out.write(f">{genome}\n" + "\n".join(protein[i:i + 60] for i in range(0, len(protein), 60)) + "\n")
PY

# B3. markers found in at least 60% of the genomes (at least four)
n_genomes=$(find "$UF/proteins" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
min_genomes=$(python3 -c "import math; print(max(4, math.ceil($n_genomes * $MIN_OCCUPANCY)))")
for marker in $(find "$UF/proteins" -name '*.faa' -exec basename {} .faa \; | sort -u); do
  files=( "$UF"/proteins/*/"$marker.faa" )
  (( ${#files[@]} >= min_genomes )) || continue
  cat "${files[@]}" > "$UF/markers/$marker.faa"
  mafft --auto --thread "$THREADS" "$UF/markers/$marker.faa" > "$UF/aligned/$marker.faa" 2>/dev/null
  trimal -in "$UF/aligned/$marker.faa" -out "$UF/trimmed/$marker.faa" -automated1 >/dev/null 2>&1 || true
  [[ -s "$UF/trimmed/$marker.faa" ]] || cp "$UF/aligned/$marker.faa" "$UF/trimmed/$marker.faa"
done
echo "UnFATE markers in at least $min_genomes of $n_genomes genomes: $(ls "$UF"/trimmed/*.faa | wc -l)"

# B4. supermatrix, partitions and maximum-likelihood tree
$AMAS concat -i "$UF"/trimmed/*.faa -f fasta -d aa -u phylip -t "$UF/supermatrix.phy" -p "$UF/partitions.amas.txt"
raxml_partitions "$UF/partitions.amas.txt" "$UF/partitions.txt"
raxml-ng --all --msa "$UF/supermatrix.phy" --model "$UF/partitions.txt" --bs-trees "$BOOTSTRAPS" \
         --outgroup "$OUTGROUP" --seed "$SEED_UNFATE" --threads "$THREADS" --prefix "$UF/UnFATE_195markers"

# =============================================================================
# C. Robinson-Foulds distance between the supermatrix and UnFATE trees
# =============================================================================
"$RSCRIPT" --vanilla - "$SM/supermatrix_ML.raxml.support" "$UF/UnFATE_195markers.raxml.support" <<'RSCRIPT' | tee "$OUT_DIR/tree_comparison.txt"
suppressPackageStartupMessages({ library(ape); library(phangorn) })
a  <- commandArgs(trailingOnly = TRUE)
t1 <- read.tree(a[1]); t2 <- read.tree(a[2])
t1$tip.label <- gsub("-", "_", t1$tip.label); t2$tip.label <- gsub("-", "_", t2$tip.label)
common <- intersect(t1$tip.label, t2$tip.label)
t1 <- keep.tip(t1, common); t2 <- keep.tip(t2, common)
rf <- as.numeric(RF.dist(t1, t2))                        # on unrooted trees
max_rf <- 2 * (length(common) - 3)
cat(sprintf("taxa: %d\nRobinson-Foulds distance: %d\nnormalized: %.3f\nshared splits: %.1f%%\n",
            length(common), as.integer(rf), rf / max_rf, 100 * (1 - rf / max_rf)))
RSCRIPT
