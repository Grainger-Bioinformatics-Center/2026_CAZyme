#!/usr/bin/env bash
# =============================================================================
# 1_metagenome_decontamination.sh
# Materials and Methods: "Metagenome decontamination"
#
# Removes contigs of contaminating organisms from genome assemblies before
# annotation. Every contig is classified with Kraken 2 against a database built
# from NCBI nt and the NCBI taxonomy; contigs assigned to Bacteria (taxid 2),
# Archaea (2157), Viruses (10239) or Viridiplantae (33090; green-algal photobiont
# and plant sequences), or to any of their descendant taxa, are removed with
# SeqKit. Unclassified contigs and contigs assigned to any other taxon
# (fungal, other eukaryotic, root) are kept.
#
# Usage:
#   KRAKEN_DB=/path/to/kraken2_nt_db  bash scripts/1_metagenome_decontamination.sh  ASSEMBLY_DIR  [OUT_DIR]
#
#   ASSEMBLY_DIR  assemblies to screen, one FASTA per genome (<genome>.fna, .fa or .fasta)
#   OUT_DIR       default results/1_decontamination
#
# Outputs (OUT_DIR):
#   <genome>.fasta                    filtered assembly (input to script 2)
#   kraken2/<genome>.kraken.out       per-contig Kraken 2 classification
#   kraken2/<genome>.kraken.report    Kraken 2 report
#   decontamination_summary.tsv       contigs and base pairs before and after filtering
#
# Settings: Kraken 2 --confidence 0.05. A database larger than the available
# memory needs --memory-mapping (set MEMORY_MAPPING=true, the default).
# Software: Kraken 2 v2.1.3, SeqKit v2.13.0 (environment: envs/kraken.yml).
#
# Note on the published data set (Table S2): 54 assemblies were filtered with
# this rule. Two earlier annotation batches used other rules: 18 assemblies were
# classified at the default confidence (0) and only contigs assigned exactly to
# taxid 2, 2157, 2759, 33208, 33090, 4773, 4751 or 10239 were removed; for five
# assemblies only unclassified contigs and contigs assigned to Fungi outside
# Ascomycota, to Eukaryota or to the root were kept. Seventeen assemblies (16
# culture-derived and the metagenome-assembled Alyxoria varia) were not screened.
# Re-applying this script's rule to the Kraken 2 classifications of all 77
# screened assemblies flagged 111 further contigs (70 kb), none of which carries
# a CAZyme gene, so no gene count in the article changes.
# =============================================================================
set -euo pipefail

ASSEMBLY_DIR=${1:?usage: KRAKEN_DB=... bash $0 ASSEMBLY_DIR [OUT_DIR]}
OUT_DIR=${2:-results/1_decontamination}
KRAKEN_DB=${KRAKEN_DB:?set KRAKEN_DB to the Kraken 2 database built from NCBI nt}
THREADS=${THREADS:-30}
CONFIDENCE=0.05
MEMORY_MAPPING=${MEMORY_MAPPING:-true}

for tool in kraken2 seqkit; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not found (conda env create -f envs/kraken.yml)"; exit 1; }
done
NODES="$KRAKEN_DB/taxonomy/nodes.dmp"
[[ -s "$NODES" ]] || { echo "ERROR: $NODES not found"; exit 1; }

MM_OPTION=""; [[ "$MEMORY_MAPPING" == true ]] && MM_OPTION=--memory-mapping
KDIR="$OUT_DIR/kraken2"
mkdir -p "$KDIR"

# -----------------------------------------------------------------------------
# 1. Taxids to remove: Bacteria, Archaea, Viruses, Viridiplantae and all their
#    descendants, from the taxonomy of the Kraken 2 database
# -----------------------------------------------------------------------------
REMOVE_TAXIDS="$KDIR/removed_taxids.txt"
awk -F'|' '
  { gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2); parent[$1] = $2 }
  END {
    split("2 2157 10239 33090", top, " ")
    for (i in top) remove[top[i]] = 1
    for (t in parent) {
      cur = t
      while (cur != 1 && (cur in parent)) {
        if (cur in remove) { print t; break }
        cur = parent[cur]
      }
    }
  }' "$NODES" | sort -u > "$REMOVE_TAXIDS"
echo "taxids removed (Bacteria, Archaea, Viruses, Viridiplantae and descendants): $(wc -l < "$REMOVE_TAXIDS" | tr -d ' ')"

# -----------------------------------------------------------------------------
# 2. Classify and filter each assembly
# -----------------------------------------------------------------------------
SUMMARY="$OUT_DIR/decontamination_summary.tsv"
printf "genome\tcontigs_in\tbp_in\tcontigs_out\tbp_out\tpct_bp_retained\n" > "$SUMMARY"

shopt -s nullglob
for fasta in "$ASSEMBLY_DIR"/*.{fna,fa,fasta}; do
  genome=$(basename "${fasta%.*}")
  out="$KDIR/$genome.kraken.out"

  kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --confidence "$CONFIDENCE" $MM_OPTION \
          --report "$KDIR/$genome.kraken.report" --output "$out" "$fasta"

  # keep unclassified (U) contigs and classified (C) contigs outside the removed taxa
  awk -v REMOVE="$REMOVE_TAXIDS" '
    BEGIN { while ((getline t < REMOVE) > 0) rm[t] = 1 }
    $1 == "U"                     { print $2; next }
    $1 == "C" && !($3 in rm)      { print $2 }' "$out" | sort -u > "$KDIR/$genome.keep.ids"

  seqkit grep -f "$KDIR/$genome.keep.ids" "$fasta" > "$OUT_DIR/$genome.fasta"

  read -r n_in bp_in   < <(seqkit stats -T "$fasta"                | awk 'NR == 2 { print $4, $5 }')
  read -r n_out bp_out < <(seqkit stats -T "$OUT_DIR/$genome.fasta" | awk 'NR == 2 { print $4, $5 }')
  awk -v g="$genome" -v a="$n_in" -v b="$bp_in" -v c="$n_out" -v d="$bp_out" \
      'BEGIN { printf "%s\t%s\t%s\t%s\t%s\t%.2f\n", g, a, b, c, d, (b > 0 ? 100 * d / b : 0) }' >> "$SUMMARY"
  echo "$genome: kept $n_out of $n_in contigs"
done

echo "Filtered assemblies and summary: $OUT_DIR"
