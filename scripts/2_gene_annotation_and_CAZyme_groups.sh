#!/usr/bin/env bash
# =============================================================================
# 2_gene_annotation_and_CAZyme_groups.sh
# Materials and Methods: "Gene identification and annotation" and
#                        "CAZyme groups identification and selection"
#
# Re-annotates every genome with the same funannotate pipeline, assesses the
# completeness of the predicted proteomes with BUSCO, counts the CAZyme genes
# per family from the dbCAN annotations and sums the families acting on
# (hemi-)cellulose, pectin and lignin (data/CAZyme_family_groups.tsv).
#
# Usage:
#   FUNANNOTATE_DB=/path/to/funannotate_db  [BUSCO=/path/to/busco]  bash scripts/2_gene_annotation_and_CAZyme_groups.sh  ASSEMBLY_DIR  [OUT_DIR]
#
#   ASSEMBLY_DIR  one FASTA per genome, named <genome>.fasta (the Name column of
#                 data/genomes.tsv); metagenome-derived assemblies after script 1
#   OUT_DIR       default results/2_annotation
#
# Outputs (OUT_DIR):
#   funannotate/<genome>/           funannotate predict and annotate results
#   busco/<genome>/                 BUSCO run on the predicted proteome
#   CAZyme_family_counts.csv        genes per CAZyme family and genome (Table S4; data/CAZyme_family_counts.csv)
#   genome_statistics.csv           BUSCO, assembly statistics, gene and CAZyme totals,
#                                   (hemi-)cellulose (cell), pectin (pec) and lignin (lign)
#                                   genes, sugar transporters (Pfam PF00083) (Table S1;
#                                   data/genome_statistics.csv)
#
# Software: funannotate v1.8.17 with AUGUSTUS v3.5.0, GeneMark-ES v4.71, SNAP 2013-11-29,
# GlimmerHMM v3.0.4, EVidenceModeler v1.1.1, tRNAscan-SE v2.0.12, tantan v49;
# annotation databases dbCAN 12.0, Pfam 36.0, MEROPS 12.0, eggNOG-mapper v2.1.12
# (environment: envs/funannotate.yml); BUSCO v6.0.0 in protein mode with
# ascomycota_odb10 (environment: envs/phy.yml).
#
# Notes on the published data set:
#  - Contigs were renamed (ctg prefix) before prediction because AUGUSTUS accepts only
#    short contig names. For the first batch of genomes (the lichens), the contigs
#    were renamed in input order without sorting; all other genomes used
#    "funannotate sort" as below.
#  - Epichloe typhina: contigs shorter than 500 bp were removed (funannotate sort
#    --minlen 500) to drop a 289-bp contig that stopped gene prediction.
#  - The first batch was annotated with the default BUSCO set of "funannotate annotate".
# =============================================================================
set -euo pipefail

ASSEMBLY_DIR=${1:?usage: FUNANNOTATE_DB=... bash $0 ASSEMBLY_DIR [OUT_DIR]}
OUT_DIR=${2:-results/2_annotation}
export FUNANNOTATE_DB=${FUNANNOTATE_DB:?set FUNANNOTATE_DB to the funannotate database folder}
THREADS=${THREADS:-30}
GENOMES=${GENOMES:-data/genomes.tsv}
GROUPS=${GROUPS:-data/CAZyme_family_groups.tsv}
BUSCO_LINEAGE=ascomycota_odb10
BUSCO=${BUSCO:-busco}       # BUSCO v6.0.0 (envs/phy.yml; the funannotate environment has an older BUSCO)

for tool in funannotate "$BUSCO"; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not found"; exit 1; }
done
FUN="$OUT_DIR/funannotate"; BUS="$OUT_DIR/busco"
mkdir -p "$FUN" "$BUS"

# -----------------------------------------------------------------------------
# 1. Gene prediction, functional annotation and BUSCO, one genome at a time
# -----------------------------------------------------------------------------
for genome in $(tail -n +2 "$GENOMES" | cut -f1); do      # genome names contain no spaces
  fasta="$ASSEMBLY_DIR/$genome.fasta"
  [[ -s "$fasta" ]] || { echo "WARNING: $fasta not found, skipped"; continue; }
  out="$FUN/$genome"
  species=$(awk -F_ '{ print $1 " " $2 }' <<< "$genome")                     # e.g. "Pyrenula aspistea"
  locus=$(awk -F_ '{ print $1 $2 }' <<< "$genome" | tr -cd '[:alnum:]' | cut -c1-24)
  minlen=0; [[ "$genome" == Epichloe_typhina_* ]] && minlen=500
  mkdir -p "$out"

  # short contig names for AUGUSTUS
  funannotate sort -i "$fasta" -o "$out/$genome.sorted.fa" -b ctg --minlen "$minlen"

  # public assemblies are used with their soft-masking; an unmasked assembly
  # (no lower-case bases; here Fitzroyomyces cyperi) is soft-masked with tantan
  if awk '!/^>/ && /[acgt]/ { found = 1; exit } END { exit !found }' "$out/$genome.sorted.fa"; then
    cp "$out/$genome.sorted.fa" "$out/$genome.masked.fa"
  else
    funannotate mask -i "$out/$genome.sorted.fa" -o "$out/$genome.masked.fa" --method tantan --cpus "$THREADS"
  fi

  funannotate predict  --input "$out/$genome.masked.fa" --out "$out" --name "$locus" --species "$species" --cpus "$THREADS"
  funannotate annotate --input "$out" --cpus "$THREADS" --busco_db ascomycota

  proteins=$(ls "$out"/predict_results/*.proteins.fa | head -n1)
  "$BUSCO" -i "$proteins" -l "$BUSCO_LINEAGE" -m proteins -o "$genome" --out_path "$BUS" -c "$THREADS" --force
done

# -----------------------------------------------------------------------------
# 2. CAZyme genes per family and genome-level statistics
#    A gene is counted once for every CAZyme family it is annotated with
#    ("CAZy:<family>" in annotate_misc/annotations.dbCAN.txt), so a gene with
#    domains of two families counts for both. Subfamily suffixes are kept
#    as written by funannotate.
# -----------------------------------------------------------------------------
python3 - "$GENOMES" "$GROUPS" "$FUN" "$BUS" "$OUT_DIR" <<'PY'
import csv, glob, gzip, os, re, sys

genomes_tsv, groups_tsv, fun_dir, busco_dir, out_dir = sys.argv[1:6]
genomes = [r["Name"] for r in csv.DictReader(open(genomes_tsv), delimiter="\t")]
groups = {}
for r in csv.DictReader(open(groups_tsv), delimiter="\t"):
    groups.setdefault(r["substrate"], set()).add(r["family"])

FAMILY = re.compile(r"^[A-Za-z]+[0-9]+(_[A-Za-z0-9_]+)?")

def first(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None

def cazy_counts(path):
    """genes per CAZyme family; one count per (gene, family) pair"""
    pairs = set()
    for line in open(path):
        gene = line.split("\t", 1)[0]
        for part in line.split("CAZy:")[1:]:
            m = FAMILY.match(part)
            if m:
                pairs.add((gene, m.group(0)))
    counts = {}
    for _, fam in pairs:
        counts[fam] = counts.get(fam, 0) + 1
    return counts

def busco_scores(path):
    """C, D, F, M as fractions from a BUSCO short summary"""
    for line in open(path):
        if "C:" in line and "[S:" in line and "F:" in line and "M:" in line:
            v = {k: float(re.search(k + r":([0-9.]+)%", line).group(1)) / 100 for k in ("C", "D", "F", "M")}
            return [f"{v[k]:.3f}" for k in ("C", "D", "F", "M")]
    return ["", "", "", ""]

def scaffold_stats(path):
    lengths, n = [], 0
    for line in open(path):
        if line.startswith(">"):
            if n: lengths.append(n)
            n = 0
        else:
            n += len(line.strip())
    if n: lengths.append(n)
    lengths.sort(reverse=True)
    half, acc, n50 = sum(lengths) / 2, 0, 0
    for x in lengths:
        acc += x
        if acc >= half:
            n50 = x
            break
    return len(lengths), n50

def fmt_n50(n):
    return f"{n / 1e6:.0f} MB" if n >= 1e6 else (f"{n / 1e3:.0f} KB" if n >= 1e3 else f"{n} bp")

per_genome, stats = {}, {}
for g in genomes:
    d = os.path.join(fun_dir, g)
    dbcan = os.path.join(d, "annotate_misc", "annotations.dbCAN.txt")
    if not os.path.exists(dbcan):
        print(f"WARNING: no dbCAN annotation for {g}", file=sys.stderr)
        continue
    counts = cazy_counts(dbcan)
    per_genome[g] = counts
    pfam = os.path.join(d, "annotate_misc", "annotations.pfam.txt")
    gff = first(os.path.join(d, "annotate_results", "*.gff3")) or first(os.path.join(d, "predict_results", "*.gff3"))
    scaf = first(os.path.join(d, "annotate_results", "*.scaffolds.fa")) or first(os.path.join(d, "predict_results", "*.scaffolds.fa"))
    busco = first(os.path.join(busco_dir, g, "**", "short_summary*.txt")) or first(os.path.join(busco_dir, g, "short_summary*.txt"))
    n_scaf, n50 = scaffold_stats(scaf) if scaf else ("", None)
    n_cds = sum(1 for line in open(gff) if line.split("\t")[2:3] == ["mRNA"])
    pf00083 = sum(1 for line in open(pfam) if "PF00083" in line) if os.path.exists(pfam) else 0
    by_class = {c: sum(v for f, v in counts.items() if f.startswith(c)) for c in ("AA", "CBM", "CE", "GH", "GT", "PL")}
    by_group = {k: sum(counts.get(f, 0) for f in fams) for k, fams in groups.items()}
    stats[g] = (busco_scores(busco) if busco else ["", "", "", ""]) + [fmt_n50(n50) if n50 else "", n_scaf, n_cds] + \
               [by_class[c] for c in ("AA", "CBM", "CE", "GH", "GT", "PL")] + \
               [by_group["(hemi-)cellulose"], by_group["pectin"], by_group["lignin"], pf00083]

# families in CAZy order (AA, CBM, CE, GH, GT, PL), then by number
def fam_key(f):
    order = ["AA", "CBM", "CE", "GH", "GT", "PL"]
    prefix = re.match(r"[A-Za-z]+", f).group(0)
    num = int(re.search(r"[0-9]+", f).group(0))
    return (order.index(prefix) if prefix in order else len(order), num, f)

families = sorted({f for c in per_genome.values() for f in c}, key=fam_key)
with open(os.path.join(out_dir, "CAZyme_family_counts.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Name"] + families)
    for g in per_genome:
        w.writerow([g] + [per_genome[g].get(f, 0) for f in families])
with open(os.path.join(out_dir, "genome_statistics.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Name", "C(%)", "D(%)", "F(%)", "M(%)", "Scaffold N50", "Scaffold Count", "Number of CDS",
                "CAZymes_AA", "CAZymes_CBM", "CAZymes_CE", "CAZymes_GH", "CAZymes_GT", "CAZymes_PL",
                "cell", "pec", "lign", "Pfam_PF00083"])
    for g in stats:
        w.writerow([g] + stats[g])
print(f"{len(per_genome)} genomes, {len(families)} CAZyme families -> {out_dir}")
PY
