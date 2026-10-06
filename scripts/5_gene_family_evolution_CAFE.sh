#!/usr/bin/env bash
# =============================================================================
# 5_gene_family_evolution_CAFE.sh
# Materials and Methods: "Gene-family evolution (CAFE)"
#
# Expansion and contraction of CAZyme gene families along the chronogram with
# CAFE 5: one genome-wide birth-death rate (lambda; -k 1), Poisson root-size
# distribution estimated from all families in the run (-p), family-wide and
# per-branch significance at p < 0.05.
#
#   1. CAFE input: genes per CAZyme family and genome (data/CAZyme_family_counts.csv);
#      families with 100 or more genes in a genome or with non-zero counts in
#      fewer than four genomes are excluded (177 families remain). In a second
#      input, the 11 pectin families are summed into one "PECTIN" family.
#   2. Per-family and aggregated-pectin runs on the final chronogram and on the
#      five sensitivity chronograms (script 4).
#   3. Other plant cell wall substrates and the root-size distribution (final chronogram):
#      - the three lignin families summed into one "LIG" family (run as in 2);
#      - the PECTIN and LIG aggregates analysed alone, with lambda fixed at the
#        genome-wide estimate and a Poisson root distribution estimated from the
#        aggregate itself (-p) or a uniform root distribution (no -p);
#      - all (hemi-)cellulose, pectin and lignin families with per-branch
#        p-values for every family (-P 1.0; lambda and root distribution of the
#        per-family run). (Hemi-)cellulose families cannot be analysed as an
#        aggregate (CAFE cannot compute the likelihood for up to 307 genes per genome),
#        so their reconstructed sizes are summed.
#   4. Summary tables.
#
# Usage (from the repository root):
#   bash scripts/5_gene_family_evolution_CAFE.sh [CHRONOGRAM] [SENSITIVITY_DIR] [OUT_DIR]
#     CHRONOGRAM       default trees/chronogram_treePL.nwk
#     SENSITIVITY_DIR  default trees/chronogram_sensitivity (one tree per file)
#     OUT_DIR          default results/CAFE
#
# Outputs (OUT_DIR):
#   cafe_input.tsv, cafe_input_pectin.tsv           CAFE inputs
#   <chronogram>/family/, <chronogram>/pectin/      CAFE output (final = the final chronogram)
#   <chronogram>/pectin_branches.tsv                aggregated pectin: count, change and p on every branch
#   significant_families.tsv                        lambda and families with p < 0.05 on every chronogram
#   key_branches.txt                                aggregated pectin on the branches discussed in the article
#   substrates/                                     LIG and fixed-lambda runs, per-family run (-P 1.0),
#                                                   key_nodes_substrates.txt, key_branches_per_family.txt
#
# Software: CAFE 5.1.0 (environment: envs/paper_cafe.yml); Python 3 (standard library).
# Note: CAFE 5 exits with status 0 even after errors such as "Invalid branch
# length", so its logs are checked below.
# =============================================================================
set -euo pipefail

CHRONOGRAM=${1:-trees/chronogram_treePL.nwk}
SENSITIVITY_DIR=${2:-trees/chronogram_sensitivity}
OUT_DIR=${3:-results/CAFE}
COUNTS=data/CAZyme_family_counts.csv
GROUPS=data/CAZyme_family_groups.tsv
GENOMES=data/genomes.tsv
CAFE=${CAFE5:-cafe5}
command -v "$CAFE" >/dev/null || { echo "ERROR: cafe5 not found (conda env create -f envs/paper_cafe.yml)"; exit 1; }
mkdir -p "$OUT_DIR/substrates"

check_log () { if grep -qiE "invalid|error|failed" "$1"; then echo "ERROR in $1"; exit 1; fi; }

# -----------------------------------------------------------------------------
# 1. CAFE inputs
# -----------------------------------------------------------------------------
python3 - "$COUNTS" "$GROUPS" "$OUT_DIR" <<'PY'
import csv, sys
counts_csv, groups_tsv, out = sys.argv[1:4]
MAX_SIZE, MIN_GENOMES = 100, 4
rows = list(csv.reader(open(counts_csv, newline="")))
names = [n.replace("-", "_") for n in (r[0] for r in rows[1:])]
groups = {}
for r in csv.DictReader(open(groups_tsv), delimiter="\t"):
    groups.setdefault(r["substrate"], []).append(r["family"])
kept = []
for j, fam in enumerate(rows[0][1:], start=1):
    v = [int(r[j]) for r in rows[1:]]
    if max(v) < MAX_SIZE and sum(x > 0 for x in v) >= MIN_GENOMES:
        kept.append([fam, fam] + v)
def write(path, table):
    with open(path, "w") as fh:
        fh.write("Desc\tFamily ID\t" + "\t".join(names) + "\n")
        for r in table:
            fh.write("\t".join(map(str, r)) + "\n")
def aggregate(table, families, label):      # sum the listed families into one row, appended at the end
    rest = [r for r in table if r[1] not in families]
    summed = [sum(col) for col in zip(*[r[2:] for r in table if r[1] in families])]
    return rest + [[label, label] + summed]
write(f"{out}/cafe_input.tsv", kept)
write(f"{out}/cafe_input_pectin.tsv", aggregate(kept, groups["pectin"], "PECTIN"))
write(f"{out}/substrates/cafe_input_lignin.tsv", aggregate(kept, groups["lignin"], "LIG"))
write(f"{out}/substrates/cafe_input_substrate_families.tsv",
      [r for r in kept if r[1] in sum(groups.values(), [])])
print(f"CAFE input: {len(kept)} of {len(rows[0]) - 1} families, {len(names)} genomes")
PY

# -----------------------------------------------------------------------------
# 2. Per-family and aggregated-pectin runs on every chronogram
# -----------------------------------------------------------------------------
run_set () {   # NAME TREE
  local d="$OUT_DIR/$1"; rm -rf "$d"; mkdir -p "$d"
  "$CAFE" -i "$OUT_DIR/cafe_input.tsv"        -t "$2" -o "$d/family" -k 1 -p -P 0.05 -c 4 > "$d/family.log" 2>&1
  "$CAFE" -i "$OUT_DIR/cafe_input_pectin.tsv" -t "$2" -o "$d/pectin" -k 1 -p -P 0.05 -c 4 > "$d/pectin.log" 2>&1
}
CHRONOGRAMS=(final)
run_set final "$CHRONOGRAM" &
for tree in "$SENSITIVITY_DIR"/*; do
  name=$(basename "${tree%.*}"); CHRONOGRAMS+=("$name")
  run_set "$name" "$tree" &
done
wait
for name in "${CHRONOGRAMS[@]}"; do check_log "$OUT_DIR/$name/family.log"; check_log "$OUT_DIR/$name/pectin.log"; done

# -----------------------------------------------------------------------------
# 3. Lignin aggregate, fixed-lambda runs and per-family branch tests (final chronogram)
# -----------------------------------------------------------------------------
S="$OUT_DIR/substrates"
LAMBDA=$(awk '/^Lambda:/ { print $2 }' "$OUT_DIR/final/family/Base_results.txt")
POISSON=$(sed -n 's/^Poisson lambda: \([0-9.e-]*\) .*/\1/p' "$OUT_DIR/final/family.log" | head -n 1)
echo "lambda $LAMBDA  root Poisson $POISSON" | tee "$S/params.txt"

"$CAFE" -i "$S/cafe_input_lignin.tsv" -t "$CHRONOGRAM" -o "$S/LIG" -k 1 -p -P 0.05 -c 4 > "$S/LIG.log" 2>&1
check_log "$S/LIG.log"

for agg in PECTIN:"$OUT_DIR/cafe_input_pectin.tsv" LIG:"$S/cafe_input_lignin.tsv"; do
  name=${agg%%:*}; input=${agg#*:}
  { head -n 1 "$input"; awk -F'\t' -v n="$name" '$2 == n' "$input"; } > "$S/cafe_input_${name}_only.tsv"
  "$CAFE" -i "$S/cafe_input_${name}_only.tsv" -t "$CHRONOGRAM" -o "$S/fixed_lambda/$name" \
          -l "$LAMBDA" -k 1 -p -P 0.05 -c 2 > "$S/fixed_lambda_$name.log" 2>&1
  "$CAFE" -i "$S/cafe_input_${name}_only.tsv" -t "$CHRONOGRAM" -o "$S/fixed_lambda_uniform_root/$name" \
          -l "$LAMBDA" -k 1 -P 0.05 -c 2 > "$S/fixed_lambda_uniform_root_$name.log" 2>&1
  check_log "$S/fixed_lambda_$name.log"; check_log "$S/fixed_lambda_uniform_root_$name.log"
done

"$CAFE" -i "$S/cafe_input_substrate_families.tsv" -t "$CHRONOGRAM" -o "$S/per_family" \
        -l "$LAMBDA" -p"$POISSON" -k 1 -P 1.0 -c 4 > "$S/per_family.log" 2>&1
check_log "$S/per_family.log"

# -----------------------------------------------------------------------------
# 4. Summary tables
# -----------------------------------------------------------------------------
python3 - "$OUT_DIR" "$GENOMES" "$GROUPS" "${CHRONOGRAMS[@]}" <<'PY'
import csv, os, re, sys

out, genomes_tsv, groups_tsv, chronograms = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
S = os.path.join(out, "substrates")
meta = {r["Name"].replace("-", "_"): r for r in csv.DictReader(open(genomes_tsv), delimiter="\t")}
groups = {}
for r in csv.DictReader(open(groups_tsv), delimiter="\t"):
    groups.setdefault(r["substrate"], []).append(r["family"])
LISTS = {"CELL": groups["(hemi-)cellulose"], "LIG": groups["lignin"], "PECTIN": groups["pectin"]}
NODE = re.compile(r"([^:,()]*)<(\d+)>(\*?)_(\d+)(?::([0-9.eE+-]+))?")

def load(run, family):
    """Reconstructed tree of one family from Base_asr.tre: node id -> count, tips, parent, branch p, significance"""
    line = [l for l in open(os.path.join(run, "Base_asr.tre")) if l.strip().startswith(f"TREE {family} ")][0]
    nwk = line.split("=", 1)[1].strip().rstrip(";")
    tab = [l.rstrip("\n").split("\t") for l in open(os.path.join(run, "Base_branch_probabilities.tab"))]
    row = [r for r in tab if r[0] == family]
    probs = {re.search(r"<(\d+)>", h).group(1): v for h, v in zip(tab[0][1:], row[0][1:])} if row else {}
    nodes, pos = {}, 0
    def parse():
        nonlocal pos
        kids = []
        if nwk[pos] == "(":
            pos += 1
            while True:
                kids.append(parse())
                if nwk[pos] == ",": pos += 1; continue
                if nwk[pos] == ")": pos += 1; break
        m = NODE.match(nwk, pos); pos = m.end(); nid = m.group(2)
        tips = [m.group(1)] if not kids else [t for k in kids for t in nodes[k]["tips"]]
        nodes[nid] = dict(id=nid, sig=bool(m.group(3)), count=int(m.group(4)), tips=tips, p=probs.get(nid, ""),
                          bl=float(m.group(5)) if m.group(5) else 0.0)
        for k in kids: nodes[k]["parent"] = nid
        return nid
    root = parse()
    return nodes, root

def mrca(nodes, sel):
    sel = set(sel)
    return min((n for n in nodes.values() if sel <= set(n["tips"])), key=lambda n: len(n["tips"]))

def parent_count(nodes, n):
    return nodes[n["parent"]]["count"] if n.get("parent") else None

def fmt(nodes, n):
    pc = parent_count(nodes, n); p = n["p"]
    ps = "%.3f" % float(p) if p not in ("", "N/A") else "NA"
    return "%s->%s p=%s%s" % ("-" if pc is None else pc, n["count"], ps, "*" if (n["sig"] and pc is not None) else "")

def branch_table(nodes, path):
    with open(path, "w") as f:
        f.write("node\tn_tips\ttips\tcount\tparent_count\tchange\tp\tsignificant\tbranch_length\n")
        for n in nodes.values():
            pc = parent_count(nodes, n)
            f.write(f"{n['id']}\t{len(n['tips'])}\t{';'.join(sorted(n['tips']))}\t{n['count']}\t"
                    f"{'' if pc is None else pc}\t{'' if pc is None else n['count'] - pc}\t{n['p']}\t"
                    f"{int(n['sig'] and pc is not None)}\t{n['bl']}\n")

def summary(nodes, root):
    sig = [n for n in nodes.values() if n["sig"] and n.get("parent")]
    gains = sum(n["count"] > parent_count(nodes, n) for n in sig); losses = sum(n["count"] < parent_count(nodes, n) for n in sig)
    return f"root {nodes[root]['count']}; significant branches {len(sig)} ({gains} gains, {losses} losses)"

# --- 4a. pectin branch tables and significant families on every chronogram ------
pectin = {}
with open(os.path.join(out, "significant_families.tsv"), "w") as f:
    f.write("chronogram\tlambda\tfamilies with p < 0.05\tpectin aggregate\n")
    for c in chronograms:
        nodes, root = load(os.path.join(out, c, "pectin"), "PECTIN")
        pectin[c] = nodes
        branch_table(nodes, os.path.join(out, c, "pectin_branches.tsv"))
        lam = [l.split()[1] for l in open(os.path.join(out, c, "family", "Base_results.txt")) if l.startswith("Lambda")][0]
        fams = [r[0] for r in (l.rstrip("\n").split("\t") for l in open(os.path.join(out, c, "family", "Base_family_results.txt")))
                if not r[0].startswith("#") and float(r[1]) < 0.05]
        f.write(f"{c}\t{lam}\t{' '.join(fams)}\t{summary(nodes, root)}\n")

# --- 4b. aggregated pectin on the branches discussed in the article --------------
final = pectin["final"]
alltips = sorted(max(final.values(), key=lambda n: len(n["tips"]))["tips"])
def tips(pred): return sorted(t for t in alltips if pred(t, meta.get(t, {})))
def genus(prefix): return lambda t, m: t.startswith(prefix)
def cls(*c): return lambda t, m: m.get("Class") in c
def order(*o): return lambda t, m: m.get("Order") in o
KEY = [
    ("Dothideomyceta crown (Dothideomycetes + Arthoniomycetes)", tips(cls("Dothideomycetes", "Arthoniomycetes"))),
    ("Dothideomycetes crown", tips(cls("Dothideomycetes"))),
    ("Trypetheliales crown", tips(order("Trypetheliales"))),
    ("Aureobasidium", tips(genus("Aureobasidium"))[:1]), ("Baudoinia", tips(genus("Baudoinia"))[:1]),
    ("Zasmidium", tips(genus("Zasmidium"))[:1]), ("Phyllosticta", tips(genus("Phyllosticta"))[:1]),
    ("Aeminium", tips(genus("Aeminium"))[:1]), ("Arthonia", tips(genus("Arthonia_"))[:1]), ("Alyxoria", tips(genus("Alyxoria"))[:1]),
    ("Arthoniomycetes crown", tips(cls("Arthoniomycetes"))),
    ("Pyrenula + Phaeomoniellales", tips(genus("Pyrenula")) + tips(genus("Phaeomoniella")) + tips(genus("Pseudophaeomoniella"))),
    ("Pyrenula crown", tips(genus("Pyrenula"))),
    ("Phaeomoniella + Pseudophaeomoniella", tips(genus("Phaeomoniella")) + tips(genus("Pseudophaeomoniella"))),
    ("Verrucariales", tips(order("Verrucariales"))),
    ("Chaetothyriomycetidae crown", tips(genus("Pyrenula")) + tips(order("Verrucariales", "Chaetothyriales"))),
    ("Eurotiomycetes crown", tips(cls("Eurotiomycetes"))),
    ("Lecanoromycetes crown", tips(cls("Lecanoromycetes"))),
    ("Umbilicariaceae", tips(genus("Lasallia")) + tips(genus("Umbilicaria"))),
    ("Trapelia", tips(genus("Trapelia"))[:1]),
    ("Acarosporales", tips(order("Acarosporales"))),
    ("OG clade (Ostropales, Graphidales, Gyalectales)", tips(order("Ostropales", "Graphidales", "Gyalectales"))),
    ("Sordariomyceta", tips(cls("Sordariomycetes", "Leotiomycetes"))),
    ("root", alltips),
]
with open(os.path.join(out, "key_branches.txt"), "w") as f:
    f.write("Aggregated pectin family: reconstructed count at the parent -> at the node, branch p (* p < 0.05)\n")
    f.write("%-58s %s\n" % ("branch leading to", " | ".join("%-18s" % c for c in chronograms)))
    for name, sel in KEY:
        f.write("%-58s %s\n" % (name, " | ".join("%-18s" % fmt(pectin[c], mrca(pectin[c], sel)) for c in chronograms)))

# --- 4c. lignin and pectin aggregates, root distributions, (hemi-)cellulose sums ---
AGG = [("LIG", os.path.join(S, "LIG")), ("PECTIN", os.path.join(out, "final", "pectin")),
       ("LIG", os.path.join(S, "fixed_lambda", "LIG")), ("PECTIN", os.path.join(S, "fixed_lambda", "PECTIN")),
       ("LIG", os.path.join(S, "fixed_lambda_uniform_root", "LIG")), ("PECTIN", os.path.join(S, "fixed_lambda_uniform_root", "PECTIN"))]
famrun = os.path.join(out, "final", "family")
famtrees = {}
for fam in sum(LISTS.values(), []):
    try: famtrees[fam] = load(famrun, fam)
    except IndexError: pass                # not in the CAFE input, or not present at the root
def fmt_sum(sub, sel):
    tot = par = 0; sigs = []
    for f in [f for f in LISTS[sub] if f in famtrees]:
        nodes, _ = famtrees[f]; n = mrca(nodes, sel); tot += n["count"]
        pc = parent_count(nodes, n)
        if pc is None: par = "-"; continue
        if par != "-": par += pc
        if n["sig"]: sigs.append("%s%s(%d->%d)" % (f, "+" if n["count"] > pc else "-", pc, n["count"]))
    return "%s->%s %s" % (par, tot, ("sig:" + ",".join(sigs)) if sigs else "")
D = {}
lines = [f"sum-{k}: {len([f for f in v if f in famtrees])} of {len(v)} listed families in the per-family run"
         for k, v in LISTS.items()]
for fam, run in AGG:
    nodes, root = load(run, fam); D[run] = (fam, nodes)
    lam = [x.strip() for x in open(os.path.join(run, "Base_results.txt")) if x.startswith("Lambda")][0]
    label = os.path.relpath(run, out)
    lines.append(f"{label} [{fam}]: {lam}; {summary(nodes, root)}")
    if not run.startswith(os.path.join(out, "final")):
        branch_table(nodes, os.path.join(run, f"{fam}_branches.tsv"))
KEYD = dict(KEY)
SUB_KEY = [(k, KEYD[k]) for k in (
    "root", "Dothideomyceta crown (Dothideomycetes + Arthoniomycetes)", "Dothideomycetes crown", "Trypetheliales crown",
    "Eurotiomycetes crown", "Chaetothyriomycetidae crown", "Pyrenula + Phaeomoniellales",
    "Phaeomoniella + Pseudophaeomoniella", "Pyrenula crown", "Verrucariales", "Lecanoromycetes crown",
    "OG clade (Ostropales, Graphidales, Gyalectales)")]
for name, sel in SUB_KEY:
    lines.append(f"\n## {name}")
    for run, (fam, nodes) in D.items():
        lines.append("   aggregate %-60s %s" % (f"{fam} {os.path.relpath(run, out)}", fmt(nodes, mrca(nodes, sel))))
    for sub in LISTS:
        lines.append("   sum       %-60s %s" % (f"{sub} (sum of the per-family reconstructions)", fmt_sum(sub, sel)))
open(os.path.join(S, "key_nodes_substrates.txt"), "w").write("\n".join(lines) + "\n")

# --- 4d. per-family branch p-values (-P 1.0) at the key branches -------------------
pf = os.path.join(S, "per_family")
new = {}
for l in open(os.path.join(pf, "Base_asr.tre")):
    if l.strip().startswith("TREE "): fam = l.split()[1]; new[fam] = load(pf, fam)[0]
tab = [l.rstrip("\n").split("\t") for l in open(os.path.join(pf, "Base_branch_probabilities.tab"))]
ids = [re.search(r"<(\d+)>", h).group(1) for h in tab[0][1:]]
P = {r[0]: dict(zip(ids, r[1:])) for r in tab[1:]}
diff = sum(1 for f in new if f in famtrees for n in new[f] if new[f][n]["count"] != famtrees[f][0][n]["count"])
lines = [f"families: {len(new)}; node counts differing from the per-family run on the final chronogram: {diff}"]
for name, sel in SUB_KEY[1:]:
    lines.append(f"\n## branch leading to {name}")
    for sub, fams in LISTS.items():
        fams = [f for f in fams if f in new]; s = []; par = tot = 0
        for f in fams:
            n = mrca(new[f], sel); pc = new[f][n["parent"]]["count"]; par += pc; tot += n["count"]
            p = P.get(f, {}).get(n["id"], "")
            if p not in ("", "N/A") and float(p) < 0.05: s.append("%s %d->%d p=%.3f" % (f, pc, n["count"], float(p)))
        lines.append(f"   {sub:7s} sum {par}->{tot}; families with branch p < 0.05: {'; '.join(s) if s else 'none'} (of {len(fams)})")
open(os.path.join(S, "key_branches_per_family.txt"), "w").write("\n".join(lines) + "\n")
print(f"summaries written to {out}")
PY
