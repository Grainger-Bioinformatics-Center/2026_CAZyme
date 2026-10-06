# 2026_CAZyme

**Plant cell wall-degrading enzymes in lichen-forming fungi are a legacy of plant-associated ancestry**

Scripts, input data and phylogenetic trees for:

> Sun Y, Lumbsch HT, Sangvichien E, Nelsen MP, Ametrano CG, Grewe F. Plant cell wall-degrading enzymes in lichen-forming fungi are a legacy of plant-associated ancestry. *MycoKeys* (submitted).

The study compares the carbohydrate-active enzyme (CAZyme) repertoires of 94 Pezizomycotina genomes, 51 of them from lichen-forming fungi. All genomes were re-annotated with one pipeline. The plant cell wall-degrading enzyme (PCWDE) families acting on (hemi-)cellulose, pectin and lignin were counted in each genome. Their evolution was then analysed on a maximum-likelihood supermatrix phylogeny with CAFE 5 and phylogenetic principal component analysis (pPCA).

## Scripts

There is one script for each section of the Materials and Methods, numbered in the order in which they are run.

| Script | Materials and Methods | Description |
|---|---|---|
| `scripts/1_metagenome_decontamination.sh` | Metagenome decontamination | Classifies contigs with Kraken 2 against NCBI nt, then removes contigs assigned to Bacteria, Archaea, Viruses or Viridiplantae with SeqKit. |
| `scripts/2_gene_annotation_and_CAZyme_groups.sh` | Gene identification and annotation; CAZyme groups identification and selection | Runs funannotate gene prediction and functional annotation (dbCAN, Pfam, MEROPS, eggNOG-mapper) and BUSCO on the predicted proteomes. Counts CAZyme genes per family and the (hemi-)cellulose, pectin and lignin genes per genome. |
| `scripts/3_phylogenetic_reconstruction.sh` | Phylogenetic reconstruction | Builds the supermatrix tree: OrthoFinder single-copy orthogroups, MAFFT, trimAl, AMAS and RAxML-NG. Roots it on Sordariomyceta. Builds the corroborating tree from the 195 UnFATE markers (miniprot) and computes the Robinson–Foulds distance between the two trees. |
| `scripts/4_ultrametric_transformation.sh` | Ultrametric transformation of the phylogeny | Dates the tree in treePL with two calibrations fixed at their median ages and the smoothing value chosen by cross-validation. Also produces five sensitivity chronograms. |
| `scripts/5_gene_family_evolution_CAFE.sh` | Gene-family evolution (CAFE) | Runs CAFE 5 per family and on the aggregated pectin family on all six chronograms. Analyses lignin and (hemi-)cellulose and tests other root-size distributions. Writes summary tables. |
| `scripts/6_phylogenetic_PCA.sh` | Phylogenetic PCA and statistical analyses | Runs the pPCA of the PCWDE families (`phytools::phyl.pca`) and estimates the phylogenetic signal of the families and of the pPCA scores. Tests BUSCO completeness against gene counts with Spearman's rank correlation. |

All scripts are bash scripts; steps that need Python 3 or R are included in them. Each script begins with a header that lists its inputs, outputs, settings and software versions. Run the scripts from the repository root.

## Repository contents

| Path | Content |
|---|---|
| `data/genomes.tsv` | The 94 analysed genomes and the second strains of two species, which were used only in the supermatrix tree. Columns: name used in all files (`Name`), name used in the article (`Species`), class, order, lifestyle, photobiont, culture or metagenome source, GenBank accession, BioProject and assembly. |
| `data/CAZyme_family_counts.csv` | Genes per CAZyme family and genome (output of script 2; Supplementary Table S4). |
| `data/genome_statistics.csv` | Per-genome statistics (output of script 2; Supplementary Table S1): BUSCO completeness, assembly statistics, number of genes, CAZyme totals per class, (hemi-)cellulose (`cell`), pectin (`pec`) and lignin (`lign`) genes, and sugar transporters (Pfam PF00083). |
| `data/CAZyme_family_groups.tsv` | The 49 PCWDE families and the substrate assigned to each (Supplementary Table S3). |
| `data/outgroup_Sordariomyceta.txt` | The 17 Sordariomycetes and Leotiomycetes on which the tree is rooted. |
| `trees/supermatrix_ML_96taxa.raxml.support` | Maximum-likelihood supermatrix tree with bootstrap support (96 genomes). |
| `trees/supermatrix_ML_rooted_94taxa.nwk` | The same tree rooted on Sordariomyceta and pruned to the 94 analysed genomes; used for the pPCA and the dating. |
| `trees/UnFATE_195markers.raxml.support` | Tree from the 195 UnFATE markers, with bootstrap support. |
| `trees/chronogram_treePL.nwk` | Final chronogram (Figure 4), used for CAFE. |
| `trees/chronogram_sensitivity/` | The five sensitivity chronograms: smoothing 0.1 and 10; calibrations spanning their 95% HPDs (`HPD_bounds`, identical to `HPD_upper`); both nodes at the lower (`HPD_lower`) or upper (`HPD_upper`) ends of their HPDs. |
| `envs/` | Conda environments with the exact versions used: `kraken` (script 1), `funannotate` (script 2), `phy` (BUSCO and script 3), `paper_cafe` (script 5), `lichen_ppca` (R steps of scripts 3, 4 and 6). |

The analysis results are not included. Scripts 4–6 recreate them in `results/` from the files in `data/` and `trees/`; the values they should reproduce are listed under Expected results.

## Requirements

| Analysis | Software |
|---|---|
| Decontamination | Kraken 2 v2.1.3, SeqKit v2.13.0 |
| Annotation | funannotate v1.8.17 with AUGUSTUS v3.5.0, GeneMark-ES v4.71, SNAP 2013-11-29, GlimmerHMM v3.0.4, EVidenceModeler v1.1.1, tRNAscan-SE v2.0.12, tantan v49; dbCAN 12.0, Pfam 36.0, MEROPS 12.0, eggNOG-mapper v2.1.12 |
| Proteome completeness | BUSCO v6.0.0 in protein mode, ascomycota_odb10 |
| Phylogeny | OrthoFinder v3.1.2, MAFFT v7.526, trimAl v1.5, AMAS v1.0, RAxML-NG v2.0.1, miniprot v0.18 |
| Dating | treePL (https://github.com/blackrim/treePL, commit f41af04, compiled with the bundled ADOL-C 2.6.3 and NLopt 2.4.2) |
| Comparative analyses | CAFE v5.1.0; R v4.5.3 with ape 5.8-1, phytools 2.5-2 and phangorn 2.12.1; Python 3 (standard library) |

GeneMark-ES requires a licence key from its distributor. treePL is not available from conda and has to be compiled.

Scripts 1–3 also need the following external data:

- **Genome assemblies.** Their accessions are listed in `data/genomes.tsv` and in Supplementary Tables S1 and S2. Ninety-three assemblies are available from NCBI GenBank. The assembly of *Fitzroyomyces cyperi* CBS 143170 is available from the JGI MycoCosm portal.
- **Kraken 2 database.** Built from NCBI nt and the NCBI taxonomy (January 2026).
- **funannotate databases.** dbCAN 12.0, Pfam 36.0, MEROPS 12.0, UniProt 2024_01, eggNOG 5.0.2 and BUSCO ascomycota_odb10.
- **UnFATE markers.** `UnFATE_markers_195.fas` from https://github.com/claudioametrano/UnFATE.

## Running the analyses

```bash
# 1  decontamination of the assemblies that need screening (Kraken 2 database in KRAKEN_DB)
KRAKEN_DB=/path/to/kraken2_db bash scripts/1_metagenome_decontamination.sh  metagenome_assemblies/

# 2  annotation of all 96 assemblies, one file per genome named <Name>.fasta as in data/genomes.tsv
#    (the filtered assemblies from step 1 and the other assemblies in one folder)
FUNANNOTATE_DB=/path/to/funannotate_db bash scripts/2_gene_annotation_and_CAZyme_groups.sh  assemblies/

# 3  supermatrix tree, UnFATE tree and their comparison
bash scripts/3_phylogenetic_reconstruction.sh  assemblies/  UnFATE_markers_195.fas

# 4-6 run from the files in data/ and trees/
bash scripts/4_ultrametric_transformation.sh
bash scripts/5_gene_family_evolution_CAFE.sh
bash scripts/6_phylogenetic_PCA.sh
```

Each script writes to its own folder under `results/`. Scripts 4–6 read the published trees and tables in `data/` and `trees/` by default. To run them on new outputs of scripts 2–4, pass the corresponding files as arguments or copy them to those folders; the header of each script gives the arguments.

## Expected results

| Analysis | Value |
|---|---|
| Supermatrix | 109 single-copy orthogroups; 96 taxa × 57,900 aligned amino-acid sites |
| Tree comparison | Robinson–Foulds distance 12 between the supermatrix and UnFATE trees (normalized 0.065; 93.5% of splits shared; 96 taxa) |
| Chronogram | root 458 Ma and Lecanoromycetes crown 248.6 Ma (fixed); smoothing 1 by cross-validation (χ² 71,272.5); identical with the seeds 20261002, 7 and 99 |
| CAFE, per family | 169 families present at the root; λ = 0.0026040; six families with p < 0.05: AA1, AA3, AA7, AA9, GH28, GH43 (the same on all six chronograms) |
| CAFE, aggregated pectin | root 9 genes; 49 branches with p < 0.05 (27 gains, 22 losses); 48–51 on the sensitivity chronograms |
| CAFE, other root distributions | pectin aggregate alone with λ fixed: root 15 genes with an estimated or a uniform root distribution; the gain on the branch leading to *Pyrenula* and the Phaeomoniellales becomes p = 0.086 and the loss leading to *Alyxoria varia* p = 0.038 |
| PCWDE pPCA | 94 genomes × 46 families; Pagel's λ = 0.895; PC1 35.6%, PC2 7.0%; 45 of 46 families with significant Pagel's λ |
| BUSCO vs gene counts (51 lichenized genomes) | pectin ρ = −0.076, p = 0.598; PCWDE ρ = −0.016, p = 0.912 |

## Notes

- **Names.** All files use the names under which the genomes were analysed. Four species are named differently in the article or at NCBI (column `Species` of `data/genomes.tsv` gives the name used in the article):
  - *Lasallia hispanica* is *Umbilicaria hispanica*.
  - *Phialophora attae* is *Cyphellophora attinorum*.
  - *Fitzroyomyces cyperi* CBS 143170 is listed by NCBI as *Ebollia carnea*.
  - *Talaromyces cellulolyticus* (GCA_000829775.1) is listed by NCBI as *Talaromyces pinophilus*.
- **Decontamination.** Script 1 applies the rule used for 54 assemblies. Two earlier annotation batches used other rules: 18 and 5 assemblies. Seventeen assemblies were not screened. The header of script 1 describes these batches. Re-applying the rule of script 1 to the classifications of all 77 screened assemblies flagged no further contigs carrying CAZyme genes.
- **Annotation batches.** The header of script 2 lists two small differences among annotation batches:
  - contig renaming before gene prediction;
  - the minimum contig length for *Epichloe typhina*.
- **Second strains.** The supermatrix tree was inferred from 96 proteomes. These include second strains of *Aureobasidium pullulans* (GCA_003336255.1) and *Trypethelium eluteriae* (GCA_051942245.1), which were removed for all comparative analyses.
- **Reproducibility of the trees.** To reproduce the published supermatrix and trees, script 3 excludes orthogroup OG0002496, which the original run did not read, and fixes the RAxML-NG seeds recorded in the logs of the published runs.
- **pPCA axes.** The signs of the pPCA axes are arbitrary. Figure 3 shows PC2 with the opposite sign.
- **CAFE logs.** CAFE 5 exits with status 0 even after errors such as "Invalid branch length". Script 5 therefore searches the CAFE logs for errors.

## Citation

If you use these scripts or data, please cite the article above. Volume, article number and DOI will be added on publication.

## License

[To be added by the authors.]
