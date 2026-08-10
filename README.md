# TARGET-WT: which layers carry the variation in Wilms tumour?

Multi-omics factor analysis of the TARGET Wilms tumour cohort (GDC open access), asking a
data-appropriateness question before a biological one: **which axes of variation are shared
between the transcriptome and the methylome, and which are private to one layer?**

The project began as a somatic mutation analysis. The mutation data is reported here, but it
is presented as evidence for why the analysis moved to expression and methylation — not as a
result in its own right.

---

## 1. The coding view is the wrong instrument for this disease

Open-access masked somatic mutation calls were available for **38 cases**, carrying **377
variant rows** in total.

| Statistic | Value |
|---|---|
| Cases | 38 |
| Median coding SNVs per case | 6.5 (IQR 5–9, range 1–21) |
| Median TMB | 0.17 mutations/Mb |
| Genes mutated in ≥1 case | 273 |
| **Genes mutated in ≥10% of cases** | **1** |

At n = 38, a 10% prevalence filter requires four cases. Exactly one gene clears it: **TP53
(13/38, 34%)**. The next most frequent genes are BCOR and TTN at 3 cases each.

The same counting applied across four TARGET cohorts, processed through the
identical GDC pipeline:

| Project | Class | n | Median SNVs | % over 50 |
|---|---|---|---|---|
| TARGET-ALL-P3 | recurrent | 3 | 269 | 100% |
| TARGET-ALL-P2 | recurrent | 21 | 35 | 29% |
| TARGET-ALL-P2 | primary | 717 | 16 | 4% |
| **TARGET-WT** | **primary** | **38** | **7.5** | **0%** |

No Wilms tumour sample reaches 50 mutations — the rough floor below which
fitting to COSMIC reference signatures becomes unstable. (Counts here include
all variant classes; the coding-only median in the table above is 6.5.)

Two observations from the frequency table matter more than the burden statistics:

- **TP53's prominence is a cohort artifact.** TARGET-WT selected on anaplastic histology *or*
  relapsed favourable histology. Anaplastic Wilms tumour is defined partly by TP53 pathway
  disruption, so this frequency reflects the inclusion criteria, not Wilms tumour generally.
- **WT1 appears in 2/38 cases and does not survive the filter, while TTN — a canonical FLAGS
  passenger gene — appears in 3.** A coding-only analysis of this cohort ranks a sequencing
  artifact above the disease's namesake gene.

Mutational signature analysis was considered and rejected on the same grounds: de novo NMF
extraction needs mutation counts in the hundreds per sample, and fitting to COSMIC reference
signatures is unstable below roughly 50 mutations per sample.

A binary gene-level mutation matrix was the original plan for a third integration view. It was
dropped because one qualifying gene cannot constitute a view.

The defining lesion in Wilms tumour — loss of imprinting at 11p15 (IGF2/H19) — is epigenetic
and never appears in a MAF at any sequencing depth. The instrument had to change.

---

## 2. Cohort construction

Layer availability, TARGET-WT, GDC open access:

| Layer combination | Cases |
|---|---|
| expression + methylation | 123 |
| all three layers | 36 |
| expression only | 2 |
| methylation only | 1 |

The 38 open mutation cases are almost entirely a subset of the other two layers — requiring all
three views would have cost 87 cases to gain the single-gene matrix described above. The
two-layer design is therefore a data-availability finding, not a compromise.

Restricting to primary tumour aliquots (sample type `01`), one aliquot per case, and
intersecting the two layers gives the analysis cohort:

**n = 122 cases, expression + methylation.**

Normal-tissue aliquots (6, expression only) were dropped rather than modelled — six pairs
cannot support a tumour/normal analysis but are enough to dominate the leading factor.
Recurrent (`02`) and metastatic (`06`) aliquots were excluded for tissue-context consistency.

Feature matrices:

| View | Features after filtering | Features used |
|---|---|---|
| Expression (STAR counts) | 14,573 protein-coding | top 5,000 by variance |
| Methylation (HM450 beta) | 384,431 probes | top 5,000 by variance |

Expression: TMM normalisation, log-CPM. Methylation: complete-case probes only, autosomes only,
beta clipped to [0.001, 0.999] and converted to M-values.

---

## 3. Factor structure

MOFA2, 2 views, 12 factors requested, 10 retained. Total variance explained: **42.3%
(expression), 46.3% (methylation)**.

Every factor was refit across 5 seeds and matched by absolute correlation. **All 10 factors
reproduce with min |r| ≥ 0.95**, so none were discarded as fitting noise.

| Factor | Expression % | Methylation % | Character |
|---|---|---|---|
| 1 | 1.9 | 22.6 | methylation-dominant |
| 2 | 6.6 | 8.8 | **shared** |
| 3 | 12.1 | 2.0 | expression-dominant |
| 4 | 4.6 | 3.7 | shared |
| 5 | 5.4 | 1.1 | expression-leaning |
| 6 | 3.8 | 1.0 | expression-leaning |
| 7 | 4.6 | 0.02 | expression-private |
| 8 | 0.3 | 3.9 | **methylation-private** |
| 9 | 1.9 | 2.2 | shared |
| 10 | 2.4 | 1.5 | shared |

Factors were selected by property rather than by index, since MOFA factor numbering and sign
are not stable across refits: the *shared* factor is the one maximising the minimum variance
explained across both views (Factor 2), and the *methylation-private* factor is the one with
the highest methylation variance and <1% expression variance (Factor 8).

### A note on the MOFA quality-control warning

MOFA's built-in QC flagged Factors 1 and 3 as correlated with the number of expressed features.
This was checked directly against detection rates computed from the **raw** data (non-zero genes
per sample; non-NA probes per sample). Maximum absolute correlation across all ten factors was
**0.33**; Factors 1 and 3 were at 0.08. The warning does not reproduce against the raw metric —
MOFA's check counts exact zeros, which are essentially absent from log-CPM and M-value matrices.
No factor was excluded as technical. TMM normalisation was added upstream before this check was
run.

---

## 4. Are the shared factors coupled through promoter methylation?

For each factor, every promoter-associated CpG-island probe (TSS200, TSS1500, 5'UTR, 1stExon)
was matched to its annotated gene in the expression view, and the two per-factor weights were
correlated across all matched pairs (n = 580). A negative correlation is the expected signature
of promoter hypermethylation tracking reduced expression.

| Factor | Expr % | Meth % | r | p | BH-adjusted p |
|---|---|---|---|---|---|
| 10 | 2.4 | 1.5 | −0.357 | 7.0e−19 | 5.6e−18 |
| 9 | 1.9 | 2.2 | −0.269 | 4.2e−11 | 1.7e−10 |
| **2** | **6.6** | **8.8** | **−0.222** | 6.7e−08 | 1.8e−07 |
| 1 | 1.9 | 22.6 | −0.194 | 2.7e−06 | 5.3e−06 |
| 5 | 5.4 | 1.1 | −0.153 | 2.2e−04 | 3.5e−04 |
| 4 | 4.6 | 3.7 | +0.147 | 3.9e−04 | 5.1e−04 |
| 6 | 3.8 | 1.0 | −0.124 | 2.7e−03 | 3.1e−03 |
| 3 | 12.1 | 2.0 | −0.003 | 0.94 | 0.94 |

**The main result: variance and coupling are decoupled.** Factor 3 explains the most expression
variance in the cohort (12.1%) and shows no promoter-methylation coupling whatsoever (r =
−0.003). The dominant transcriptional axis is driven by something other than promoter
methylation. Meanwhile Factor 2 — the only factor that is both substantially shared and
meaningfully coupled — is the axis where the two layers co-vary *and* the mechanism is visible.

Factors 9 and 10 show stronger correlations but explain little variance in either view; a strong
correlation on a small factor is easy to over-read.

Factor 4's positive coefficient runs opposite to the repressive direction. With eight tests and
non-independent pairs, this is more plausibly noise than biology and is reported without
interpretation.

Genomic context of the top 100 probes:

| | Island | N_Shore | S_Shore | N_Shelf | S_Shelf | OpenSea |
|---|---|---|---|---|---|---|
| Factor 2 (shared) | 65 | 11 | 3 | 3 | 1 | 17 |
| Factor 8 (meth-private) | 55 | 8 | 6 | 6 | 4 | 21 |

The island enrichment for the shared factor is modest. The coupling result above is the stronger
evidence, since it uses complete weight vectors across 580 matched pairs rather than a truncated
top-100 list.

---

## 5. Limitations

- **The cohort is selected, not representative.** Inclusion required anaplastic histology or
  relapsed favourable histology, plus ≥80% tumour nuclei. Findings apply to high-risk Wilms
  tumour only. Relapse is partly an inclusion criterion, so it cannot be modelled as an outcome.
- **Effect sizes are modest.** r = −0.222 corresponds to roughly 5% of the variance in
  expression weights. This is coupling that exists, not coupling that dominates.
- **The 580 probe–gene pairs are not independent.** Probes map to multiple genes, genes carry
  multiple probes, and neighbouring CpGs co-vary. The p-values are anti-conservative; the effect
  sizes carry the claims.
- **Variance-based feature selection is not neutral.** It favours bimodally distributed features,
  which for methylation means imprinted and polymorphic loci are preferentially retained.
- **Gene symbol matching via `UCSC_RefGene_Name` is approximate.** The HM450 annotation predates
  current gene models and includes multi-gene probes.
- **No survival analysis yet.** Clinical data is aligned for all 122 cases with 50 deaths; the
  factor–histology and factor–survival tests are not part of this release.
- **Nothing here is a discovery.** Wilms tumour subtypes and the 11p15 imprinting mechanism are
  established. This repository is a study of layer concordance and modality choice, not of novel
  biology.

---

## 6. Repository

```
R/barcodes.R              TARGET barcode field parsing (shared)
00_data_inventory.R       layer availability, intersections, sample types
02_download_layers.R      GDC download, expression + methylation
03_build_matrices.R       aliquot selection, normalisation, feature selection
04_mutation_layer.R       burden, gene frequency, oncoplot (descriptive)
05_mofa_fit.R             MOFA2 fit, 5-seed stability
06_clinical_align.R       clinical join, endpoint construction
07_factor_features.R      detection QC, factor selection, coupling test
results/                  all tables and figures
```

Run in numeric order from the repository root. `data/`, `models/` and `GDCdata/` are
gitignored; rerunning `02` regenerates them.

TARGET barcodes place the sample type in the fourth hyphen-delimited field, and the patient code
is six characters rather than TCGA's four — fixed character offsets do not transfer. All barcode
handling goes through `R/barcodes.R` for this reason.

Package versions are recorded in `sessionInfo.txt`. The GDC data release used for download is
recorded in `data/download_provenance.txt` (not committed); the GDC portal displayed a
repository-under-review notice at the time of download.

**Data:** NCI TARGET Wilms Tumour (phs000471), accessed via the GDC Data Portal. Controlled-access
mutation calls (~593 additional cases) were not used; this analysis is reproducible from open
data alone.
