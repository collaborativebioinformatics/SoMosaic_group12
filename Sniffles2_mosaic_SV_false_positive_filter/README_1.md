# Sniffles2 mosaic SV false-positive filter

A random forest that scores structural-variant calls from [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) and flags likely false positives, using only fields already present in the VCF — no realignment, no extra passes over the BAM.

It was trained on a HapMap6 ONT whole-genome callset labelled by [truvari](https://github.com/ACEnglish/truvari) against the MIMS truth set. On held-out chromosomes it keeps **99.2% of true positives while removing 32% of false positives**, cutting the callset's FP rate from 12.4% to 8.9%.

The filter is **annotative, not destructive**: every record gets an `INFO/RF_PROB` score and a `PASS` / `RF_FAIL` tag, and nothing is dropped unless you ask for it.


---

## Quick start

```bash
pip install -r requirements.txt

python sv_rf_filter.py \
    -m models/sv_rf_filter_final.joblib \
    -i sample.sniffles.vcf \
    -o sample.rf_filtered.vcf
```

```
sample.rf_filtered.vcf: 24184 records | PASS 21908 (90.6%) | RF_FAIL 2276 (9.4%) at RF_PROB >= 0.37
```

The output is the input VCF with two additions per record:

| field | meaning |
|---|---|
| `INFO/RF_PROB` | forest probability the call is a true positive, 0–1 |
| `FILTER` | `PASS` if `RF_PROB >= threshold`, else `RF_FAIL` |

Useful flags:

```bash
-t 0.50            # override the stored threshold (default 0.37)
--drop-fail        # actually remove failing records instead of tagging them
--neighbour-vcf X  # measure neighbourhood density against a different callset (see below)
```

Then, for a hard-filtered file:

```bash
bcftools view -f PASS sample.rf_filtered.vcf -Oz -o sample.pass.vcf.gz
```

### From Python

```python
import sv_rf_filter as rf

bundle = rf.load_model("models/sv_rf_filter_final.joblib")
calls  = rf.read_vcf("sample.sniffles.vcf")

# nb_ref must be the WHOLE callset -- see "Neighbourhood features" below
prob = rf.score_calls(calls, bundle, nb_ref=calls)

calls["RF_PROB"] = prob
calls["RF_PASS"] = prob >= bundle["thresh"]
print(calls.loc[~calls.RF_PASS, ["CHROM", "POS", "SVTYPE", "SVLEN", "RF_PROB"]].head())
```

### Input requirements

- **Single-sample Sniffles2 VCF** (v2.8.0 or compatible), uncompressed or bgzipped. Split multi-sample files first: `bcftools view -s SAMPLE`.
- **GRCh38 contig naming** (`chr1`…`chr22`, `chrX`). Calls on other contigs are still scored, with a warning.
- The INFO/FORMAT keys the features are built from: `SVTYPE`, `SVLEN`, `END`, `VAF`, `NM`, `STDEV_POS`, `STDEV_LEN`, `STRAND`, `COVERAGE`, `SUPPORT_SA`, `SUPPORT_LONG`, and `FORMAT` `DR`, `DV`, `GQ`. If a key is missing the tool warns which features went dead rather than failing silently.

---

## Choosing a threshold

`0.37` is the default because it was the highest cutoff still retaining 99% of true positives in cross-validation. It is one point on a curve, not a magic number — pick from the exchange rate you're willing to accept. Held-out chromosomes (chr12 + chr22, 1,219 TP / 173 FP):

| threshold | TP kept | TP lost | FP removed | resulting FP rate |
|---:|---:|---:|---:|---:|
| 0.25 | 1219 (100.0%) | 0 | 29 (16.8%) | 10.6% |
| **0.37** | **1209 (99.2%)** | **10** | **55 (31.8%)** | **8.9%** |
| 0.50 | 1186 (97.3%) | 33 | 86 (49.7%) | 6.8% |
| 0.65 | 1135 (93.1%) | 84 | 121 (69.9%) | 4.4% |
| 0.80 | 984 (80.7%) | 235 | 151 (87.3%) | 2.2% |

Baseline with no filter: 12.4% FP.

`RF_PROB` is a **ranking score, not a calibrated probability** — the forest was fit with `class_weight="balanced_subsample"`, so 0.37 does not mean "37% chance of being real". Compare scores against each other, not against an absolute prior.

---

## Features

Twenty features, all derived from the VCF itself. Absolute read-depth fields (`SUPPORT`, `DR`, `DV`, raw `COVERAGE`) were deliberately **excluded**: they are tied to one sample's sequencing depth, and a model that splits on `coverage <= 85` is meaningless on a callset sequenced four times deeper. Only ratios and shape statistics survive.

Importances are from the held-out fit; `corr(y)` is Spearman against the label on the full callset (positive = associated with true positives).

### Call quality and allele support

| feature | source | importance | corr(y) |
|---|---|---:|---:|
| `vaf` | `INFO/VAF` — variant allele fraction | 0.174 | +0.24 |
| `nm` | `INFO/NM` — mean edit distance of supporting reads | 0.078 | −0.01 |
| `depth` | `FORMAT/DR + DV` — total reads at the site | 0.068 | −0.08 |
| `cov_ratio` | `COVERAGE[center] / mean(COVERAGE[up], COVERAGE[down])` — local depth dip or spike | 0.063 | −0.08 |
| `qual` | `QUAL` | 0.025 | −0.00 |
| `gq` | `FORMAT/GQ` | 0.013 | −0.01 |
| `support_long` | `INFO/SUPPORT_LONG`, 0 when absent | 0.002 | −0.03 |
| `support_sa` | `INFO/SUPPORT_SA`, 0 when absent | 0.002 | +0.05 |

`vaf` is the single strongest feature, which is the expected story: a call supported by a small fraction of reads at a well-covered site is usually noise.

### Breakpoint precision

| feature | source | importance | corr(y) |
|---|---|---:|---:|
| `stdev_len_norm` | `STDEV_LEN / SVLEN` | 0.074 | −0.03 |
| `stdev_len` | `INFO/STDEV_LEN` — spread of the size estimate across reads | 0.057 | −0.07 |
| `stdev_pos` | `INFO/STDEV_POS` — spread of the breakpoint estimate | 0.041 | −0.08 |
| `stdev_pos_norm` | `STDEV_POS / SVLEN` | 0.039 | −0.06 |

Reads that disagree about where an event starts, or how long it is, are the signature of an alignment artifact. The `_norm` variants exist because 50 bp of disagreement means something different on a 200 bp insertion than on a 20 kb deletion.

### Size, type, orientation

| feature | source | importance | corr(y) |
|---|---|---:|---:|
| `svlen` | `INFO/SVLEN` — signed | 0.095 | −0.06 |
| `gc` | GC fraction of the inserted (ALT) or deleted (REF) allele; `NaN` for symbolic ALTs | 0.123 | +0.11 |
| `strand_both` | 1 if `INFO/STRAND == "+-"` — supported from both orientations | 0.008 | +0.10 |
| `is_DUP` / `is_DEL` / `is_INS` | one-hot `SVTYPE`. `INV` and `BND` fall through as all-zero | 0.007 / 0.005 / 0.005 | −0.19 / +0.07 / −0.02 |

### Neighbourhood features

| feature | meaning | importance | corr(y) |
|---|---|---:|---:|
| `nb_n_1kb` | other calls in the callset within ±1000 bp | 0.080 | **−0.30** |
| `nb_n_0kb` | other calls within ±500 bp (see naming note below) | 0.042 | −0.27 |

These two carry the strongest *negative* association with the label in the whole set: a call with neighbours crowded around it is much more likely to be a false positive. Artifacts cluster — tandem repeat arrays, segmental duplication edges, and any region where the aligner keeps producing inconsistent representations of one event.

**Both are measured against the entire callset**, TP and FP together. This matters. If density were computed inside `tp.vcf` alone, a true positive would look crowded only because the calls around it were also in `tp.vcf` — "number of neighbours" would become a partial restatement of the label. The union is what the caller actually emitted, which is also exactly what you have at inference time. A call is never its own neighbour (self-matches drop out by VCF ID).

At inference `sv_rf_filter.py` uses the input VCF as its own reference by default. **Point `--neighbour-vcf` at an equivalently filtered callset if your input is not comparable to the training reference** — the model learned densities over truvari's comparison set (19,503 records, PASS-only, inside the benchmark regions). A raw Sniffles VCF containing every non-PASS call and every BND will produce systematically higher counts than the model was trained on.

> **Naming note.** `nb_n_0kb` is the **500 bp** window. The training code generated column names with `f"nb_n_{w // 1000}kb"` over windows `(100, 500, 1000)`, so the 100 bp and 500 bp windows both resolved to `nb_n_0kb` and the 500 bp value overwrote the 100 bp one. The name is preserved here because the saved model's feature list uses it. Rename both sides together when retraining — and reclaim the 100 bp window, which is currently computed and thrown away.

---

## Feature correlation

<img width="1213" height="977" alt="feature_correlation" src="https://github.com/user-attachments/assets/83150a61-f8f6-49a7-81d3-41a7a8d7d2fb" />


Rows and columns are ordered by hierarchical clustering, so redundancy shows up as blocks on the diagonal. The right-hand panel is each feature's Spearman correlation with the label, in the same row order. Generated by the last cell of the notebook.

What to read out of it:

- **Two near-duplicate pairs.** `is_DEL` / `is_INS` at r = −0.98 (near-complementary once DUP is rare) and `nb_n_0kb` / `nb_n_1kb` at r = 0.95. Each pair carries roughly one feature's worth of information, and the forest splits importance between the members — which is why `nb_n_1kb` at 0.080 understates how much the neighbourhood block actually contributes.
- **The size/precision block.** `svlen`, `stdev_len`, `stdev_len_norm`, `stdev_pos_norm` and `is_INS` cluster together: insertions are shorter and their length estimates vary more.
- **Weak label correlation ≠ useless feature.** `gc` has a Spearman of only +0.11 with the label but the second-highest importance, because the forest uses it non-monotonically and in interaction with SV type. The notebook's ranking table reports per-feature AUC and mutual information alongside correlation for exactly this reason.

Because the plot covers all 22 chromosomes including the held-out ones, it is **descriptive only** — dropping a feature on the strength of it would leak the test set. Rerun on `X.iloc[tr_idx]` / `y_tr` for feature selection.

---

## How the model was built

**Data.** Sniffles2 v2.8.0 on `HapMap6_ONT_uwsc_SMAFIUU2RO1T_WGS`, GRCh38, median coverage ~65×. Truvari labels: 17,174 TP + 2,329 FP = 19,503 calls, an 11.9% FP rate. Truvari's own scoring keys (`TruScore`, `PctSeqSimilarity`, `MatchId`, …) are stripped at parse time so they cannot leak into the features.

**Split.** chr12 and chr22 held out entirely — 1,392 calls (1,219 TP / 173 FP), never touched during training, feature work, or threshold selection. Training used the remaining 20 chromosomes (18,111 calls).

**Validation.** Leave-one-chromosome-out cross-validation on the training chromosomes, so every call is scored by a forest that never saw its chromosome. This mirrors what the final held-out test does, which is what makes the OOF-selected threshold transfer.

**Model.** `RandomForestClassifier(n_estimators=500, min_samples_leaf=2, class_weight="balanced_subsample", random_state=42)`.

| | AUC | AP |
|---|---:|---:|
| Leave-one-chromosome-out (20 training chromosomes) | 0.892 | 0.982 |
| Held-out chr12 + chr22 | 0.913 | 0.984 |

Per-chromosome AUC across the CV folds ranged **0.851 – 0.921**. Any single-chromosome result inside that band is indistinguishable from chromosome-to-chromosome noise — worth remembering before reading much into one test chromosome.

**Threshold.** 0.37, the highest cutoff retaining 99% of TPs out-of-fold (99.04% TP kept, 27.3% FP removed). Stability checks: dropping any one chromosome moves it to 0.370–0.375; bootstrap 5–95% interval 0.355–0.385.

**Production model.** After the held-out test was read once and the threshold locked, the shipped model was refit on all 22 chromosomes. The bundle records its provenance:

```python
{"model": ..., "feats": [...], "thresh": 0.37,
 "trained_on": ["chr1", ..., "chr22"], "test_label": ["chr22", "chr12"],
 "oof_auc": 0.892, "heldout_auc": 0.913}
```

---

## Limitations

**One sample, one depth.** Everything here comes from a single ~65× ONT run. Chromosome holdout tests generalisation across genomic position; it says nothing about generalisation across samples, coverage, basecaller, or Sniffles version — which is what actually breaks in practice. `depth` is still an absolute read count and is the feature most likely to misbehave on a callset at a different depth. Before trusting this on a new sample, compare feature distributions against the training set; the notebook has a `shift_report` / `adversarial_auc` pair for exactly this. An adversarial AUC above ~0.65 means the transfer is not safe.

**The labels are truvari's, not ground truth.** "FP" means "did not match the truth set", which bundles together real artifacts, real variants the truth set is missing, and representation differences truvari could not reconcile. To the extent the third category is common, the filter is partly learning to reproduce truvari's bookkeeping rather than to detect wrong calls. Manual review of a subset is the only way to measure that gap.

**Sniffles2-specific.** The features read Sniffles INFO keys directly. Other callers will need remapping.

**INV and BND are effectively unmodelled.** They one-hot to all zeros and have no `SVLEN`, so they are scored but were never a target of feature design. Treat their scores with suspicion.

**Not calibrated.** See the threshold section.

---

## Retraining

Open `sv_fp_filter.ipynb`, point `TP_VCF_WGS` / `FP_VCF_WGS` at your truvari output, and run through. The notebook covers the split, leave-one-chromosome-out CV, threshold selection with stability checks, the held-out test, the correlation analysis, and a set of optional cells that compare model classes, test coverage-normalised features, and benchmark the forest against simple hand rules.

Two things to preserve if you change the feature set:

1. `sv_rf_filter.py` must emit exactly the columns in `bundle["feats"]`, by the same names. It raises rather than silently scoring on `-1` placeholders — pass `--no-strict` to override, but read the warning first.
2. Neighbourhood features must be computed against the full callset in a single pass, at training *and* inference. Scoring `tp.vcf` and `fp.vcf` separately halves every count.

> **Known issue in the notebook.** The in-notebook `apply_filter()` calls `featurize(df)` without `nb_ref=`, so `nb_n_0kb` and `nb_n_1kb` arrive as the `-1` sentinel — roughly 12% of the model's importance is inert in any VCF written by that function. `sv_rf_filter.py` does not have this bug. Fix the notebook by threading `nb_ref` through:
>
> ```python
> p_model = model.predict_proba(
>     featurize(df, nb_ref=nb_ref).reindex(columns=feats).fillna(-1))[:, 1]
> ```

**Version pinning matters.** The `.joblib` is a pickle of a scikit-learn estimator. Loading it under a different scikit-learn version may warn, misbehave, or fail — keep the pin in `requirements.txt` in step with whatever wrote the model.
