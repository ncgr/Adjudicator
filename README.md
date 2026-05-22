# Adjudicator

> **Version:** 0.1.0 · **License:** Apache-2.0 · **Status:** Beta  
> **PyPI:** `pip install adjudicator` · **Source:** [ncgr/Adjudicator](https://github.com/ncgr/Adjudicator)  
> **Author:** Connor Cameron · ctc@ncgr.org · National Center for Genome Resources

---

## Overview

Adjudicator is a command-line tool for collapsing and filtering structural genome annotations across multiple sources. It uses best-hit HMM domain scores from gene family assignments (via the [Legume Information System](https://data.legumeinfo.org/LEGUMES/Fabaceae/genefamilies/legume.fam3.VLMQ/)) to compare overlapping gene models and select the best-supported annotation for all overlapping models.

**Two commands are provided:**

- **`collapse`** — Merge overlapping gene models from two or more annotators, selecting the best model per region.
- **`repeat-filter`** — Remove gene models that overlap known repeat or transposon regions beyond a configurable coverage threshold.

---

## Requirements

| Requirement | Version |
|---|---|
| Python | ≥ 3.10 |
| click | ≥ 8.1 |
| intervaltree | ≥ 3.2.1 |
| sortedcontainers | ≥ 2.4.0 |

---

## Installation

```bash
pip install adjudicator
```

```bash
adjudicator --version
# adjudicator, version 0.1.0
```

---

## Input File Formats

### TSV Sample Sheet (`--input-tsv`)

Tab-separated. Lines beginning with `#` and blank lines are skipped.

| Column | Type | Description |
|--------|------|-------------|
| 1 | string | Unique label for this evidence set. |
| 2 | path | Path to the `.gff3` structural annotation file. |
| 3 | path | Path to the `.gfa` LIS gene family assignment file. |

Row order determines precedence when gene models have equivalent scores.

```tsv
# label	gff3_path	gfa_path
maker	/data/ann/maker.gff3	/data/fam/maker.gfa
helixer	/data/ann/helixer.gff3	/data/fam/helixer.gfa
stringtie	/data/ann/stringtie.gff3	/data/fam/stringtie.gfa
```

### GFF3 (`.gff3`)

Standard GFF3 format with a three-level hierarchy: `gene` → `mRNA` → `exon`. See the [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

### GFA (`.gfa`)

LIS gene family assignment files produced by the [Legume Information System gene family pipeline](https://data.legumeinfo.org/LEGUMES/Fabaceae/genefamilies/legume.fam3.VLMQ/).

---

## Commands

### `collapse`

Collapses overlapping structural annotations across all entries in the TSV. Processing is hierarchical: Row 1 vs. Row 2 produces an intermediate result, which is then compared against Row 3, and so on.

#### Synopsis

```bash
adjudicator collapse --input-tsv <FILE> [OPTIONS]
```

#### Options

| Option | Short | Type | Default | Valid range | Description |
|--------|-------|------|---------|-------------|-------------|
| `--input-tsv` | `-i` | path | *(required)* | — | Tab-separated sample sheet. |
| `--min-overlap` | `-m` | float | `0.00001` | 0.0 – 1.0 | Minimum fractional overlap of feature A by feature B to consider them overlapping. |
| `--no-orphans` | `-n` | flag | `False` | — | Exclude genes with no gene family assignment from the output. |
| `--output-dir` | `-o` | path | `.` | — | Directory to write output files. Created if it does not exist. |
| `--strict` / `--no-strict` | | flag | `False` | — | Exit with error if any referenced input file does not exist on disk. |
| `--verbose` | `-v` | flag | `False` | — | Print per-sample file paths and processing steps to stdout. |

#### Output Files

```
<output-dir>/
├── A_B.wao.gff3                       # Overlap intersections
├── A_B.unique_b.gff3                  # Gene models unique to annotator B
├── A_B.final.gff3                     # Adjudicated gene IDs
├── A_B.gfa                            # Merged gene family assignments
└── A_B.final.wsubfeatures.gff3        # ✅ Primary output
```

#### Examples

```bash
adjudicator collapse \
    --input-tsv samples.tsv \
    --output-dir results/collapse/
```

```bash
adjudicator collapse \
    --input-tsv samples.tsv \
    --no-orphans \
    --min-overlap 0.4 \
    --output-dir results/collapse/ \
    --verbose
```

---

### `repeat-filter`

Filters gene models from each entry in the TSV against a reference repeat annotation. Gene models whose exons exceed `--max-coverage` overlap with a repeat region are removed.

#### Synopsis

```bash
adjudicator repeat-filter --input-tsv <FILE> --annotation <FILE> [OPTIONS]
```

#### Options

| Option | Short | Type | Default | Valid range | Description |
|--------|-------|------|---------|-------------|-------------|
| `--input-tsv` | `-i` | path | *(required)* | — | Tab-separated sample sheet. |
| `--annotation` | `-a` | path | *(required)* | — | GFF3 file of repeat regions to filter against. |
| `--max-coverage` | `-m` | float | `0.4` | 0.0 – 1.0 | Maximum fractional overlap between a gene's exons and a repeat region before the model is removed. |
| `--output-dir` | `-o` | path | `.` | — | Directory to write output files. Created if it does not exist. |
| `--strict` / `--no-strict` | | flag | `False` | — | Exit with error if any referenced input file does not exist on disk. |
| `--verbose` | `-v` | flag | `False` | — | Print per-sample file paths and processing steps to stdout. |

#### Output Files

```
<output-dir>/
├── <label>_repeat_filter.wao.gff3                  # Overlap intersections
└── <label>_repeat_filter.final.wsubfeatures.gff3   # ✅ Primary output
```

#### Examples

```bash
adjudicator repeat-filter \
    --input-tsv samples.tsv \
    --annotation repeats.gff3 \
    --output-dir results/filtered/
```

```bash
adjudicator repeat-filter \
    --input-tsv samples.tsv \
    --annotation transposons.gff3 \
    --max-coverage 0.3 \
    --output-dir results/filtered/ \
    --strict \
    --verbose
```

---

## Workflow

```bash
# Step 1: Filter repeat regions
adjudicator repeat-filter \
    --input-tsv raw_samples.tsv \
    --annotation repeats.gff3 \
    --output-dir step1_filtered/

# Step 2: Rewrite TSV to point to filtered outputs (GFA paths unchanged)
awk -F'\t' 'OFS="\t" { $2="step1_filtered/"$1"_repeat_filter.final.wsubfeatures.gff3"; print }' \
    raw_samples.tsv > filtered_samples.tsv

# Step 3: Collapse filtered annotations
adjudicator collapse \
    --input-tsv filtered_samples.tsv \
    --output-dir step2_collapsed/
```

---

## Error Reference

| Condition | `--strict` off | `--strict` on |
|-----------|---------------|---------------|
| Input file not found | Warning to stderr | `Error: The following files were not found: ...` |
| TSV row has wrong column count | `BadParameter: Line N: expected 3 columns, got N.` | Same |
| Label (column 1) is empty | `BadParameter: Line N: column 1 (label) must not be empty.` | Same |
| GFF3 path does not end in `.gff3` | `BadParameter: Line N: column 2 must end in '.gff3'` | Same |
| GFA path does not end in `.gfa` | `BadParameter: Line N: column 3 must end in '.gfa'` | Same |
| TSV contains no data rows | `Error: No data rows found in '<file>'.` | Same |

---

## Glossary

| Term | Definition |
|------|------------|
| **Gene model** | A predicted gene structure represented as a `gene` → `mRNA` → `exon` hierarchy in GFF3. |
| **GFF3** | Generic Feature Format version 3. Tab-delimited format for genomic features and their hierarchical relationships. |
| **GFA** | Gene Family Assignment file from the LIS pipeline, containing HMM domain scores used to rank competing gene models. |
| **Adjudication** | Selection of one gene model from a set of overlapping candidates based on HMM score evidence. |
| **Orphan gene** | A gene model with no gene family assignment in the GFA file. |
| **WAO intersection** | A bedtools-style "write all overlaps" operation reporting fractional overlap between features across two GFF3 files. |
