# Primer BLAST Filter & Amplicon Visualizer

This tool parses BLAST results of primer pairs against genome assemblies, identifies viable amplicons, and generates tabular and graphical outputs. It is designed for researchers validating primer specificity and distribution across genomes.

---

## 🔧 Installation (Conda Environment)

Clone the repo:

```bash
git clone https://github.com/YourUser/primer-blast-filter.git
cd primer-blast-filter
```

Create a Conda environment:

```bash
conda create -n primerfilter python=3.10 biopython matplotlib
conda activate primerfilter
```

Optionally install BLAST if you want to run it from the same environment:

```bash
conda install -c bioconda blast
```

---

## 📑 Input Files

### 1. Primer FASTA (`primer.fasta`)
Contains primer sequences in FASTA format. Primer IDs must end in `F` or `R` to mark forward/reverse orientation.

Example:
```
>MyGeneF
ATGCGTACGTTAGC
>MyGeneR
CGTACGACTTACGA
```

### 2. Genome FASTA (`genome.fasta`)
Reference genome(s) to test primers against.

### 3. BLAST Results (`blast.tsv`)
Tabular BLAST output (format 6 + sstrand + sseq).  
Run like this:

```bash
makeblastdb -in genome.fasta -dbtype nucl
blastn -query primer.fasta -db genome.fasta \
       -outfmt "6 qseqid sseqid sstart send sstrand pident length mismatch gapopen evalue bitscore sseq" \
       -out blast.tsv
```

---

## ▶️ Usage

```bash
python primer_filter.py \
  --primers primer.fasta \
  --blast blast.tsv \
  --genome genome.fasta \
  --out_prefix results \
  --min_len 80 \
  --max_len 3000 \
  --tick_units Mb \
  --tick_step 1
```

---

## 📂 Outputs

- `results_amplicons.tsv` → Viable amplicons with sequences, mismatches, % identity, orientation.
- `results_failed.tsv` → Primer hits that failed filters (with reasons).
- `results_plot.pdf` → Visual representation:
  - Black bar = chromosome
  - Green bar = viable amplicon (dark = normal, light = inverted)
  - Red tick = failed primer hit
  - Numbered labels tied to primer names in legend
 
![Example Output](PrimerBlast_Fusarium_odoratissimum_GCA_000260195.2_FO_II5_V1_plot.png)

---

## ⚙️ Parameters

### Required
- `--primers` : Primer FASTA file.
- `--blast` : BLAST results (format 6 + sstrand + sseq).
- `--genome` : Genome FASTA file.

### General
- `--out_prefix` : Prefix for output files (default `results`).

### Filtering
- `--min_len` : Minimum allowed amplicon length (default: 80 bp).
- `--max_len` : Maximum allowed amplicon length (default: 3000 bp).
- `--require_3p` : Number of 3′ bp that must match exactly (default: 3).
- `--max_mismatches` : Max mismatches allowed in primer hits (default: 10).
- `--min_pident` : Minimum % identity required (default: 0.0).
- `--len_tolerance` : Allowable difference (± bp) between primer length and BLAST hit length (default: 5).
- `--min_fail_len_frac` : Fraction of primer length required for a failed hit to be reported (default: 0.8).
- `--min_fail_pident` : Minimum % identity required for failed hits to be reported (default: 70.0).

### Plotting
- `--tick_units` : Units for x-axis ticks (`bp`, `kb`, `Mb`, `Gb`) (default: Mb).
- `--tick_step` : Spacing between ticks in chosen units (default: 1.0).

---

## 📊 Example Run

```bash
python primer_filter.py \
  --primers primer.fasta \
  --blast blast.tsv \
  --genome genome.fasta \
  --min_len 100 \
  --max_len 2000 \
  --tick_units kb \
  --tick_step 100
```

---

## 🔬 Features
- Handles **normal** and **inverted primer orientation**.
- Filters **partial alignments** (keeps only hits ~full primer length).
- Exports **amplicons** and **failed hits** separately.
- Publication-ready **PDF plots** with scaled chromosome axes.
- Configurable **axis units** (bp / kb / Mb / Gb) and tick spacing.
- Annotated legends mapping primers → numbers.

---

## 📖 Citation
Conrad R. (2025). BLAST_primer_filter. GitHub repository. Available at: https://github.com/rotheconrad/BLAST_primer_filter

---

# Overview: the filtering parameters and what they control

### A) Filters that affect *which BLAST hits are considered at all*

* `--len_tolerance`
* `--max_mismatches`
* `--min_pident`

These happen during BLAST TSV parsing. 

### B) Filters that affect *primer binding validity* (3′ end check)

* `--require_3p`

Applied after parsing, before pairing.

### C) Filters that affect *amplicon viability*

* `--min_len`
* `--max_len`

Applied after pairing F and R hits on same contig with valid orientations. 

### D) Filters that affect *what “failed hits” get reported*

* `--min_fail_len_frac`
* `--min_fail_pident`

These do **not** control viability. They only control whether a failed hit is written to `*_failed.tsv`.

---

# Walkthrough: what each filter does, in the order they’re applied

## 1) Primer ingest and pairing expectations (`load_primers`)

* Reads FASTA, requires IDs ending in `F` or `R`
* Uses the prefix as the “pair name” (e.g., `SIX1F` + `SIX1R` belong together)

If an ID doesn’t end in F/R, it hard-fails. 

**Why it matters:** every downstream filter assumes a primer pair is named `pairname + "F"` and `pairname + "R"`.

---

## 2) BLAST hit prefiltering (`parse_blast`)

This stage reads your BLAST outfmt line-by-line, creates hit records, and assigns a `fail_reason` for some hits.

### 2.1 `--len_tolerance` (silent discard of partial alignments)

```python
if abs(length - primer_len) > len_tolerance:
    continue
```

**What it means:** BLAST alignments that are too short/long compared to the primer length are **thrown away completely** and never appear in *either* amplicons or failed reports.

* This is the biggest “why didn’t I see it?” filter.
* It’s meant to exclude partial primer matches (common in noisy BLAST output).

**Example:** primer length 20, `len_tolerance=5`
Accepted alignment lengths: 15–25
Anything outside that is discarded.

---

### 2.2 `--max_mismatches` and `--min_pident` (mark as failed, but keep the hit)

```python
if mism > max_mismatches:
    reason = f"mismatches>{max_mismatches}"
elif pident < min_pident:
    reason = f"pident<{min_pident}"
```

**Key detail:** these do **not** discard the hit. They store it with `fail_reason`, and later `find_amplicons()` decides whether to record it in `*_failed.tsv` (using the “failed reporting” thresholds below).

So hits fall into two buckets coming out of `parse_blast()`:

* **Passing hits:** `fail_reason = None`
* **Failing hits:** `fail_reason` contains mismatch/pident reason

---

## 3) 3′-end exact matching (`three_prime_match` + `--require_3p`)

This is applied in `find_amplicons()` *only* to hits that didn’t already fail mismatch/pident.

```python
if strand == "plus":
    return primer_seq[-n:] == hit_seq[-n:]
else:
    return primer_seq[-n:] == hit_seq[:n]
```

**Interpretation:**

* For a **plus-strand** hit, it checks the last `n` bases of the primer against the **end** of the aligned subject string.
* For a **minus-strand** hit, it checks the last `n` bases of the primer against the **start** of the aligned subject string.

If it fails, the script converts it into a failed hit with:
`fail_reason = "3prime_mismatch (require n)"`

**Practical note:** this is a *hard exact-match* requirement for those last bases. If you’re looking for very distant primer-like matches, `require_3p=3` may be strict.

---

## 4) “Failed hit reporting” filters (`--min_fail_len_frac`, `--min_fail_pident`)

These apply only when a hit is already marked failed (either mismatch/pident from parsing, or 3′ mismatch later).

For a failed hit to be written to `*_failed.tsv`, it must pass:

### 4.1 `--min_fail_len_frac`

```python
if F["length"] < len(fwd_seq) * min_fail_len_frac:
    continue
```

This prevents tiny junk alignments from cluttering the failed report.

### 4.2 `--min_fail_pident`

```python
if F["pident"] < min_fail_pident:
    continue
```

This prevents super-low identity hits from cluttering the failed report.

**Again:** these do **not** decide if an amplicon is viable. They only decide if you *see* the failure in the output.

---

## 5) Pairing logic + orientation filter

After F and R hits individually pass the non-failure checks (and 3′ check), the script tries all F×R combinations and enforces:

### 5.1 Same contig

```python
if F["contig"] != R["contig"]:
    continue
```

### 5.2 Orientation must be “normal” or “inverted”

```python
valid_normal   = (F["strand"] == "plus"  and R["strand"] == "minus")
valid_inverted = (F["strand"] == "minus" and R["strand"] == "plus")
```

So it allows:

* **Normal PCR layout:** F on +, R on –
* **Inverted layout:** F on –, R on + (this can happen depending on how primers/hits are represented)

If neither, the pair is rejected silently.

---

## 6) Amplicon length filter (`--min_len`, `--max_len`)

Amplicon span is computed from the outermost hit coordinates:

```python
start = min(F["start"], F["end"], R["start"], R["end"])
end   = max(F["start"], F["end"], R["start"], R["end"])
length = end - start + 1
```

Then:

```python
if min_len <= length <= max_len:
    ... record amplicon ...
else:
    ... record failure amplicon_len=... out_of_bounds ...
```

If it passes, it extracts sequence from the genome FASTA and writes it to the amplicon table. 

---

# The “differences” between the filtering parameters in one sentence each

* **`len_tolerance`**: hard, silent discard of partial alignments (never shows up anywhere). 
* **`max_mismatches` / `min_pident`**: mark hits as failed (may still be reported and plotted as failed). 
* **`require_3p`**: hard exact-match at the primer’s 3′ end; failing hits become “3prime_mismatch”.
* **`min_len` / `max_len`**: filter on the *product size* after pairing; out-of-range becomes a pair-level failure. 
* **`min_fail_len_frac` / `min_fail_pident`**: only controls whether failed hits get *written* to the failed report (verbosity control). 

---

## Two common “gotchas” in this script (worth knowing)

1. If you’re hunting distant/partial primer matches, **`len_tolerance` can hide them completely** (because partials are silently dropped). 
2. You can have lots of “failed hits” in reality, but see none because **`min_fail_pident` is 70% by default** and you’re looking for much lower identity signals. 

---

