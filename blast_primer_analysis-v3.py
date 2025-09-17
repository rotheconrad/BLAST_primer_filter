#!/usr/bin/env python3
"""
Primer BLAST Filter & Amplicon Visualizer
=========================================

This tool parses BLAST results of primer pairs against genome assemblies,
filters valid hits, identifies amplicons, and generates tables and plots.

Installation
------------
1. Clone the repo
   git clone https://github.com/YourUser/primer-blast-filter.git
   cd primer-blast-filter

2. Create Conda environment
   conda create -n primerfilter python=3.10 biopython matplotlib
   conda activate primerfilter

3. (Optional) Install BLAST
   conda install -c bioconda blast

Inputs
------
- primer.fasta : primers, IDs ending in F/R
- genome.fasta : reference genome(s)
- blast.tsv    : BLAST results (outfmt 6 + sstrand + sseq)

   Example BLAST:
   makeblastdb -in genome.fasta -dbtype nucl
   blastn -query primer.fasta -db genome.fasta \
          -outfmt "6 qseqid sseqid sstart send sstrand pident length mismatch gapopen evalue bitscore sseq" \
          -out blast.tsv

Outputs
-------
- *_amplicons.tsv : viable products with sequences
- *_failed.tsv    : failed primer hits (with reasons)
- *_plot.pdf      : graphical overview

Parameters
----------
--primers            Primer FASTA (required)
--blast              BLAST results in format 6 (required)
--genome             Genome FASTA (required)
--out_prefix         Output prefix [default: results]

Filtering:
--min_len            Minimum amplicon length [default: 80]
--max_len            Maximum amplicon length [default: 3000]
--require_3p         Bases at 3′ that must match [default: 3]
--max_mismatches     Max mismatches allowed [default: 10]
--min_pident         Min % identity [default: 0.0]
--len_tolerance      ±bp tolerance for primer vs hit length [default: 5]
--min_fail_len_frac  Fraction length for failed hits [default: 0.8]
--min_fail_pident    Min % identity for failed hits [default: 70.0]

Plotting:
--tick_units         Units for x-axis ticks: bp, kb, Mb, Gb [default: Mb]
--tick_step          Spacing between ticks in chosen units [default: 1.0]

Features
--------
- Normal + inverted orientations
- Filters partial hits
- Summary TSVs + PDF plots
- Axis ticks in bp/kb/Mb/Gb
- Shaded background for readability

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Sep 2025
License :: GNU GPLv3
Copyright 2025 Roth Conrad
All rights reserved
-------------------------------------------
"""

import argparse
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt


def load_primers(fasta_file):
    primers = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        if record.id.endswith(("F", "R")):
            pairname = record.id[:-1]
            orient = record.id[-1]
            primers.setdefault(pairname, {})[orient] = seq
        else:
            raise ValueError(f"Primer ID {record.id} must end with F or R")
    return primers


def three_prime_match(primer_seq, hit_seq, strand, require_n=3):
    if require_n == 0:
        return True
    primer_seq, hit_seq = primer_seq.upper(), hit_seq.upper()
    if strand == "plus":   # forward
        return primer_seq[-require_n:] == hit_seq[-require_n:]
    else:                  # reverse
        return primer_seq[-require_n:] == hit_seq[:require_n]


def parse_blast(blast_file, max_mismatches, min_pident,
                primer_lengths, len_tolerance=5):
    """
    Parse BLAST results, prefiltering:
      - discard partial matches (alignment length not within ±len_tolerance of primer length)
      - flag hits exceeding mismatches or falling below % identity
    """
    hits = defaultdict(list)
    with open(blast_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            qid, sid, sstart, send, strand, pident, length, mism, gaps, evalue, bits, sseq = line.strip().split("\t")
            sstart, send, length, mism = map(int, [sstart, send, length, mism])
            pident = float(pident)

            primer_len = primer_lengths.get(qid, None)
            if primer_len is None:
                continue

            # Prefilter: discard partial alignments silently
            if abs(length - primer_len) > len_tolerance:
                continue

            reason = None
            if mism > max_mismatches:
                reason = f"mismatches>{max_mismatches}"
            elif pident < min_pident:
                reason = f"pident<{min_pident}"

            hit = {
                "contig": sid, "start": sstart, "end": send, "strand": strand,
                "pident": pident, "length": length, "mismatches": mism,
                "sseq": sseq, "fail_reason": reason
            }
            hits[qid].append(hit)
    return hits


def find_amplicons(primers, hits, genome_dict,
                   min_len, max_len, require_3p,
                   min_fail_len_frac, min_fail_pident):
    results, failed_hits = [], []
    amp_counter, fail_counter = 1, 1

    for pairname, pair in primers.items():
        fwd_seq, rev_seq = pair.get("F"), pair.get("R")
        f_hits, r_hits = hits.get(pairname+"F", []), hits.get(pairname+"R", [])

        for F in f_hits:
            if F["fail_reason"]:
                if F["length"] < len(fwd_seq) * min_fail_len_frac:
                    continue
                if F["pident"] < min_fail_pident:
                    continue
                F["label"] = f"Fail{fail_counter}"
                failed_hits.append(("F", pairname, F))
                fail_counter += 1
                continue
            if not three_prime_match(fwd_seq, F["sseq"], F["strand"], require_3p):
                F["label"] = f"Fail{fail_counter}"
                F["fail_reason"] = f"3prime_mismatch (require {require_3p})"
                failed_hits.append(("F", pairname, F))
                fail_counter += 1
                continue

            for R in r_hits:
                if R["fail_reason"]:
                    if R["length"] < len(rev_seq) * min_fail_len_frac:
                        continue
                    if R["pident"] < min_fail_pident:
                        continue
                    R["label"] = f"Fail{fail_counter}"
                    failed_hits.append(("R", pairname, R))
                    fail_counter += 1
                    continue
                if not three_prime_match(rev_seq, R["sseq"], R["strand"], require_3p):
                    R["label"] = f"Fail{fail_counter}"
                    R["fail_reason"] = f"3prime_mismatch (require {require_3p})"
                    failed_hits.append(("R", pairname, R))
                    fail_counter += 1
                    continue

                if F["contig"] != R["contig"]:
                    continue

                # Orientation check
                valid_normal = (F["strand"] == "plus" and R["strand"] == "minus")
                valid_inverted = (F["strand"] == "minus" and R["strand"] == "plus")
                if not (valid_normal or valid_inverted):
                    continue
                orientation = "normal" if valid_normal else "inverted"

                start = min(F["start"], F["end"], R["start"], R["end"])
                end   = max(F["start"], F["end"], R["start"], R["end"])
                length = end - start + 1
                if min_len <= length <= max_len:
                    seq = str(genome_dict[F["contig"]].seq[start-1:end])
                    results.append({
                        "pair": pairname, "contig": F["contig"],
                        "amplicon_start": start, "amplicon_end": end, "amplicon_len": length,
                        "f_mismatches": F["mismatches"], "f_pident": F["pident"],
                        "r_mismatches": R["mismatches"], "r_pident": R["pident"],
                        "sequence": seq,
                        "label": f"Amp{amp_counter}",
                        "orientation": orientation
                    })
                    amp_counter += 1
                else:
                    fail = {"contig": F["contig"], "start": start, "end": end,
                            "strand": ".", "mismatches": ".", "pident": ".",
                            "sseq": ".", "fail_reason": f"amplicon_len={length} out_of_bounds"}
                    fail["label"] = f"Fail{fail_counter}"
                    failed_hits.append(("Pair", pairname, fail))
                    fail_counter += 1
    return results, failed_hits


def plot_hits(results, failed, out_pdf, genome_dict, primers,
              tick_units="Mb", tick_step=1.0):
    contig_lengths = {cid: len(rec.seq) for cid, rec in genome_dict.items()}
    contigs = sorted(contig_lengths.keys(), key=lambda c: contig_lengths[c], reverse=True)
    fig_height = max(6, len(contigs) * 0.5)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    contig_positions = {c: i for i, c in enumerate(contigs)}

    # Primer map
    primer_names = sorted(set(r["pair"] for r in results))
    primer_map = {name: idx+1 for idx, name in enumerate(primer_names)}

    # Chromosome bars
    for contig in contigs:
        y = contig_positions[contig]
        length = contig_lengths[contig]
        ax.plot([0, length], [y, y], color="black", linewidth=2)

    # Amplicons
    for r in results:
        y = contig_positions[r["contig"]] + 0.3
        color = "green" if r["orientation"] == "normal" else "lightgreen"
        ax.plot([r["amplicon_start"], r["amplicon_end"]], [y, y],
                color=color, linewidth=4)
        label_num = primer_map[r["pair"]]
        ax.text((r["amplicon_start"]+r["amplicon_end"])/2, y+0.2,
                str(label_num), color=color, ha="center", fontsize=8)

    # Failed hits
    for t, pair, hit in failed:
        if hit["contig"] not in contig_positions:
            continue
        y = contig_positions[hit["contig"]] - 0.3
        x = (hit["start"]+hit["end"])//2 if str(hit["start"]).isdigit() else 0
        ax.plot(x, y, marker="|", color="red", markersize=8)

    # Axis ticks with unit conversion
    unit_divisors = {"bp": 1, "kb": 1e3, "Mb": 1e6, "Gb": 1e9}
    divisor = unit_divisors[tick_units]
    max_length = max(contig_lengths.values())
    tick_interval = tick_step * divisor
    xticks = list(range(0, int(max_length) + 1, int(tick_interval)))
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{int(x/divisor)}" for x in xticks])
    ax.set_xlabel(f"Genomic position ({tick_units})")

    # Background shading every other interval
    for i in range(0, int(max_length/divisor), int(tick_step*2)):
        ax.axvspan(i*divisor, (i+tick_step)*divisor, facecolor="lightgray", alpha=0.2)

    # Y-axis
    ax.set_yticks(list(contig_positions.values()))
    ax.set_yticklabels(contigs)
    ax.set_ylabel("Chromosome / contig")
    ax.set_title("Primer Matches and Amplicons")

    # Legend
    legend_handles = [
        plt.Line2D([0],[0], color="green", lw=4, label="Normal viable amplicon"),
        plt.Line2D([0],[0], color="lightgreen", lw=4, label="Inverted viable amplicon"),
        plt.Line2D([0],[0], marker="|", color="red", markersize=8, linestyle="None", label="Failed primer hit")
    ]
    for name, idx in primer_map.items():
        legend_handles.append(plt.Line2D([0],[0], color="white", label=f"{idx}. {name}"))
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()


def main():
    ap = argparse.ArgumentParser(
                        description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter
                    )
    ap.add_argument("--primers", required=True)
    ap.add_argument("--blast", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--min_len", type=int, default=80)
    ap.add_argument("--max_len", type=int, default=3000)
    ap.add_argument("--require_3p", type=int, default=3)
    ap.add_argument("--max_mismatches", type=int, default=10)
    ap.add_argument("--min_pident", type=float, default=0.0)
    ap.add_argument("--len_tolerance", type=int, default=5,
                    help="Allowed difference between hit length and primer length")
    ap.add_argument("--min_fail_len_frac", type=float, default=0.8,
                    help="Minimum fraction of primer length required for failed hits to be reported")
    ap.add_argument("--min_fail_pident", type=float, default=70.0,
                    help="Minimum %% identity required for failed hits to be reported")
    ap.add_argument("--tick_units", choices=["bp","kb","Mb","Gb"], default="Mb",
                    help="Units for x-axis ticks (bp, kb, Mb, Gb)")
    ap.add_argument("--tick_step", type=float, default=1.0,
                    help="Spacing between tick marks in chosen units")
    ap.add_argument("--out_prefix", default="results")
    args = ap.parse_args()

    primers = load_primers(args.primers)
    primer_lengths = {pairname+orient: len(seq)
                      for pairname, pair in primers.items()
                      for orient, seq in pair.items()}
    hits = parse_blast(args.blast, args.max_mismatches, args.min_pident,
                       primer_lengths, args.len_tolerance)
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

    results, failed = find_amplicons(primers, hits, genome_dict,
                                     args.min_len, args.max_len, args.require_3p,
                                     args.min_fail_len_frac, args.min_fail_pident)

    with open(args.out_prefix+"_amplicons.tsv", "w") as out:
        header = ["pair","contig","amplicon_start","amplicon_end","amplicon_len",
                  "f_mismatches","f_pident","r_mismatches","r_pident","sequence","label","orientation"]
        out.write("\t".join(header)+"\n")
        for r in results:
            out.write("\t".join(str(r[h]) for h in header)+"\n")

    with open(args.out_prefix+"_failed.tsv", "w") as out:
        out.write("primer_type\tpair\tcontig\tstart\tend\tstrand\tmismatches\tpident\tlabel\tfail_reason\n")
        for t, pair, hit in failed:
            out.write("\t".join(map(str, [
                t, pair, hit["contig"], hit["start"], hit["end"], hit["strand"],
                hit["mismatches"], hit["pident"], hit.get("label",""), hit.get("fail_reason","")
            ])) + "\n")

    plot_hits(results, failed, args.out_prefix+"_plot.pdf", genome_dict, primers,
              args.tick_units, args.tick_step)

    print(f"Wrote {len(results)} amplicons, {len(failed)} filtered failed hits, and plot to {args.out_prefix}_plot.pdf")


if __name__ == "__main__":
    main()
