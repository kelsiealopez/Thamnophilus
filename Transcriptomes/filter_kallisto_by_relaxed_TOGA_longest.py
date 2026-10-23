#!/usr/bin/env python3

import sys
import csv

# ---------- CONSTANT CONFIG (same for all species) ----------------

# Relaxed chicken reference transcripts (>= 8 species, one per line)
relaxed_tx_list = (
    "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/"
    "02_toga_galGal/toga_transcripts_o2o_I_PI_UL_in_atleast_8_species.lst"
)

# TOGA galGal-on-thaDol orthology table
orthology_file = (
    "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/"
    "02_toga_galGal/thaDol/toga_project/toga_galGal_on_thaDol/orthology_classification.tsv"
)

# Base directory that contains all per-species subdirs
top_kallisto_dir = (
    "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/"
    "nf_pipeline/galGal_to_thaDol_anno"
)

# ------------------------------------------------------------------


def read_relaxed_transcripts(path):
    """Read 1 chicken transcript ID per line into a set."""
    tx = set()
    with open(path) as f:
        for line in f:
            t = line.strip()
            if t:
                tx.add(t)
    return tx


def build_qtx_map(orthology_path, relaxed_tset):
    """
    From orthology_classification.tsv, build:

      q2info[q_transcript] = {
          't_transcript': <chicken transcript>,
          't_gene': <chicken gene>
      }

    only for rows where:
      - orthology_class == 'one2one'
      - t_transcript ∈ relaxed_tset
      - q_transcript != 'None'
    """
    q2info = {}
    with open(orthology_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("orthology_class", "") != "one2one":
                continue

            t_tx = row.get("t_transcript", "")
            if t_tx not in relaxed_tset:
                continue

            q_tx = row.get("q_transcript", "")
            if not q_tx or q_tx == "None":
                continue

            t_gene = row.get("t_gene", "") or "NA"

            q2info[q_tx] = {
                "t_transcript": t_tx,
                "t_gene": t_gene,
            }

    return q2info


def read_tx_lengths(lengths_path):
    """
    Read kallisto.merged.transcript_lengths.tsv

    Format:
      tx    gene_id    sample1    sample2    ...

    We compute the mean of the length columns (3rd onward) for each tx.
    """
    tx_len = {}
    with open(lengths_path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)  # skip header
        for row in reader:
            if not row:
                continue
            tx_id = row[0]
            # columns 2+ are gene_id, then per-sample lengths
            numeric = []
            for v in row[2:]:
                try:
                    numeric.append(float(v))
                except (ValueError, TypeError):
                    continue
            if not numeric:
                continue
            avg_len = sum(numeric) / len(numeric)
            tx_len[tx_id] = avg_len
    return tx_len


def choose_longest_qtx_per_tgene(q2info, tx_lengths):
    """
    For each chicken gene (t_gene), among all its projected q_transcripts
    that we have length information for, keep the q_transcript with
    the largest length.

    Returns:
        keep_q: set of q_transcripts to retain
    """
    # t_gene -> (best_q_tx, best_length)
    best = {}

    for q_tx, info in q2info.items():
        t_gene = info["t_gene"]
        length = tx_lengths.get(q_tx)
        if length is None:
            continue

        if t_gene not in best or length > best[t_gene][1]:
            best[t_gene] = (q_tx, length)

    keep_q = {q_tx for (q_tx, _) in best.values()}
    return keep_q


def main():
    if len(sys.argv) != 2:
        sys.stderr.write(
            "Usage: python filter_kallisto_by_relaxed_TOGA_longest.py <species_tag>\n"
            "Example: python filter_kallisto_by_relaxed_TOGA_longest.py drySqu_to_thaDol\n"
        )
        sys.exit(1)

    species_tag = sys.argv[1]
    base_dir = f"{top_kallisto_dir}/{species_tag}/output_1/kallisto"

    kallisto_gene_counts = f"{base_dir}/kallisto.merged.gene_counts.tsv"
    kallisto_tx_lengths = f"{base_dir}/kallisto.merged.transcript_lengths.tsv"
    out_file = f"{base_dir}/kallisto.merged.gene_counts.relaxed_TOGA_longest.tsv"

    print(f"Species: {species_tag}")
    print(f"Kallisto dir: {base_dir}")
    print(f"Input gene_counts: {kallisto_gene_counts}")
    print(f"Input transcript_lengths: {kallisto_tx_lengths}")
    print(f"Output: {out_file}\n")

    print("Reading relaxed TOGA chicken transcripts (>= 8 species)...")
    relaxed_t = read_relaxed_transcripts(relaxed_tx_list)
    print(f"  Relaxed transcripts: {len(relaxed_t)}")

    print("\nBuilding q_transcript -> (chicken_transcript, chicken_gene) map...")
    q2info = build_qtx_map(orthology_file, relaxed_t)
    print(f"  Projections (q_transcripts) for relaxed transcripts: {len(q2info)}")

    print("\nReading Kallisto transcript lengths...")
    tx_len = read_tx_lengths(kallisto_tx_lengths)
    print(f"  Transcript lengths read: {len(tx_len)}")

    print("\nSelecting longest projected transcript per chicken gene...")
    keep_q = choose_longest_qtx_per_tgene(q2info, tx_len)
    print(f"  q_transcripts kept as longest per chicken gene: {len(keep_q)}")

    print("\nFiltering Kallisto gene_counts to those longest projections...")
    kept = 0
    total = 0

    with open(kallisto_gene_counts) as inp, open(out_file, "w") as out:
        reader = csv.reader(inp, delimiter="\t")
        header = next(reader)

        # original header: gene_id, gene_name, sample1, ...
        new_header = (
            header[:2]
            + ["chicken_transcript", "chicken_gene"]
            + header[2:]
        )
        out.write("\t".join(new_header) + "\n")

        for row in reader:
            if not row:
                continue
            total += 1
            q_tx = row[0]  # gene_id column is the projected q_transcript ID

            if q_tx not in keep_q:
                continue

            info = q2info.get(q_tx, {})
            chicken_tx = info.get("t_transcript", "NA")
            chicken_gene = info.get("t_gene", "NA")

            new_row = row[:2] + [chicken_tx, chicken_gene] + row[2:]
            out.write("\t".join(new_row) + "\n")
            kept += 1

    print("Done.")
    print(f"  Total kallisto rows read: {total}")
    print(f"  Rows kept (relaxed TOGA galGal projections, longest per gene): {kept}")
    print(f"  Output written to: {out_file}")


if __name__ == "__main__":
    main()
