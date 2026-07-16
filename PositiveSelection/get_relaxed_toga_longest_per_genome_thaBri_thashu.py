#!/usr/bin/env python3
import sys
import csv

def read_relaxed_transcripts(path):
    tx = set()
    with open(path) as f:
        for line in f:
            t = line.strip()
            if t:
                tx.add(t)
    return tx

def read_orthology(orthology_path, relaxed_tset):
    """
    Returns: t_gene -> set(q_transcripts)
    Only for rows:
      - orthology_class == 'one2one'
      - t_transcript in relaxed_tset
      - q_transcript != 'None'
    """
    tgene_to_qtx = {}
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
            tgene_to_qtx.setdefault(t_gene, set()).add(q_tx)
    return tgene_to_qtx

def read_fasta_lengths(fasta_path):
    """
    Parse FASTA, return: id -> seq_length
    Assumes ID is first token in header line.
    """
    lengths = {}
    with open(fasta_path) as f:
        seq_id = None
        seq_chunks = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
                if seq_id is not None:
                    lengths[seq_id] = sum(len(c) for c in seq_chunks)
                header = line[1:].split()[0]
                seq_id = header
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if seq_id is not None:
            lengths[seq_id] = sum(len(c) for c in seq_chunks)
    return lengths

def choose_longest_per_tgene(tgene_to_qtx, fasta_lengths):
    """
    For each t_gene, pick q_transcript with longest sequence in this genome.
    """
    keep_q = set()
    for t_gene, qset in tgene_to_qtx.items():
        best_q = None
        best_len = -1
        for q in qset:
            L = fasta_lengths.get(q)
            if L is None:
                continue
            if L > best_len:
                best_len = L
                best_q = q
        if best_q is not None:
            keep_q.add(best_q)
    return keep_q

def main():
    if len(sys.argv) != 5:
        sys.stderr.write(
            "Usage: get_relaxed_toga_longest_per_genome_thaBri_thashu.py "
            "<relaxed_chicken_tx.lst> <orthology_classification.tsv> "
            "<nucleotide.fasta> <out_qtranscript_ids.txt>\n"
        )
        sys.exit(1)

    relaxed_path = sys.argv[1]
    orthology_path = sys.argv[2]
    fasta_path = sys.argv[3]
    out_path = sys.argv[4]

    relaxed = read_relaxed_transcripts(relaxed_path)
    tgene_to_qtx = read_orthology(orthology_path, relaxed)
    fasta_lengths = read_fasta_lengths(fasta_path)
    keep_q = choose_longest_per_tgene(tgene_to_qtx, fasta_lengths)

    with open(out_path, "w") as out:
        for q in sorted(keep_q):
            out.write(q + "\n")

if __name__ == "__main__":
    main()
