#!/usr/bin/env python3

"""
Find reference transcripts (and genes) that are:

  (a) orthology_class == "one2one" in TOGA, AND
  (b) have at least one projection with loss_summ_data status in {I, PI, UL},

and summarize:

  - Strict intersection: present (by that definition) in ALL species.
  - Relaxed: present in ≥ MIN_SPECIES_PRESENT species (e.g. 8 of 9).

Then map these transcript sets to unique chicken genes using the isoforms file.

Assumes directory structure:

  base_dir/
    species_list_11.txt              # one species name per line, e.g. drySqu
    <species>/
      toga_project/
        toga_galGal_on_<species>/
          orthology_classification.tsv
          loss_summ_data.tsv
    thaDol/reference/isoforms.tsv   # or another isoforms file mapping gene -> transcript

Adjust CONFIG below as needed.
"""

import os
import csv

# ---------------- CONFIGURATION ----------------

base_dir = "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/annotation/02_toga_galGal"
species_list_file = os.path.join(base_dir, "species_list_11.txt")

# Isoforms file mapping gene -> transcript (two-column TSV)
isoforms_file = os.path.join(
    base_dir,
    "thaDol", "reference", "isoforms.tsv"
)

# Orthology classes to keep (NOW: only one2one)
KEEP_CLASSES = {"one2one"}

# Projection-level TOGA classes to keep:
# I  = intact
# PI = partially intact
# UL = uncertain loss
ALLOWED_STATUSES = {"I", "PI", "UL"}

# Relaxed threshold: minimum #species in which a transcript must be "good"
MIN_SPECIES_PRESENT = 10

# Output files (transcripts)
strict_tx_out = os.path.join(
    base_dir,
    "toga_transcripts_o2o_I_PI_UL_in_all_species_11.lst"
)
relaxed_tx_out = os.path.join(
    base_dir,
    f"toga_transcripts_o2o_I_PI_UL_in_atleast_{MIN_SPECIES_PRESENT}_species_11.lst"
)

# Output files (genes)
strict_gene_out = os.path.join(
    base_dir,
    "toga_genes_o2o_I_PI_UL_in_all_species_11.lst"
)
relaxed_gene_out = os.path.join(
    base_dir,
    f"toga_genes_o2o_I_PI_UL_in_atleast_{MIN_SPECIES_PRESENT}_species_11.lst"
)

# ------------------------------------------------


def read_species_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def load_projection_status(toga_dir):
    """
    Read PROJECTION-level statuses from loss_summ_data.tsv.

    Returns:
        proj_status[q_transcript] = status
    """
    loss_file = os.path.join(toga_dir, "loss_summ_data.tsv")
    if not os.path.exists(loss_file):
        raise FileNotFoundError(f"Missing loss_summ_data.tsv: {loss_file}")

    proj_status = {}
    with open(loss_file) as lf:
        for line in lf:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            level, ident, status = parts[0], parts[1], parts[2]
            if level == "PROJECTION":
                proj_status[ident] = status
    return proj_status


def transcripts_for_species(sp, base_dir, keep_classes, allowed_statuses):
    """
    For a single species, return the set of t_transcript IDs that satisfy:

      - orthology_class ∈ keep_classes (e.g. {"one2one"})
      - AND at least one projection (q_transcript) has status in allowed_statuses
        (e.g. {"I","PI","UL"}).

    This matches: one2one + I/PI/UL filter per species.
    """
    toga_dir = os.path.join(base_dir, sp, "toga_project", f"toga_galGal_on_{sp}")
    proj_status = load_projection_status(toga_dir)

    orth_file = os.path.join(toga_dir, "orthology_classification.tsv")
    if not os.path.exists(orth_file):
        raise FileNotFoundError(f"Missing orthology_classification.tsv for {sp}: {orth_file}")

    good_t = set()

    with open(orth_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = row.get("orthology_class", "")
            if cls not in keep_classes:
                continue

            t_tx = row.get("t_transcript", "")
            q_tx = row.get("q_transcript", "")

            if not t_tx or t_tx in ("None", ""):
                continue
            if not q_tx or q_tx in ("None", ""):
                continue

            status = proj_status.get(q_tx)
            if status not in allowed_statuses:
                continue

            # This reference transcript has at least one projection with
            # orthology_class in KEEP_CLASSES and status in ALLOWED_STATUSES
            good_t.add(t_tx)

    return good_t


def read_isoforms_tx2gene(path):
    """
    Read isoforms file mapping gene -> transcripts, invert to tx -> gene.

    Expects 2-column TSV: gene_id \t transcript_id
    with or without a header line.
    """
    tx2gene = {}
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
        parts = first.split("\t")
        first_is_header = (len(parts) != 2) or (
            not parts[0].replace("_", "").replace(".", "").isalnum()
        )
        if not first_is_header:
            gene, tx = parts
            tx2gene[tx] = gene

        for line in fh:
            line = line.strip()
            if not line:
                continue
            gene, tx = line.split("\t")
            tx2gene[tx] = gene
    return tx2gene


def main():
    species_list = read_species_list(species_list_file)
    n_species = len(species_list)
    print("Species:", species_list)

    # Per-species transcript sets (filtered by class + I/PI/UL)
    species_to_tset = {}

    for sp in species_list:
        print(f"\nProcessing species: {sp}")
        tset = transcripts_for_species(sp, base_dir, KEEP_CLASSES, ALLOWED_STATUSES)
        species_to_tset[sp] = tset
        print(f"{sp}: {len(tset)} reference transcripts with "
              f"orthology_class in {KEEP_CLASSES} and ≥1 projection with status in {ALLOWED_STATUSES}")

    # Strict intersection: present in ALL species
    strict_tx = None
    for sp in species_list:
        if strict_tx is None:
            strict_tx = set(species_to_tset[sp])
        else:
            strict_tx &= species_to_tset[sp]

    if strict_tx is None:
        strict_tx = set()

    print("\n--------------------------------------------------")
    print(f"Transcripts with {KEEP_CLASSES} + {ALLOWED_STATUSES} in ALL {n_species} species: "
          f"{len(strict_tx)}")

    # Relaxed: present in ≥ MIN_SPECIES_PRESENT species
    tx_counts = {}
    for sp, tset in species_to_tset.items():
        for t in tset:
            tx_counts[t] = tx_counts.get(t, 0) + 1

    relaxed_tx = {t for t, c in tx_counts.items() if c >= MIN_SPECIES_PRESENT}

    print(f"Transcripts with {KEEP_CLASSES} + {ALLOWED_STATUSES} in ≥{MIN_SPECIES_PRESENT} "
          f"of {n_species} species: {len(relaxed_tx)}")

    # Identify "worst" species (fewest transcripts by this criterion)
    worst_sp = min(species_to_tset, key=lambda s: len(species_to_tset[s]))
    print(f"Worst-annotated species by this combined criterion appears to be: "
          f"{worst_sp} ({len(species_to_tset[worst_sp])} transcripts)")

    # Write transcript lists
    with open(strict_tx_out, "w") as out:
        for t in sorted(strict_tx):
            out.write(t + "\n")
    print(f"Strict transcript list written to: {strict_tx_out}")

    with open(relaxed_tx_out, "w") as out:
        for t in sorted(relaxed_tx):
            out.write(t + "\n")
    print(f"Relaxed (≥{MIN_SPECIES_PRESENT} species) transcript list written to: {relaxed_tx_out}")

    # ---- Map transcripts to genes and count orthologous genes ----

    print("\nReading isoforms to map transcripts to genes...")
    tx2gene = read_isoforms_tx2gene(isoforms_file)
    print(f"Transcripts with gene mapping: {len(tx2gene)}")

    # Strict genes
    strict_genes = set()
    missing_tx_strict = 0
    for tx in strict_tx:
        g = tx2gene.get(tx)
        if g is None:
            missing_tx_strict += 1
        else:
            strict_genes.add(g)

    print(f"\nStrict set: {len(strict_tx)} transcripts "
          f"-> {len(strict_genes)} unique genes")
    if missing_tx_strict:
        print(f"(Warning: {missing_tx_strict} strict transcripts had no gene mapping)")

    with open(strict_gene_out, "w") as out:
        for g in sorted(strict_genes):
            out.write(g + "\n")
    print(f"Strict gene list written to: {strict_gene_out}")

    # Relaxed genes
    relaxed_genes = set()
    missing_tx_relaxed = 0
    for tx in relaxed_tx:
        g = tx2gene.get(tx)
        if g is None:
            missing_tx_relaxed += 1
        else:
            relaxed_genes.add(g)

    print(f"\nRelaxed set: {len(relaxed_tx)} transcripts "
          f"-> {len(relaxed_genes)} unique genes")
    if missing_tx_relaxed:
        print(f"(Warning: {missing_tx_relaxed} relaxed transcripts had no gene mapping)")

    with open(relaxed_gene_out, "w") as out:
        for g in sorted(relaxed_genes):
            out.write(g + "\n")
    print(f"Relaxed gene list written to: {relaxed_gene_out}")


if __name__ == "__main__":
    main()

