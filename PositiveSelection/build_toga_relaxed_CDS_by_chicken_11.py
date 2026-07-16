#!/usr/bin/env python3
import sys, os, csv

def read_cds_fasta(path):
    """
    Return dict: transcript_id_without_protein_suffix -> CDS sequence

    TransDecoder CDS headers look like:
      >XM_003641458.6.385.p1 GENE.XM_003641458.6.385~~XM_003641458.6.385.p1 ...

    We want to index by 'XM_003641458.6.385' so that it matches the q IDs
    in relaxed_toga_orthologs_9sp.tsv (which do NOT have the .p1 suffix).
    """
    d = {}
    with open(path) as f:
        name = None
        buf = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    # use first token, drop .p1/.p2/... suffix if present
                    base = name.split()[0]
                    if ".p" in base:
                        base = base.rsplit(".p", 1)[0]
                    d[base] = "".join(buf)
                # new header
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            base = name.split()[0]
            if ".p" in base:
                base = base.rsplit(".p", 1)[0]
            d[base] = "".join(buf)
    return d

def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: build_toga_relaxed_CDS_by_chicken_11.py "
            "<relaxed_toga_orthologs_11sp.tsv> <out_dir>\n"
        )
        sys.exit(1)

    ortho_tsv = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    POS_BASE = "/n/netscratch/edwards_lab/Lab/kelsielopez/Thamnophilus/positive_sel"

    # adjust names if needed
    species_cds = {
        "drySqu": os.path.join(POS_BASE, "drySqu_toga_relaxed_TD_11", "drySqu_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "sakCan": os.path.join(POS_BASE, "sakCan_toga_relaxed_TD_11", "sakCan_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "sakCri": os.path.join(POS_BASE, "sakCri_toga_relaxed_TD_11", "sakCri_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "sakLuc": os.path.join(POS_BASE, "sakLuc_toga_relaxed_TD_11", "sakLuc_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "thaAtr": os.path.join(POS_BASE, "thaAtr_toga_relaxed_TD_11", "thaAtr_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "thaBer": os.path.join(POS_BASE, "thaBer_toga_relaxed_TD_11", "thaBer_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "thaCae": os.path.join(POS_BASE, "thaCae_toga_relaxed_TD_11", "thaCae_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "thaDol": os.path.join(POS_BASE, "thaDol_toga_relaxed_TD_11", "thaDol_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "thaRuf": os.path.join(POS_BASE, "thaRuf_toga_relaxed_TD_11", "thaRuf_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "thaBri": os.path.join(POS_BASE, "thaBri_toga_relaxed_TD_11", "thaBri_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
        "thaShu": os.path.join(POS_BASE, "thaShu_toga_relaxed_TD_11", "thaShu_relaxed_longest_transcripts_11.fa.transdecoder.cds"),
    }

    cds_by_sp = {}
    for sp, path in species_cds.items():
        sys.stderr.write(f"Loading CDS for {sp} from {path}\n")
        cds_by_sp[sp] = read_cds_fasta(path)

    with open(ortho_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        species_cols = [c for c in reader.fieldnames if c.endswith("_q")]
        species = [c.replace("_q", "") for c in species_cols]

        written = 0

        for row in reader:
            t_gene = row["t_gene"]
            t_tx  = row["t_transcript"]

            seqs = []
            for sp, col in zip(species, species_cols):
                q = row[col].strip()
                if not q:
                    continue
                seq = cds_by_sp[sp].get(q)
                if not seq:
                    continue
                seqs.append((sp, seq))

            # require at least 3 species with CDS
            if len(seqs) < 3:
                continue

            # filename based on chicken gene and transcript
            base = f"{t_gene}__{t_tx}".replace("/", "_")
            out_path = os.path.join(out_dir, f"{base}.cds.fa")

            with open(out_path, "w") as out_fh:
                for sp, s in seqs:
                    out_fh.write(f">{sp}\n")
                    for i in range(0, len(s), 60):
                        out_fh.write(s[i:i+60] + "\n")

            written += 1
            if written % 1000 == 0:
                sys.stderr.write(f"Wrote {written} CDS files so far...\n")

    sys.stderr.write(f"Total CDS files written: {written}\n")

if __name__ == "__main__":
    main()
