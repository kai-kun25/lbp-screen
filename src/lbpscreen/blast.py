"""
lbpscreen.blast - Parsing and filtering tabular BLAST output.

Author: Francis John Raymundo
"""

import pandas as pd

__all__ = ['BLAST_COLUMNS', 'read_blast', 'filter_hits']


BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovs",
]


def read_blast(path, columns=None):
    """Read tabular BLAST output into a DataFrame.

    Assumes the extended format:
      -outfmt "6 qseqid sseqid pident length mismatch gapopen
               qstart qend sstart send evalue bitscore qlen slen qcovs"

    Plain -outfmt 6 lacks qlen, slen and qcovs, so coverage cannot be
    computed from it.
    """
    cols = columns or BLAST_COLUMNS
    df = pd.read_csv(path, sep="\t", header=None, names=cols)

    # Coverage from the alignment itself, as a cross-check on qcovs.
    # They differ when a pair has several separate alignments.
    df["qcov_hsp"] = (100 * (df["qend"] - df["qstart"] + 1)
                      / df["qlen"]).round(2)
    df["scov_hsp"] = (100 * (df["send"] - df["sstart"] + 1)
                      / df["slen"]).round(2)
    return df


def filter_hits(df, min_pident=30.0, min_qcov=70.0, max_evalue=1e-20,
                min_bitscore=50.0, coverage_col="qcovs"):
    """Apply threshold criteria and label each hit.

    All four criteria must pass. Identity alone is unreliable: a
    12-residue chance alignment reached 50 percent identity in the
    Day 9 calibration, while genuine homologues ranged from 58 to 99
    percent. Coverage is what separates signal from noise.

    Returns the frame with a boolean 'passes' column and a 'failed'
    column naming which criteria failed. Failures are kept rather
    than dropped, so the filtering is auditable.
    """
    checks = {
        "evalue":   df["evalue"]     <= max_evalue,
        "coverage": df[coverage_col] >= min_qcov,
        "identity": df["pident"]     >= min_pident,
        "bitscore": df["bitscore"]   >= min_bitscore,
    }
    out = df.copy()
    out["passes"] = pd.concat(checks.values(), axis=1).all(axis=1)

    reasons = []
    for i in out.index:
        failed = [name for name, ok in checks.items() if not ok[i]]
        reasons.append(", ".join(failed) if failed else "")
    out["failed"] = reasons
    return out
