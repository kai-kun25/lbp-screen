"""
lbpscreen.genes - Decoding AGORA2 gene identifiers.

Author: Francis John Raymundo
"""

import re

import pandas as pd

__all__ = ['classify_gene_id', 'gene_id_census', 'decode_bvbrc']


GENE_ID_PATTERNS = [
    ("refseq_protein",      r"WP_\d+\.\d+"), # WP_002356251.1
    ("bvbrc_peg_canonical", r"\d+\.\d+\.peg\.\d+"), # 349741.23.peg.853
    ("bvbrc_peg_feature",   r"\d+\.\d+\.\d+\.peg"), # 1158653.5.620.peg
    ("bvbrc_peg_reordered", r"\d+\.\d+\.[0-9_]+\.peg"), # 479437.5.492_46.peg
    ("short_peg",           r"g\.\d+\.peg\.\d+"), # g.1738.peg.2293
    ("ncbi_locus_tag",      r"[A-Za-z0-9]+_RS\d+"), # EFMMH594_RS13230
    ("legacy_locus_tag",    r"[A-Za-z][A-Za-z0-9]*_\d+"), # b4025 or ECOLIN_01234
    ("unknown_literal",     r"Unknown"),
]


def classify_gene_id(label):
    """Name the namespace a gene label belongs to.

    AGORA2 gene labels are NOT one namespace. Seven formats occur, and
    they are mixed within a single model. The retrieval route depends
    on which format a gene carries.
    """
    for name, pattern in GENE_ID_PATTERNS:
        if re.fullmatch(pattern, label):
            return name
    return "unclassified"


def gene_id_census(label_maps):
    """Count gene label formats per model.

    Takes the output of load_gene_labels, not the models themselves.
    """
    rows = []
    for strain, labels in label_maps.items():
        counts = pd.Series(
            [classify_gene_id(l) for l in labels.values()]
        ).value_counts()
        counts.name = strain
        rows.append(counts)
    return pd.DataFrame(rows).fillna(0).astype(int)


def decode_bvbrc(label):
    """Resolve an AGORA2 peg label to a BV-BRC feature identifier.

    Handles both orderings:
        349741.23.peg.853   -> fig|349741.23.peg.853   (canonical)
        1158653.5.620.peg   -> fig|1158653.5.peg.620   (reordered)
        479437.5.492_46.peg -> fig|479437.5.peg.492_46 (reordered)

    Returns a dict, or None if the label is not a peg format.

    The reordered case is a HYPOTHESIS. Verify it against BV-BRC
    before relying on it.
    """
    label = label.strip()

    m = re.fullmatch(r"(\d+\.\d+)\.peg\.(\d+)", label)
    if m:
        gid, feat = m.group(1), m.group(2)
        return {"genome_id": gid, "feature": f"peg.{feat}",
                "fig_id": f"fig|{gid}.peg.{feat}", "ordering": "canonical"}

    m = re.fullmatch(r"(\d+\.\d+)\.([0-9_]+)\.peg", label)
    if m:
        gid, feat = m.group(1), m.group(2)
        return {"genome_id": gid, "feature": f"peg.{feat}",
                "fig_id": f"fig|{gid}.peg.{feat}", "ordering": "reordered"}

    return None
