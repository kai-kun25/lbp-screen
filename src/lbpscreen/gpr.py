"""
lbpscreen.gpr - Gene-protein-reaction rules.

Author: Francis John Raymundo
"""

import re

import pandas as pd

__all__ = ['split_gpr', 'classify_gpr', 'gpr_summary']


def split_gpr(rule):
    """Individual gene identifiers in a GPR rule, logic discarded."""
    if not rule or not rule.strip():
        return []
    cleaned = rule.replace("(", " ").replace(")", " ")
    cleaned = re.sub(r"\b(and|or)\b", " ", cleaned)
    return [t for t in re.findall(r"[A-Za-z0-9_.\-]+", cleaned) if t]


def classify_gpr(rule):
    """Boolean structure of a GPR rule."""
    if not rule or not rule.strip():
        return "none"
    has_and, has_or = " and " in rule, " or " in rule
    if has_and and has_or:
        return "complex + isozymes"
    if has_and:
        return "complex (and)"
    if has_or:
        return "isozymes (or)"
    return "single gene"


def gpr_summary(models):
    """Counts of each GPR type per model. 'none' reactions can never be
    knocked out, so they appear falsely non-essential."""
    rows = []
    for label, m in models.items():
        counts = pd.Series(
            [classify_gpr(r.gene_reaction_rule) for r in m.reactions]
        ).value_counts()
        counts.name = label
        rows.append(counts)
    out = pd.DataFrame(rows).fillna(0).astype(int)
    out["total"] = out.sum(axis=1)
    out["pct_unannotated"] = (
        100 * out.get("none", 0) / out["total"]
    ).round(1)
    return out
