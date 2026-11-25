"""
gemtools.py - reusable functions for the GEM projects

Author: Francis John Raymundo
"""

from pathlib import Path
import re
import warnings

import cobra
import pandas as pd

__all__ = [
    "load_models", "model_summary",
    "search_metabolites", "reactions_containing", "exchange_audit",
    "split_gpr", "classify_gpr", "gpr_summary",
    "vmh_to_cobra", "load_diet", "apply_medium", "diet_coverage",
    "assert_present", "sbml_gene_labels", "load_gene_labels", "classify_gene_id", "gene_id_census", "decode_bvbrc", "BLAST_COLUMNS", "read_blast", "filter_hits"
]

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

# --------------------------------------------------------------------
# Model loading
# --------------------------------------------------------------------

def load_models(model_dir, file_map=None, quiet=True):
    """Load SBML models into a dict {label: cobra.Model}."""
    model_dir = Path(model_dir)
    if file_map is None:
        file_map = {p.stem: p.name for p in sorted(model_dir.glob("*.xml"))}
    models = {}
    for label, fname in file_map.items():
        path = model_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"{path} does not exist")
        with warnings.catch_warnings():
            if quiet:
                warnings.simplefilter("ignore")
            models[label] = cobra.io.read_sbml_model(str(path))
    return models


def model_summary(models):
    """One row per model: sizes and objective."""
    rows = []
    for label, m in models.items():
        obj = [r.id for r in m.reactions if r.objective_coefficient != 0]
        rows.append({
            "strain": label,
            "reactions": len(m.reactions),
            "metabolites": len(m.metabolites),
            "genes": len(m.genes),
            "exchanges": len(m.exchanges),
            "objective": obj[0] if obj else None,
        })
    return pd.DataFrame(rows).set_index("strain")

# --------------------------------------------------------------------
# Searching
# --------------------------------------------------------------------

def assert_present(models, met_id, strain):
    """Positive control. Raise if a known identifier is absent.

    A substring search that finds nothing returns an empty result, not
    an error. Pair every audit with this so a misspelled search term
    fails loudly instead of looking like a true absence.
    """
    if strain not in models:
        raise KeyError(f"no model labelled {strain!r}")
    if met_id not in {x.id for x in models[strain].metabolites}:
        raise AssertionError(
            f"POSITIVE CONTROL FAILED: {met_id!r} absent from {strain}. "
            "Check the spelling against the model, not the literature."
        )
    return True


def search_metabolites(models, substring):
    """Metabolites whose id or name contains substring (case-insensitive)."""
    s = substring.lower()
    rows = []
    for label, m in models.items():
        for met in m.metabolites:
            if s in met.id.lower() or s in (met.name or "").lower():
                rows.append({
                    "strain": label, "met_id": met.id, "met_name": met.name,
                    "compartment": met.compartment,
                    "n_reactions": len(met.reactions),
                })
    return pd.DataFrame(rows)


def reactions_containing(models, met_id_fragments):
    """Reactions resolved by METABOLITE CONTENT, never by name.

    Reaction names are unreliable. AGORA2 names EX_dopa_e_ "Exchange of
    Dopamine" while EX_34dhphe_e_ carries levodopa. Filtering on names
    or ids silently misclassifies both.
    """
    frags = [f.lower() for f in met_id_fragments]
    rows = []
    for label, m in models.items():
        for rxn in m.reactions:
            mids = [met.id for met in rxn.metabolites]
            if any(any(f in mid.lower() for f in frags) for mid in mids):
                rows.append({
                    "strain": label, "rxn_id": rxn.id, "rxn_name": rxn.name,
                    "equation": rxn.reaction,
                    "lb": rxn.lower_bound, "ub": rxn.upper_bound,
                    "reversible": rxn.reversibility,
                    "gpr": rxn.gene_reaction_rule,
                    "subsystem": rxn.subsystem,
                })
    return pd.DataFrame(rows)


def exchange_audit(models, met_bases, suffix="_e"):
    """Does each extracellular metabolite have an exchange reaction?

    An extracellular metabolite with no exchange cannot cross the system
    boundary, so the model can neither consume nor produce it.
    """
    rows = []
    for label, m in models.items():
        model_mets = {x.id for x in m.metabolites}
        ex_mets = {met.id for r in m.exchanges for met in r.metabolites}
        for base in met_bases:
            mid = f"{base}{suffix}"
            present, has_ex = mid in model_mets, mid in ex_mets
            rows.append({
                "strain": label, "metabolite": mid,
                "in_model": present, "has_exchange": has_ex,
                "status": ("absent" if not present
                           else "OK" if has_ex else "GAP: no exchange"),
            })
    return pd.DataFrame(rows)

# --------------------------------------------------------------------
# Gene-protein-reaction rules
# --------------------------------------------------------------------

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

# --------------------------------------------------------------------
# Media and diets
# --------------------------------------------------------------------

def vmh_to_cobra(ex_id):
    """VMH 'EX_glc_D[e]' -> AGORA2/COBRApy 'EX_glc_D_e_'.

    The TRAILING underscore is correct. AGORA2 SBML writes
    R_EX_glc_D_e_; COBRApy strips R_ and leaves the rest.
    """
    return ex_id.strip().replace("[e]", "_e_")


def load_diet(path, flux_col="Flux in mmol/(human*day)",
              ex_col="Exchange reaction", convert=True, scale=1.0):
    """Parse a VMH diet CSV into {exchange_id: uptake}.

    Reads the EXCHANGE column, not the metabolite column. Ten rows of
    the VMH diets disagree between the two (vitA -> EX_retinol[e]).

    Units are mmol per person per day; model bounds are mmol per gram
    dry weight per hour. `scale` applies a conversion factor. The
    default of 1.0 applies NONE, so record whatever factor you use and
    state it in the Methods.
    """
    df = pd.read_csv(path)
    for col in (flux_col, ex_col):
        if col not in df.columns:
            raise KeyError(f"column {col!r} not in {list(df.columns)}")
    df = df.dropna(subset=[ex_col])
    diet = {}
    for _, r in df.iterrows():
        ex = vmh_to_cobra(str(r[ex_col])) if convert else str(r[ex_col])
        try:
            val = float(r[flux_col])
        except (TypeError, ValueError):
            continue
        diet[ex] = val * scale
    return diet


def apply_medium(model, composition, close_others=True):
    """Set exchange bounds from {exchange_id: uptake}.

    Returns the requested exchanges ABSENT from the model. Always
    inspect it: a silently dropped nutrient gives a wrong answer with
    no error.
    """
    model_rxns = {r.id for r in model.reactions}
    missing = sorted(k for k in composition if k not in model_rxns)
    present = {k: abs(v) for k, v in composition.items() if k in model_rxns}
    if close_others:
        model.medium = present
    else:
        for k, v in present.items():
            model.reactions.get_by_id(k).lower_bound = -abs(v)
    return missing


def diet_coverage(models, diet):
    """What fraction of a diet can each model actually take up?

    Applying one diet dict to several models without this check gives
    each a DIFFERENT effective medium, invalidating any comparison.
    """
    rows = []
    for label, m in models.items():
        model_ex = {r.id for r in m.exchanges}
        present = [e for e in diet if e in model_ex]
        rows.append({
            "strain": label,
            "diet_components": len(diet),
            "present_in_model": len(present),
            "coverage_pct": round(100 * len(present) / len(diet), 1),
            "missing_examples": ", ".join(
                sorted(set(diet) - model_ex)[:5]),
        })
    return pd.DataFrame(rows).set_index("strain")


# --------------------------------------------------------------------
# Gene identifier decoding
# --------------------------------------------------------------------

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


def sbml_gene_labels(path):
    """Map {cobra_gene_id: original fbc:label} by parsing the SBML.

    COBRApy discards fbc:label on load and replaces dots with
    underscores in gene.id. That substitution is NOT reversible: some
    underscores were already present. WP_002356251.1 becomes
    WP_002356251_1, while EFMMH594_RS13230 is unchanged.

    The original label is therefore only available from the file.
    """
    text = Path(path).read_text(encoding="utf-8")
    pattern = r'<fbc:geneProduct\b[^>]*?fbc:id="([^"]+)"[^>]*?fbc:label="([^"]*)"'
    labels = {}
    for fid, lab in re.findall(pattern, text):
        key = fid[2:] if fid.startswith("G_") else fid
        labels[key] = lab
    return labels


def load_gene_labels(model_dir, file_map):
    """Label maps for several models: {strain: {gene_id: label}}.

    Takes the same arguments as load_models, so the two stay in step.
    """
    model_dir = Path(model_dir)
    return {label: sbml_gene_labels(model_dir / fname)
            for label, fname in file_map.items()}


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

# --------------------------------------------------------------------
# BLAST output parsing  (Day 10)
# --------------------------------------------------------------------

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