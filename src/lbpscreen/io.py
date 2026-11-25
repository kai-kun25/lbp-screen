"""
lbpscreen.io - reading AGORA2 models in any of their distributed formats.

Author: Francis John Raymundo
"""

from pathlib import Path
import gzip
import re
import warnings

import cobra
import pandas as pd

__all__ = ["SUPPORTED_FORMATS", "detect_format", "sbml_safe_id",
           "harmonize_ids", "read_model", "load_models",
           "model_summary", "sbml_gene_labels", "load_gene_labels"]

# File ending -> format name. Longer endings first so ".xml.gz" wins over ".gz".
SUPPORTED_FORMATS = {
    ".xml.gz": "sbml_gz",
    ".sbml.gz": "sbml_gz",
    ".xml": "sbml",
    ".sbml": "sbml",
    ".mat": "matlab",
    ".json": "json",
    ".yml": "yaml",
    ".yaml": "yaml",
}


def detect_format(path):
    """Return the format name for a model file, judged by its ending."""
    name = Path(path).name.lower()
    for ending, fmt in SUPPORTED_FORMATS.items():
        if name.endswith(ending):
            return fmt
    raise ValueError(
        f"Unrecognised model file {Path(path).name!r}. "
        f"Supported endings: {', '.join(SUPPORTED_FORMATS)}"
    )


def _strip_ending(filename):
    """File name without its model-format ending.

    "X_AGORA2_v2.01.xml" -> "X_AGORA2_v2.01". Only the matched ending is
    removed, so dots inside the name (version numbers) survive.
    """
    lower = filename.lower()
    for ending in SUPPORTED_FORMATS:
        if lower.endswith(ending):
            return filename[: -len(ending)]
    raise ValueError(f"Unrecognised model file {filename!r}")


def sbml_safe_id(identifier):
    """Spell an identifier the way an AGORA2 SBML file does.

    SBML identifiers may contain only letters, digits and underscores,
    so every other character becomes "_":
        EX_glc_D[e]      -> EX_glc_D_e_
        DM_dhptd(c)      -> DM_dhptd_c_
        316435.3.115.peg -> 316435_3_115_peg
    """
    return re.sub(r"[^A-Za-z0-9_]", "_", identifier)


def harmonize_ids(model):
    """Rename reactions, metabolites and genes to SBML-style ids, in place.

    This lets a JSON or MATLAB model be used with code (diets, tests)
    written for the SBML spelling. The original id of every renamed
    object is kept in obj.notes["original_id"], so nothing is lost.
    Returns the number of objects renamed.
    """
    renamed = 0
    for attr in ("reactions", "metabolites"):
        objs = getattr(model, attr)
        existing = {o.id for o in objs}
        for obj in list(objs):
            new = sbml_safe_id(obj.id)
            if new == obj.id:
                continue
            if new in existing:
                raise ValueError(f"renaming {obj.id!r} -> {new!r} would "
                                 f"collide with an existing {attr[:-1]}")
            obj.notes["original_id"] = obj.id
            existing.discard(obj.id)
            existing.add(new)
            obj.id = new
            renamed += 1

    gene_map = {g.id: sbml_safe_id(g.id) for g in model.genes
                if sbml_safe_id(g.id) != g.id}
    if len(set(gene_map.values())) != len(gene_map):
        raise ValueError("two genes would receive the same SBML-style id")
    if gene_map:
        originals = dict(gene_map)
        # rename_genes also rewrites every GPR rule that mentions them
        cobra.manipulation.rename_genes(model, gene_map)
        for old, new in originals.items():
            model.genes.get_by_id(new).notes["original_id"] = old
        renamed += len(gene_map)
    return renamed


def read_model(path, quiet=True, harmonize=True):
    """Read ONE model file, whatever its format, into a cobra.Model.

    SBML is the reference format for this package. JSON and MATLAB
    files spell some ids differently (EX_glc_D[e] vs EX_glc_D_e_;
    WP_1.1 vs WP_1_1; checked on Day 2). With harmonize=True (default)
    they are renamed to the SBML spelling, so every format gives the
    same ids. For SBML, the original gene names are attached to
    gene.notes["original_id"].
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")
    fmt = detect_format(path)

    with warnings.catch_warnings():
        if quiet:
            warnings.simplefilter("ignore")
        if fmt == "sbml":
            model = cobra.io.read_sbml_model(str(path))
        elif fmt == "sbml_gz":
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                model = cobra.io.read_sbml_model(handle.read())
        elif fmt == "matlab":
            try:
                model = cobra.io.load_matlab_model(str(path))
            except NotImplementedError as err:
                # scipy cannot read MATLAB v7.3 (HDF5) files.
                raise ValueError(
                    f"{path.name} looks like a MATLAB v7.3 file, which "
                    "scipy cannot read. Use the SBML version instead."
                ) from err
        elif fmt == "json":
            model = cobra.io.load_json_model(str(path))
        else:  # yaml
            model = cobra.io.load_yaml_model(str(path))

    if fmt in ("sbml", "sbml_gz"):
        # SBML ids are already SBML-style. Attach the true gene names,
        # which COBRApy dropped, so every format carries them.
        for gid, label in sbml_gene_labels(path).items():
            if gid in model.genes:
                model.genes.get_by_id(gid).notes["original_id"] = label
    elif harmonize:
        model.notes["ids_renamed"] = harmonize_ids(model)

    # Record where the model came from; useful in every later report.
    model.notes["source_file"] = path.name
    model.notes["source_format"] = fmt
    return model


def load_models(model_dir, file_map=None, quiet=True):
    """Load several models into a dict {label: cobra.Model}.

    file_map : {label: filename}. If None, every supported model file in
               model_dir is loaded and labelled by its file name.
    """
    model_dir = Path(model_dir)
    if file_map is None:
        file_map = {}
        for p in sorted(model_dir.iterdir()):
            try:
                detect_format(p)
            except ValueError:
                continue  # not a model file; skip it
            label = _strip_ending(p.name)
            if label in file_map:
                raise ValueError(
                    f"Two files share the label {label!r} (e.g. .xml and "
                    ".mat of one strain). Pass an explicit file_map."
                )
            file_map[label] = p.name
    return {label: read_model(model_dir / fname, quiet=quiet)
            for label, fname in file_map.items()}


def model_summary(models):
    """One row per model: sizes, objective and source format."""
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
            "format": m.notes.get("source_format"),
        })
    return pd.DataFrame(rows).set_index("strain")


def sbml_gene_labels(path):
    """Map {cobra_gene_id: original fbc:label} by reading the SBML text.

    SBML identifiers may not contain dots, so the file itself stores
    WP_002356251.1 as fbc:id="G_WP_002356251_1" and keeps the true name
    only in fbc:label. COBRApy keeps the id and drops the label. The
    dot -> underscore change is NOT reversible (EFMMH594_RS13230 already
    had an underscore), so the label must be read from the file.

    Attributes are found one at a time, so their order inside the tag
    does not matter.
    """
    path = Path(path)
    if detect_format(path) == "sbml_gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            text = handle.read()
    else:
        text = path.read_text(encoding="utf-8")

    labels = {}
    for tag in re.findall(r"<fbc:geneProduct\b[^>]*>", text):
        fid = re.search(r'\bfbc:id="([^"]+)"', tag)
        lab = re.search(r'\bfbc:label="([^"]*)"', tag)
        if fid and lab:
            key = fid.group(1)
            key = key[2:] if key.startswith("G_") else key
            labels[key] = lab.group(1)
    return labels


def load_gene_labels(model_dir, file_map):
    """Label maps for several SBML models: {strain: {gene_id: label}}.

    Takes the same arguments as load_models, so the two stay in step.
    Non-SBML files are skipped: .mat files keep the original gene names.
    """
    model_dir = Path(model_dir)
    out = {}
    for label, fname in file_map.items():
        if detect_format(model_dir / fname) in ("sbml", "sbml_gz"):
            out[label] = sbml_gene_labels(model_dir / fname)
    return out
