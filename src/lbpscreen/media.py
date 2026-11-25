"""
lbpscreen.media - Diets and growth media.

Author: Francis John Raymundo
"""

import pandas as pd

__all__ = ['vmh_to_cobra', 'load_diet', 'apply_medium', 'diet_coverage', 'anaerobic', 'apply_diet', 'diet_report',
           'essential_diet_components', 'growth_on_diets']


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
            "coverage_pct": (round(100 * len(present) / len(diet), 1) if diet else 0.0),
            "missing_examples": ", ".join(sorted(set(diet) - model_ex)[:5]),
        })
    return pd.DataFrame(rows).set_index("strain")

OXYGEN_EXCHANGES = ("EX_o2_e_", "EX_o2(e)")

def anaerobic(model):
    """Closes oxygen uptake in the model, since the colon lumen is anaerobic. It returns 
    the IDs of the oxygen exchanges that were closed (empty if the model has none)."""
    closed = []
    for ex_id in OXYGEN_EXCHANGES:
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = 0
            closed.append(ex_id)
        return closed

def apply_diet(model, diet, scale=1.0, aerobic=False, o2_uptake=10.0):
    """Set the model's medium using the diet applied. Returns the diet items
    that are NOT in the model.

        A diet defines the WHOLE MEDIUM: apply_medium closes every exchange not 
    listed in the diet. Oxygen is set here for aerobic conditions, as VMH diets
    never lists oxygen.

    scale      : factor multiplied to every uptake. VMH diets are in mmol/person/day;
                 model bounds are mmol/gDW/h. Always record the factor you set.
    aerobic    : False (default) -> keeps oxygen closed (anaerobic condition)
                 True -> opens oxygen uptake (equals 10 as default)"""

    scaled = {k: v * scale for k, v in diet.items()}
    missing = apply_medium(model, scaled)

    for ex_id in OXYGEN_EXCHANGES:
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -o2_uptake if aerobic else 0

    return missing


def diet_report(model, diet, scale=1.0, aerobic=False):
    """One row per diet item represents exchange metabolites that is in the model
    and is ESSENTIAL (removing the metabolite stop growth). Essential exchange
    metabolites uptake by the cell are the nutrients that actually drive the strain's
    growth"""
    from .analysis import _grows
    with model:
        apply_diet(model, diet, scale=scale, aerobic=aerobic)
        base = model.slim_optimize()
        rows = []
        for ex_id in diet:
            in_model = ex_id in model.reactions
            essential = False
            if in_model and _grows(base):
                with model:
                    model.reactions.get_by_id(ex_id).lower_bound = 0
                    essential = not _grows(model.slim_optimize())
            rows.append({
                "exchange": ex_id,
                "in_model": in_model,
                "requested_uptake": diet[ex_id] * scale,
                "essential": essential
            })
    out = pd.DataFrame(rows)
    out.attrs["growth_on_diet"] = base
    return out.sort_values(["in_model", "essential"], ascending=False)

def essential_diet_components(model, diet, scale=1.0, aerobic=False):
    """List of diet items whose removal stops growth."""
    report = diet_report(model, diet, scale=scale, aerobic=aerobic)
    return report[report["essential"]]["exchange"].tolist()

def growth_on_diets(models, diets, scale=1.0, aerobic=False):
    """Growth of every model on every diet. Each model is restored to its default
    bounds every after application of the diet. Rows represent the model strains. 
    Columns represent each kind of diet applied.
    diets: {diet_name: diet_dict}"""
    table = {}
    for diet_name, diet in diets.items():
        col = {}
        for label, model in models.items():
            with model:
                apply_diet(model, diet, scale=scale, aerobic=aerobic)
                col[label] = model.slim_optimize()
        table[diet_name] = col
    return pd.DataFrame(table)
