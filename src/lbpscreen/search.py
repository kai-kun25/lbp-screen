"""
lbpscreen.search - Searching models for metabolites and reactions.

Author: Francis John Raymundo
"""

import pandas as pd

__all__ = ['assert_present', 'search_metabolites', 'reactions_containing', 'exchange_audit']


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
