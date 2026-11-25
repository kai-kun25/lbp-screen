"""
lbpscreen.screen - single-strain analyses used for screening.

Author: Francis John Raymundo
"""
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion

from .analysis import objective_reaction, _grows

__all__ = ["growth_table", "flux_range", "essential_genes"]

def growth_table(models):
    """Growth rate and solver status for each model"""
    rows = []
    for label, model in models.items():
        sol = model.optimize()
        rows.append({
            "Strain": label,
            "Objective": objective_reaction(model).id,
            "Growth": sol.objective_value,
            "Status": sol.status
        })
    return pd.DataFrame(rows).set_index("Strain")

def flux_range(model, reaction_ids, fraction=0.9):  # 0.9 as default fraction
    """
    Gives the difference between the lowest and highest flux each reaction can carry 
    while growth stays at >= fraction of its maximum (flux variability analysis).
    """
    rxns = [model.reactions.get_by_id(r) for r in reaction_ids]
    fva = flux_variability_analysis(model, rxns, fraction_of_optimum=fraction)
    fva["fixed"] = (fva["maximum"] - fva["minimum"]).abs() < 1e-6
    return fva

def essential_genes(model, fraction=0.01):      # fraction = 0.01 by default
    """Genes whose single knockout drops growth below fraction of wild-type."""
    wild = model.slim_optimize()
    results = single_gene_deletion(model)
    results["gene"] = [list(ids)[0] for ids in results["ids"]]
    results["wild_type"] = wild
    results["ko_fraction"] = results["growth"] / wild
    results["essential"] = results["ko_fraction"] < fraction   # True means KO kills growth
    return results[["gene", "growth", "wild_type", "ko_fraction", 
                    "essential"]].set_index("gene").sort_values("growth")
