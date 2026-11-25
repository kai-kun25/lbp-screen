"""
lbpscreen.analysis - growth diagnostics and core flux analyses.

Author: Francis John Raymundo
"""
import math

import pandas as pd
from cobra import Reaction

__all__ = ["objective_reaction", "boundary_census", "growth_diagnosis", "biomass_precursor_check", "rescue_test",
           "producing_reactions", "apply_curation"]

TOL = 1e-6    # Tolerance: fluxes smaller than 1e-6 count as zero

def _grows(value):
    """True if value is a real growth rate above zero."""
    return not math.isnan(value) and value > TOL

def _add_drain(model, met):
    """
    Add a temporary drain reaction that removes the objective function (biomass) reactant from the model.
    """
    drain = Reaction("TEST_DRAIN_" + met.id, lower_bound=0, upper_bound=1000)
    drain.add_metabolites({met: -1})
    model.add_reactions([drain])
    return drain

def _add_source(model, met):
    """
    Add a temporary reaction to the model that supplies
    the metabolite for free.
    """
    source = Reaction("TEST_SOURCE_" + met.id, lower_bound=0, upper_bound=1000)
    source.add_metabolites({met: 1})
    model.add_reactions([source])
    return source

def objective_reaction(model):
    """
    Return the one reaction the model is told to maximize (usually biomass)
    """
    found = [r for r in model.reactions if r.objective_coefficient != 0]
    if len(found) != 1:
        raise ValueError(f"Expected one (1) objective reaction, found {len(found)}")
    return found[0]

def boundary_census(model):
    """
    Table of every boundary reaction, the type COBRApy gives it, and its bounds.
    """
    exchange_ids = {r.id for r in model.exchanges}      # List of exchange reactions
    demand_ids = {r.id for r in model.demands}          # List of demand reaction
    sink_ids = {r.id for r in model.sinks}              # List of sink reactions

    rows = []
    # Evaluates if each reaction is in the list of reaction type given above
    for r in model.boundary:
        met = list(r.metabolites)[0]
        if r.id in exchange_ids:
            kind = "exchange"
        elif r.id in demand_ids:
            kind = "demand"
        elif r.id in sink_ids:
            kind = "sink"
        else:
            kind = "other"
        rows.append({
            "Reaction ID": r.id,
            "Reaction type": kind,
            "Metabolite": met.id,
            "Compartment": met.compartment,
            "Lower bound": r.lower_bound,
            "Upper bound": r.upper_bound,
        })
    return pd.DataFrame(rows).set_index("Reaction ID")

def growth_diagnosis(model):
    """
    Growth under experiments that each remove one possible cause of zero growth.
    """
    rows = []

    # Case 1: If the as loaded model does not grow
    g = model.slim_optimize()
    rows.append({"Experiment": "1. as loaded", "Growth": g, "Does it grow?": _grows(g)})

    # Case 2: If the objective direction was set to minimum (there is growth when direction is maximize)
    with model:
        model.objective_direction = "max"
        g = model.slim_optimize()
    rows.append({"Experiment": "2. maximize", "Growth": g, "Does it grow?": _grows(g)})

    # Case 3: If there is a problem with bounds of boundary reactions
    with model:
        for r in model.boundary:
            r.bounds = (-1000, 1000)
        g = model.slim_optimize()
    rows.append({"Experiment": "3. all boundary open", "Growth": g, "Does it grow?": _grows(g)})

    # Case 4: If Case 2 and 3 are both issues
    with model:
        model.objective_direction = "max"
        for r in model.boundary:
            r.bounds = (-1000, 1000)
        g = model.slim_optimize()
    rows.append({"Experiment": "4. both 2 and 3", "Growth": g, "Does it grow?": _grows(g)})

    return pd.DataFrame(rows).set_index("Experiment")

def biomass_precursor_check(model):
    """
    Check if the model can make each biomass reactant (i.e., has an alternative pathway
    that produces the biomass reactant. Reactants that cannot be produced are listed first
    """
    bio = objective_reaction(model)                 # Reaction of the objective function
    rows = []
    for met in bio.reactants:
        # Temporary addition of the sink reaction of met to the model
        with model:
            drain = _add_drain(model, met)
            model.objective = drain
            model.objective_direction = "max"
            v = model.slim_optimize()
        rows.append({
            "Metabolite": met.id,
            "Name": met.name,
            "Coefficient": bio.metabolites[met],
            "Maximum production": v,
            "Producible?": _grows(v),
        })
    return pd.DataFrame(rows).sort_values("Producible?")    # False before True, so problematic metabolites comes first

def rescue_test(model, met_ids):
    """
    Check the growth when the listed metabolites are supplied for free.
    Growth means that the metabolite indeed blocks the biomass reaction.
    """
    with model:
        for met_id in met_ids:          # met_ids contains possible metabolite blockers
            met = model.metabolites.get_by_id(met_id)
            _add_source(model, met)
        g = model.slim_optimize()
    return g

def producing_reactions(model, met_id):
    """
    Lists the reactions that can make a metabolite, given their current bounds
    """
    met = model.metabolites.get_by_id(met_id)      # met_id --> retrieves metabolite
    rows = []
    for r in met.reactions:
        coef = r.metabolites[met]                  # Retrieves the coefficient of metabolite in a reaction
        makes_forward = coef > 0 and r.upper_bound > 0
        makes_backward = coef < 0 and r.lower_bound < 0
        if makes_forward or makes_backward:
            rows.append({
                "Reaction ID": r.id,
                "Name": r.name,
                "Equation": r.reaction,
                "Lower bound": r.lower_bound,
                "Upper bound": r.upper_bound,
            })
    return pd.DataFrame(rows)

def apply_curation(model, changes):
    """
    Apply a list of changes to model (in place). Returns a table logging each change.
    If action is "add_sink" -> provide metabolite and lower bound
    If action is "set_bounds" -> provide reaction ID, and the upper + lower bounds
    """
    for change in changes:
        if change["action"] == "add_sink":        # If curation is adding a sink reaction
            met = model.metabolites.get_by_id(change["metabolite"])
            sink = Reaction("sink_" + met.id + "_curated", lower_bound=change["lb"], upper_bound=1000)
            sink.add_metabolites({met: -1})
            model.add_reactions([sink])
        elif change["action"] == "set_bounds":    # If curation is setting bounds
            rxn = model.reactions.get_by_id(change["reaction"])
            rxn.bounds = (change["lb"], change["ub"])
        else:
            raise ValueError(f"Unknown action: {change['action']}")
    return pd.DataFrame(changes)
