"""Tests for lbpscreen.screen, on a tiny hand-made model."""
from cobra import Metabolite, Model, Reaction

import lbpscreen as lbp


def parallel_model():
    """a enters, becomes b by EITHER AB or AB2, then b -> biomass. Grows at 10."""
    a_e = Metabolite("a_e", compartment="e")
    b_c = Metabolite("b_c", compartment="c")
    bio_c = Metabolite("biomass_c", compartment="c")
    m = Model("parallel")

    def rxn(rid, stoich, lb=0, ub=1000, gpr=""):
        r = Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(stoich)
        r.gene_reaction_rule = gpr
        return r

    m.add_reactions([
        rxn("EX_a_e_", {a_e: -1}, lb=-10),
        rxn("AB",  {a_e: -1, b_c: 1}, gpr="geneA"),    # one of two routes in
        rxn("AB2", {a_e: -1, b_c: 1}, gpr="geneA2"),   # the backup route
        rxn("BIO", {b_c: -1, bio_c: 1}, gpr="geneB"),  # only biomass route: essential
        rxn("EX_biomass_c_", {bio_c: -1}),
    ])
    m.objective = "BIO"
    return m


def test_growth_table():
    table = lbp.growth_table({"x": parallel_model()})
    assert table.loc["x", "Status"] == "optimal"
    assert abs(table.loc["x", "Growth"] - 10) < 1e-6


def test_flux_range():
    m = parallel_model()
    fva = lbp.flux_range(m, ["BIO", "AB"], fraction=1.0)
    assert abs(fva.loc["BIO", "minimum"] - 10) < 1e-6   # biomass flux is forced to 10
    assert bool(fva.loc["BIO", "fixed"])                # no freedom
    assert not bool(fva.loc["AB", "fixed"])             # AB varies 0..10 (AB2 covers rest)


def test_essential_genes():
    ess = lbp.essential_genes(parallel_model())
    assert bool(ess.loc["geneB", "essential"])          # only biomass route
    assert not bool(ess.loc["geneA", "essential"])      # AB2 is a backup
    assert not bool(ess.loc["geneA2", "essential"])