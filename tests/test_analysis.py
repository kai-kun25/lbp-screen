"""Tests for lbpscreen.analysis, using tiny hand-made models."""
import pytest                                   # the testing tool
from cobra import Metabolite, Model, Reaction    # building blocks for test models

import lbpscreen as lbp                         # our package

def test_objective_reaction():
    m = Model("tiny")
    m.add_reactions([Reaction("BIO"), Reaction("OTHER")])
    m.objective = "BIO"
    assert lbp.objective_reaction(m).id == "BIO"

def test_objective_reaction_non_raises():
    m = Model("no_objective")
    m.add_reactions([Reaction("R1")])
    with pytest.raises(ValueError):
        lbp.objective_reaction(m)

def test_boundary_census():
    a_e = Metabolite("a_e", compartment="e")             # outside the cell
    bio_c = Metabolite("biomass_c", compartment="c")     # inside the cell
    ex_a = Reaction("EX_a_e_", lower_bound=-10, upper_bound=1000)
    ex_a.add_metabolites({a_e: -1})                      # a normal exchange
    ex_bio = Reaction("EX_biomass_c_", lower_bound=0, upper_bound=1000)
    ex_bio.add_metabolites({bio_c: -1})                  # named EX_ but on a "c" metabolite
    m = Model("tiny")
    m.add_reactions([ex_a, ex_bio])

    census = lbp.boundary_census(m)
    assert census.loc["EX_a_e_", "Reaction type"] == "exchange"
    assert census.loc["EX_biomass_c_", "Reaction type"] == "other"   # COBRApy does NOT call it an exchange
    assert census.loc["EX_a_e_", "Lower bound"] == -10


def _tiny_growing_model():
    """a_e is taken up, turned into biomass_c, and biomass_c leaves. Grows at 10."""
    a_e = Metabolite("a_e", compartment="e")
    bio_c = Metabolite("biomass_c", compartment="c")
    ex_a = Reaction("EX_a_e_", lower_bound=-10, upper_bound=1000)
    ex_a.add_metabolites({a_e: -1})                       # uptake of up to 10 a_e
    bio = Reaction("BIO", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({a_e: -1, bio_c: 1})              # a_e -> biomass
    ex_bio = Reaction("EX_biomass_c_", lower_bound=0, upper_bound=1000)
    ex_bio.add_metabolites({bio_c: -1})                   # biomass leaves
    m = Model("tiny")
    m.add_reactions([ex_a, bio, ex_bio])
    m.objective = "BIO"
    return m


def test_growth_diagnosis_finds_minimise():
    m = _tiny_growing_model()
    m.objective_direction = "min"                         # break it: minimise growth
    diag = lbp.growth_diagnosis(m)
    assert not diag.loc["1. as loaded", "Does it grow?"]          # no growth as loaded
    assert diag.loc["2. maximize", "Does it grow?"]               # experiment 2 fixes it
    assert m.objective_direction == "min"                 # the model itself is unchanged


def test_growth_diagnosis_finds_closed_boundary():
    m = _tiny_growing_model()
    m.reactions.EX_a_e_.lower_bound = 0                   # break it: no food can enter
    diag = lbp.growth_diagnosis(m)
    assert not diag.loc["2. maximize", "Does it grow?"]           # maximising does not help
    assert diag.loc["3. all boundary open", "Does it grow?"]      # opening the boundary does
    assert m.reactions.EX_a_e_.lower_bound == 0           # the model itself is unchanged


def test_precursor_check_all_producible():
    m = _tiny_growing_model()
    pre = lbp.biomass_precursor_check(m)
    assert pre["Producible?"].all()                  # .all() = every value is True
    assert "TEST_DRAIN_a_e" not in m.reactions      # the temporary drain is gone


def test_precursor_check_finds_missing_ingredient():
    m = _tiny_growing_model()
    m.reactions.EX_a_e_.lower_bound = 0             # break it: no food enters
    pre = lbp.biomass_precursor_check(m).set_index("Metabolite")
    assert not pre.loc["a_e", "Producible?"]         # the check must find the missing ingredient


def test_rescue_test():
    m = _tiny_growing_model()
    m.reactions.EX_a_e_.lower_bound = 0             # break it: no food enters
    assert lbp.rescue_test(m, []) == 0              # supplying nothing: still no growth
    assert lbp.rescue_test(m, ["a_e"]) > 0          # supplying the missing food: growth
    assert "TEST_SOURCE_a_e" not in m.reactions     # the temporary source is gone


def test_producing_reactions():
    m = _tiny_growing_model()
    assert lbp.producing_reactions(m, "biomass_c")["Reaction ID"].tolist() == ["BIO"]


def test_apply_curation():
    m = _tiny_growing_model()
    m.reactions.EX_a_e_.lower_bound = 0              # break it
    changes = [{"action": "add_sink", "metabolite": "a_e", "lb": -10, "reason": "test"}]
    log = lbp.apply_curation(m, changes)
    assert m.slim_optimize() > 0                     # the sink repairs growth
    assert len(log) == 1                             # one change, one log row
    with pytest.raises(ValueError):                  # an unknown action must raise
        lbp.apply_curation(m, [{"action": "delete_everything"}])
