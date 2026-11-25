"""Tests for lbpscreen.media diet handling, on a tiny model."""
from cobra import Metabolite, Model, Reaction

import lbpscreen as lbp


def gut_model():
    """Eats sugar (ess.) and iron (ess.); oxygen is OPTIONAL. Grows at 10."""
    sug_e = Metabolite("sug_e", compartment="e")
    fe_e = Metabolite("fe2_e", compartment="e")
    o2_e = Metabolite("o2_e", compartment="e")
    bio_c = Metabolite("biomass_c", compartment="c")
    m = Model("gut")

    def rxn(rid, stoich, lb=0, ub=1000):
        r = Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(stoich)
        return r

    m.add_reactions([
        rxn("EX_sug_e_", {sug_e: -1}, lb=-10),
        rxn("EX_fe2_e_", {fe_e: -1}, lb=-10),
        rxn("EX_o2_e_", {o2_e: -1}, lb=-10),          # oxygen available but not required
        rxn("BIO", {sug_e: -1, fe_e: -1, bio_c: 1}),  # needs sugar AND iron
        rxn("EX_biomass_c_", {bio_c: -1}),
    ])
    m.objective = "BIO"
    return m


DIET = {"EX_sug_e_": 10.0, "EX_fe2_e_": 10.0, "EX_notinmodel_e_": 5.0}


def test_anaerobic_closes_oxygen():
    m = gut_model()
    assert lbp.anaerobic(m) == ["EX_o2_e_"]
    assert m.reactions.EX_o2_e_.lower_bound == 0


def test_apply_diet_reports_missing_and_closes_oxygen():
    m = gut_model()
    missing = lbp.apply_diet(m, DIET)
    assert missing == ["EX_notinmodel_e_"]              # flagged, not silently dropped
    assert m.reactions.EX_o2_e_.lower_bound == 0        # anaerobic by default
    assert m.slim_optimize() > 0


def test_apply_diet_aerobic_opens_oxygen():
    m = gut_model()
    lbp.apply_diet(m, DIET, aerobic=True, o2_uptake=5.0)
    assert m.reactions.EX_o2_e_.lower_bound == -5.0


def test_apply_diet_scale():
    m = gut_model()
    lbp.apply_diet(m, DIET, scale=0.5)
    assert m.reactions.EX_sug_e_.lower_bound == -5.0    # 10 * 0.5


def test_diet_report_flags_essentials():
    m = gut_model()
    report = lbp.diet_report(m, DIET).set_index("exchange")
    assert report.loc["EX_sug_e_", "essential"]
    assert report.loc["EX_fe2_e_", "essential"]
    assert not report.loc["EX_notinmodel_e_", "in_model"]
    assert set(lbp.essential_diet_components(m, DIET)) == {"EX_sug_e_", "EX_fe2_e_"}


def test_growth_on_diets():
    m = gut_model()
    table = lbp.growth_on_diets(
        {"gut": m},
        {"rich": {"EX_sug_e_": 10.0, "EX_fe2_e_": 10.0},
         "no_iron": {"EX_sug_e_": 10.0}})
    assert table.loc["gut", "rich"] > 0
    assert table.loc["gut", "no_iron"] == 0
    assert m.reactions.EX_sug_e_.lower_bound == -10      # original restored


def test_diet_coverage_empty_diet_no_crash():
    cov = lbp.diet_coverage({"gut": gut_model()}, {})    # the zero-division fix
    assert cov.loc["gut", "coverage_pct"] == 0.0
