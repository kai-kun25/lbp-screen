"""Tests on the real AGORA2 models. Skipped if the files are not present."""
from pathlib import Path

import pytest

import lbpscreen as lbp

ROOT = Path(__file__).resolve().parent.parent
MODEL_DIR = ROOT / "data/processed/core_models"
DIETS = ROOT / "data/raw/diets"

FILE_MAP = {
    "E_faecalis_MMH594":    "Enterococcus_faecalis_MMH594_AGORA2_v2.01.xml",
    "E_lenta_DSM2243":      "Eggerthella_lenta_DSM_2243_AGORA2_v2.01.xml",
    "E_coli_Nissle1917":    "Escherichia_coli_Nissle_1917_AGORA2_v2.01.xml",
    "A_muciniphila_BAA835": "Akkermansia_muciniphila_ATCC_BAA_835_AGORA2_v2.01.xml",
}

needs_models = pytest.mark.skipif(
    not all((MODEL_DIR / f).exists() for f in FILE_MAP.values()),
    reason="AGORA2 model files not found in data/processed/core_models",
)


@pytest.fixture(scope="module")
def models():
    return lbp.load_models(MODEL_DIR, FILE_MAP)


@needs_models
def test_sizes(models):
    assert len(models) == 4
    assert len(models["E_faecalis_MMH594"].reactions) == 1558
    assert len(models["E_lenta_DSM2243"].reactions) == 1035


@needs_models
def test_positive_and_negative_controls(models):
    lbp.assert_present(models, "levdopa_c", "E_faecalis_MMH594")
    lbp.assert_present(models, "mtym_c", "E_lenta_DSM2243")
    with pytest.raises(AssertionError):
        lbp.assert_present(models, "levodopa_c", "E_faecalis_MMH594")  # typo


@needs_models
def test_audit(models):
    audit = lbp.reactions_containing(models, ["levdopa", "dpmn", "tym_", "mtym"])
    assert len(audit) == 13


@needs_models
def test_gene_labels():
    ef = lbp.load_gene_labels(MODEL_DIR, FILE_MAP)["E_faecalis_MMH594"]
    assert len(ef) == 1140
    assert ef["1158653_5_620_peg"] == "1158653.5.620.peg"
    assert ef["WP_002356251_1"] == "WP_002356251.1"


@pytest.mark.skipif(not (DIETS / "high-fiber.csv").exists(),
                    reason="VMH diet file not found")
def test_diet():
    diet = lbp.load_diet(DIETS / "high-fiber.csv")
    assert len(diet) == 91
    assert abs(diet["EX_glc_D_e_"] - 120.211) < 1e-3


CURATED = ROOT / "data/processed/curated_models"


@pytest.mark.skipif(not CURATED.exists(), reason="curated models not built yet")
def test_curated_models_grow():
    for path in CURATED.glob("*.json"):              # every saved curated model
        assert lbp.read_model(path).slim_optimize() > 1e-6, path.name
