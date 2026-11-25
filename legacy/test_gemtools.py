"""Smoke test for gemtools. Run: python test_gemtools.py"""
import sys
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import gemtools as gt

MODEL_DIR = Path(__file__).parent.parent / "data/processed/core_models"
DIETS     = Path(__file__).parent.parent / "data/raw/diets"

FILE_MAP = {
    "E_faecalis_MMH594":    "Enterococcus_faecalis_MMH594_AGORA2_v2.01.xml",
    "E_lenta_DSM2243":      "Eggerthella_lenta_DSM_2243_AGORA2_v2.01.xml",
    "E_coli_Nissle1917":    "Escherichia_coli_Nissle_1917_AGORA2_v2.01.xml",
    "A_muciniphila_BAA835": "Akkermansia_muciniphila_ATCC_BAA_835_AGORA2_v2.01.xml",
}

def main():
    models = gt.load_models(MODEL_DIR, FILE_MAP)

    assert len(models) == 4
    assert len(models["E_faecalis_MMH594"].reactions) == 1558
    assert len(models["E_lenta_DSM2243"].reactions) == 1035

    gt.assert_present(models, "levdopa_c", "E_faecalis_MMH594")
    gt.assert_present(models, "mtym_c",    "E_lenta_DSM2243")

    audit = gt.reactions_containing(models, ["levdopa", "dpmn", "tym_", "mtym"])
    assert len(audit) == 13, f"expected 13 audit rows, got {len(audit)}"

    assert gt.vmh_to_cobra("EX_glc_D[e]") == "EX_glc_D_e_"
    assert gt.split_gpr("a or b or c") == ["a", "b", "c"]
    assert gt.classify_gpr("") == "none"

    diet = gt.load_diet(DIETS / "high-fiber.csv")
    assert len(diet) == 91
    assert abs(diet["EX_glc_D_e_"] - 120.211) < 1e-3

    assert gt.classify_gene_id("WP_002356251.1") == "refseq_protein"
    assert gt.classify_gene_id("349741.23.peg.853") == "bvbrc_peg_canonical"
    assert gt.classify_gene_id("1158653.5.620.peg") == "bvbrc_peg_feature"
    assert gt.classify_gene_id("479437.5.492_46.peg") == "bvbrc_peg_reordered"
    assert gt.classify_gene_id("EFMMH594_RS13230") == "ncbi_locus_tag"

    d = gt.decode_bvbrc("1158653.5.620.peg")
    assert d["fig_id"] == "fig|1158653.5.peg.620"
    assert d["ordering"] == "reordered"
    assert gt.decode_bvbrc("EFMMH594_RS13230") is None

    labels = gt.load_gene_labels(MODEL_DIR, FILE_MAP)
    ef = labels["E_faecalis_MMH594"]
    assert len(ef) == 1140
    assert ef["1158653_5_620_peg"] == "1158653.5.620.peg"
    assert ef["WP_002356251_1"]    == "WP_002356251.1"
    assert ef["EFMMH594_RS13230"]  == "EFMMH594_RS13230"

    fake = pd.DataFrame({
        "qseqid": ["real_full", "noise_short", "real_divergent"],
        "sseqid": ["x", "x", "x"],
        "pident": [99.9, 50.0, 58.1],
        "length": [1017, 12, 289],
        "mismatch": [1, 6, 121], "gapopen": [0, 0, 0],
        "qstart": [1, 51, 1], "qend": [1017, 62, 289],
        "sstart": [1, 487, 1], "send": [1017, 498, 289],
        "evalue": [0.0, 3.1, 2e-108],
        "bitscore": [2117.0, 16.2, 326.0],
        "qlen": [1017, 1017, 289], "slen": [1017, 1017, 289],
        "qcovs": [100, 1, 100],
    })
    res = gt.filter_hits(fake)
    assert bool(res.loc[0, "passes"])          # full-length match
    assert not bool(res.loc[1, "passes"])      # short chance alignment
    assert "coverage" in res.loc[1, "failed"]
    assert bool(res.loc[2, "passes"])          # divergent but real

    # The negative control. A misspelled identifier must RAISE.
    try:
        gt.assert_present(models, "levodopa_c", "E_faecalis_MMH594")
    except AssertionError:
        pass
    else:
        raise AssertionError("negative control failed: typo was not caught")

    print("All tests passed.")


if __name__ == "__main__":
    main()