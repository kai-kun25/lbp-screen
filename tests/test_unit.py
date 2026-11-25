"""Run: pytest, does not need a model for testing"""
import pandas as pd
import pytest

import lbpscreen as lbp


def test_detect_format():
    assert lbp.detect_format("a/b/Strain_AGORA2.xml") == "sbml"
    assert lbp.detect_format("Strain.XML.GZ") == "sbml_gz"
    assert lbp.detect_format("Strain.mat") == "matlab"
    with pytest.raises(ValueError):
        lbp.detect_format("notes.txt")


def test_vmh_to_cobra():
    assert lbp.vmh_to_cobra("EX_glc_D[e]") == "EX_glc_D_e_"


def test_gpr():
    assert lbp.split_gpr("a or b or c") == ["a", "b", "c"]
    assert lbp.classify_gpr("") == "none"
    assert lbp.classify_gpr("(a and b) or c") == "complex + isozymes"


@pytest.mark.parametrize("label, expected", [
    ("WP_002356251.1",      "refseq_protein"),
    ("349741.23.peg.853",   "bvbrc_peg_canonical"),
    ("1158653.5.620.peg",   "bvbrc_peg_feature"),
    ("479437.5.492_46.peg", "bvbrc_peg_reordered"),
    ("EFMMH594_RS13230",    "ncbi_locus_tag"),
])
def test_classify_gene_id(label, expected):
    assert lbp.classify_gene_id(label) == expected


def test_decode_bvbrc():
    d = lbp.decode_bvbrc("1158653.5.620.peg")
    assert d["fig_id"] == "fig|1158653.5.peg.620"
    assert d["ordering"] == "reordered"
    assert lbp.decode_bvbrc("EFMMH594_RS13230") is None


def test_sbml_gene_labels_any_attribute_order(tmp_path):
    xml = tmp_path / "toy.xml"
    xml.write_text(
        '<fbc:geneProduct fbc:id="G_WP_1_1" fbc:label="WP_1.1"/>\n'
        '<fbc:geneProduct fbc:label="X_RS9" fbc:id="G_X_RS9"/>\n',
        encoding="utf-8")
    assert lbp.sbml_gene_labels(xml) == {"WP_1_1": "WP_1.1", "X_RS9": "X_RS9"}


def test_filter_hits():
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
    res = lbp.filter_hits(fake)
    assert res["passes"].tolist() == [True, False, True]
    assert "coverage" in res.loc[1, "failed"]

def test_load_models_auto_labels_keep_version_dots(tmp_path):
    import cobra
    toy = cobra.Model("toy")
    cobra.io.save_json_model(toy, str(tmp_path / "Strain_AGORA2_v2.01.json"))
    (tmp_path / "notes.txt").write_text("not a model")
    models = lbp.load_models(tmp_path)
    assert list(models) == ["Strain_AGORA2_v2.01"]



