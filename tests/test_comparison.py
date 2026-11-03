import numpy as np
import pandas as pd
import pytest
from scipy import stats

from imagestream_analysis.comparison import COMPARISON_FILE, MEMBRANE, SUBSETS, bh, compare, write_comparison
from imagestream_analysis.config import Pool
from imagestream_analysis.project import DatasetEntry, Project, ZipEntry

NAMES = ["BF", "MOI", "Immune", "Tumor", "A"]
SUBSET_COLUMNS = ["cell type", "neighbour type", "contacts", "expected (subtypes shuffled)",
                  "observed / expected (subtypes shuffled)", "95% low", "95% high"]


def project(tmp_path, datasets, experiment=lambda name: "E") -> tuple[Project, pd.DataFrame]:
    """A project with one zip; datasets: (name, condition, replicate, {(cell type, partner): (contacts, expected)})."""
    p = Project(tmp_path)
    ch = p.settings.channels
    ch.brightfield, ch.marker_of_interest, ch.tumor = "BF", ["MOI"], "Tumor"
    ch.immune, ch.functional = ["Immune"], ["A"]
    entries, folders = [], []
    for name, condition, replicate, contacts in datasets:
        entries.append(DatasetEntry(name, experiment(name), True, condition, replicate))
        folder = tmp_path / experiment(name) / name
        folder.mkdir(parents=True)
        rows = [{"cell type": c, "neighbour type": n, "contacts": o, "expected (subtypes shuffled)": e}
                for (c, n), (o, e) in contacts.items()]
        pd.DataFrame(rows).reindex(columns=SUBSET_COLUMNS).to_csv(folder / "Subset preference.tsv", sep="\t",
                                                                   index=False)
        folders.append(f"{experiment(name)}/{name}")
    p.zips.append(ZipEntry(tmp_path / "export.zip", 0.5, [1, 2, 3, 4, 5], list(NAMES), {}, entries))
    return p, pd.DataFrame({"folder": folders})


A_TUMOR, NEG_TUMOR, A_A = ("Immune A+", "Tumor"), ("Immune negative", "Tumor"), ("Immune A+", "Immune A+")


def subsets(result):
    return result[result["measure"] == SUBSETS].set_index(["cell type", "partner"])


def test_paired_over_replicates(tmp_path):
    ctrl = [(40, 40.0), (50, 48.0), (30, 31.0)]
    treated = [(52, 40.0), (58, 46.0), (40, 32.0)]
    data = []
    for k, (c, t) in enumerate(zip(ctrl, treated, strict=True)):
        data.append((f"D{k}_c", "control", f"D{k}", {A_TUMOR: c, NEG_TUMOR: (20, 20.0), A_A: (2, 2.0)}))
        data.append((f"D{k}_t", "treated", f"D{k}", {A_TUMOR: t, NEG_TUMOR: (18, 19.0), A_A: (3, 2.5)}))
    p, table = project(tmp_path, data)
    result, notes, conditions, reference = compare(p, table)
    assert conditions == ["control", "treated"] and reference == "control"
    assert any("paired over 3 replicates" in n for n in notes)
    row = subsets(result).loc[A_TUMOR]
    y = lambda o, e: np.log((o + 0.5) / (e + 0.5))  # noqa: E731
    d = np.array([y(*t) - y(*c) for c, t in zip(ctrl, treated, strict=True)])
    assert row["effect"] == pytest.approx(np.exp(d.mean()))
    assert row["p"] == pytest.approx(stats.ttest_1samp(d, 0).pvalue)
    assert row["95% low"] < row["effect"] < row["95% high"]
    assert row["design"] == "paired" and row["replicates"] == 3
    assert row["per replicate"].startswith("D0 1.00 → 1.30")
    # under 5 expected contacts: not tested, and why
    small = subsets(result).loc[A_A]
    assert np.isnan(small["p"]) and "too few contacts" in small["note"]
    # q: Benjamini-Hochberg over the two tested pairs
    tested = subsets(result).dropna(subset=["p"])
    assert list(tested["q"]) == pytest.approx(list(bh(list(tested["p"]))))


def test_immune_pairs_once(tmp_path):
    both = {("Immune A+", "Immune negative"): (10, 9.0), ("Immune negative", "Immune A+"): (10, 9.0)}
    p, table = project(tmp_path, [("c", "control", "D1", both), ("t", "treated", "D1", both)])
    result = compare(p, table)[0]
    assert len(result) == 1 and result.iloc[0]["cell type"] == "Immune A+"


def test_unpaired_welch(tmp_path):
    ctrl, treated = [(40, 40.0), (50, 48.0), (30, 31.0)], [(52, 40.0), (58, 46.0)]
    data = [(f"c{k}", "control", f"donor {k}", {A_TUMOR: v}) for k, v in enumerate(ctrl)]
    data += [(f"t{k}", "treated", f"donor {k + 5}", {A_TUMOR: v}) for k, v in enumerate(treated)]
    p, table = project(tmp_path, data)
    result, notes, _, _ = compare(p, table)
    y = lambda o, e: np.log((o + 0.5) / (e + 0.5))  # noqa: E731
    expect = stats.ttest_ind([y(*v) for v in treated], [y(*v) for v in ctrl], equal_var=False)
    row = subsets(result).loc[A_TUMOR]
    assert row["design"] == "unpaired" and row["replicates"] == "3 + 2"
    assert row["p"] == pytest.approx(expect.pvalue)
    assert any("unpaired" in n for n in notes)


def test_no_test_without_replicates(tmp_path):
    c = {A_TUMOR: (40, 40.0)}
    p, table = project(tmp_path, [("c", "control", "D1", c), ("t", "treated", "D1", c), ("x", "treated", "", c)])
    result, notes, _, _ = compare(p, table)
    assert np.isnan(subsets(result).loc[A_TUMOR]["p"])
    assert any("two or more replicates" in n for n in notes)
    assert any("Without a replicate, left out: x" in n for n in notes)


def test_technical_datasets_are_pooled(tmp_path):
    data = [("c1", "control", "D1", {A_TUMOR: (20, 20.0)}), ("c2", "control", "D1", {A_TUMOR: (20, 20.0)}),
            ("t1", "treated", "D1", {A_TUMOR: (30, 20.0)}), ("c3", "control", "D2", {A_TUMOR: (40, 40.0)}),
            ("t3", "treated", "D2", {A_TUMOR: (50, 40.0)})]
    p, table = project(tmp_path, data)
    row = subsets(compare(p, table)[0]).loc[A_TUMOR]
    assert row["per replicate"] == "D1 1.00 → 1.50; D2 1.00 → 1.25"


def test_reference_and_pools(tmp_path):
    c = {A_TUMOR: (40, 40.0)}
    data = [(f"{d}_{k}", k, d, c) for d in ("D1", "D2") for k in ("control", "treated")]
    p, table = project(tmp_path, data, experiment=lambda name: name.split("_")[1])
    p.reference = "treated"
    result, notes, _, reference = compare(p, table)
    assert reference == "treated" and set(result["condition"]) == {"control"}
    assert any("different pools" in n for n in notes) and "classified in different pools" in result.iloc[0]["note"]
    p.settings.classification.thresholds = {"A": 10.0}  # thresholds by hand: the same in every pool
    assert not any("pools" in n for n in compare(p, table)[1])
    p.settings.classification.thresholds = {}
    p.settings.classification.pool = Pool.DATASET
    assert any("on its own" in n for n in compare(p, table)[1])


def test_membrane_pairs(tmp_path):
    data = [(f"{d}_{k}", k, d, {}) for d in ("D1", "D2", "D3") for k in ("control", "treated")]
    p, table = project(tmp_path, data)
    shift = {"control": 0.0, "treated": 1.0}
    for k, (name, cond, _rep, _) in enumerate(data):
        rows = [{"membrane of": "Tumor", "interface with": "Immune A+", "cells with an interface": 30,
                 "median stdDevs above mean": 1.0 + shift[cond] + 0.1 * k, "marker check": "mainly on the measured "
                 "cells"},
                {"membrane of": "Tumor", "interface with": "Immune negative", "cells with an interface": 4,
                 "median stdDevs above mean": 0.5, "marker check": "partners carry part"},
                {"membrane of": "Immune A+", "interface with": "Tumor", "cells with an interface": 30,  # not this side
                 "median stdDevs above mean": 0.0, "marker check": ""},
                {"membrane of": "Tumor", "interface with": "Tumor", "cells with an interface": 30,  # a control pair
                 "median stdDevs above mean": 0.0, "marker check": ""}]
        pd.DataFrame(rows).assign(**{"marker": "MOI", "marker on": "Tumor"}).to_csv(
            tmp_path / "E" / name / "Membrane pairs.tsv", sep="\t", index=False)
    result = compare(p, table)[0]
    mem = result[result["measure"] == f"{MEMBRANE} (MOI)"].set_index(["cell type", "partner"])
    assert list(mem.index) == [("Tumor", "Immune A+"), ("Tumor", "Immune negative")]
    assert mem.loc[("Tumor", "Immune A+"), "effect"] == pytest.approx(1.1)
    low = mem.loc[("Tumor", "Immune negative")]
    assert np.isnan(low["p"]) and "too few interfaces" in low["note"] and "partners carry part" in low["note"]


def _pairs(marker, on, value, n=30):
    """Membrane pairs.tsv rows of one marker: its carrier with the other type, and the other way round."""
    other = "Immune" if on == "Tumor" else "Tumor"
    common = {"marker": marker, "marker on": on, "cells with an interface": n}
    return [common | {"membrane of": on, "interface with": other, "median stdDevs above mean": value,
                      "marker check": "mainly on the measured cells"},
            common | {"membrane of": other, "interface with": on, "median stdDevs above mean": -value,
                      "marker check": ""}]


def test_membrane_pairs_per_marker_of_interest(tmp_path):
    """Two markers of interest: a measure each, on the cell type that carries it (the column marker on); when the
    datasets disagree, the majority, with a note; an old file (one marker, no marker on) is not compared."""
    data = [(f"{d}_{k}", k, d, {}) for d in ("D1", "D2", "D3") for k in ("control", "treated")]
    p, table = project(tmp_path, data)
    for k, (name, cond, _rep, _) in enumerate(data):
        shift = 1.0 if cond == "treated" else 0.0
        on_b = "Tumor" if name == "D3_treated" else "Immune"  # one dataset says otherwise
        rows = _pairs("MOI", "Tumor", 1.0 + shift + 0.1 * k) + _pairs("B", on_b, 2.0 - shift + 0.1 * k)
        pd.DataFrame(rows).to_csv(tmp_path / "E" / name / "Membrane pairs.tsv", sep="\t", index=False)
    result, notes, _, _ = compare(p, table)
    first = result[result["measure"] == f"{MEMBRANE} (MOI)"].set_index(["cell type", "partner"])
    assert list(first.index) == [("Tumor", "Immune")] and first.iloc[0]["effect"] == pytest.approx(1.1)
    second = result[result["measure"] == f"{MEMBRANE} (B)"].set_index(["cell type", "partner"])
    assert list(second.index) == [("Immune", "Tumor")]  # the majority: on the immune cells
    assert "majority" in second.iloc[0]["note"] and any("compared on the immune cells" in n for n in notes)
    # files written before 2026-11-03: no marker columns, so no known carrier
    for name, *_ in data:
        pd.DataFrame(_pairs("x", "Tumor", 1.0)).drop(columns=["marker", "marker on"]).to_csv(
            tmp_path / "E" / name / "Membrane pairs.tsv", sep="\t", index=False)
    result, notes, _, _ = compare(p, table)
    assert not result["measure"].str.startswith(MEMBRANE).any() and any("run the analysis again" in n for n in notes)


def test_write_and_round_trip(tmp_path):
    c = {A_TUMOR: (40, 40.0)}
    data = [(f"{d}_{k}", k, d, c) for d in ("D1", "D2") for k in ("control", "treated")]
    p, table = project(tmp_path, data)
    p.reference = "control"
    p.save()
    q = Project.load(tmp_path)
    assert q.reference == "control"
    assert [(d.condition, d.replicate) for d in q.zips[0].datasets] == [(d.condition, d.replicate)
                                                                         for d in p.zips[0].datasets]
    write_comparison(q, table)
    assert (tmp_path / COMPARISON_FILE).exists()
    for d in q.zips[0].datasets:
        d.condition = ""
    write_comparison(q, table)  # nothing to compare any more: no stale file
    assert not (tmp_path / COMPARISON_FILE).exists()


def test_bh():
    assert list(bh([0.01, 0.04, 0.03])) == pytest.approx([0.03, 0.04, 0.04])
    assert list(bh([0.5])) == [0.5]
