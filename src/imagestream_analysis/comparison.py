"""Comparing conditions across replicates (manual: analysis/comparing-conditions.md).

Every analysed dataset with a condition (control, treated, ...) and a replicate (the independent unit: a donor, an
independent experiment) takes part. The datasets of one replicate and condition are pooled (their contacts summed,
their membranes weighted by cells). Per condition against the reference condition, for every subset preference and
membrane pair:

- paired, when two or more replicates have both conditions: a paired t-test on the per-replicate differences
  (log observed/expected for the subset preference, median SD above mean for the membrane enrichment); replicates
  with only one of the two conditions are left out;
- otherwise unpaired, when each condition has two or more replicates: Welch's t-test;
- otherwise no test (the per-replicate values only).

The variation between replicates is the error term: a test of one replicate's two samples against each other (the
event-sampling interval of each dataset) would also count every difference between two samples as an effect, and
was not calibrated in simulations (calibration in the manual). The p values of one measure and condition are adjusted
for the number of pairs tested (Benjamini-Hochberg: q)."""

from __future__ import annotations

import logging
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CellOfInterest, Pool
from .interactions import MIN_EXPECTED
from .io import save_table

if TYPE_CHECKING:
    from .project import Project

log = logging.getLogger(__name__)
COMPARISON_FILE = "condition comparison.tsv"
MIN_INTERFACES = 10  # fewer cells with an interface in a replicate: the membrane pair is not tested
SUBSETS = "subset preference"
MEMBRANE = "membrane enrichment"
OK_CHECK = "mainly on the measured cells"  # the marker check that needs no note (membrane.CHECK_LEVELS)
COLUMNS = ["measure", "cell type", "partner", "condition", "reference", "design", "replicates", "effect", "95% low",
           "95% high", "p", "q", "per replicate", "note"]


def labelled(project: Project, table: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """The analysed datasets with a condition and a replicate (folder, experiment, condition, replicate), and notes
    on the labelled datasets that cannot take part."""
    errors = table["error"].fillna("").astype(str) if "error" in table else pd.Series("", index=table.index)
    done = {Path(str(f)).as_posix() for f, e in zip(table["folder"], errors, strict=True) if not e}
    rows, notes, unlabelled = [], [], []
    for z, d, s in project.samples():
        if not d.condition:
            continue
        if not d.replicate:
            unlabelled.append(d.name)
            continue
        folder = (Path(s.input.output) / d.name).relative_to(project.folder).as_posix()
        if folder not in done:
            notes.append(f"{d.name}: no results in {folder} (run the analysis first)")
            continue
        rows.append({"folder": folder, "dataset": d.name, "experiment": d.experiment or Path(z.path).stem,
                     "condition": d.condition, "replicate": d.replicate})
    if unlabelled:
        notes.append("Without a replicate, left out: " + ", ".join(unlabelled))
    return pd.DataFrame(rows, columns=["folder", "dataset", "experiment", "condition", "replicate"]), notes


def reference_of(project: Project, conditions: list[str]) -> str:
    """The reference condition: as set in the project, else the first condition in the project's order."""
    return project.reference if project.reference in conditions else (conditions[0] if conditions else "")


def compare(project: Project, table: pd.DataFrame) -> tuple[pd.DataFrame, list[str], list[str], str]:
    """The comparison table (COLUMNS), notes, the conditions (project order) and the reference."""
    units, notes = labelled(project, table)
    conditions = list(dict.fromkeys(units["condition"]))
    reference = reference_of(project, conditions)
    if len(conditions) < 2:
        if conditions:
            notes.append(f"Only one condition ({conditions[0]}): nothing to compare")
        return pd.DataFrame(columns=COLUMNS), notes, conditions, reference
    pool_note = _pool_note(project, units)
    if pool_note:
        notes.append(pool_note)
    out = []
    for condition in conditions:
        if condition == reference:
            continue
        a, b = units[units["condition"] == reference], units[units["condition"] == condition]
        design, used, note = _design(sorted(set(a["replicate"])), sorted(set(b["replicate"])))
        notes.append(f"{condition} vs {reference}: {note}")
        for measure, values in ((SUBSETS, _subsets), (MEMBRANE, _membranes)):
            rows = _compare_measure(project, measure, values, a, b, design, used)
            for r in rows:
                r |= {"condition": condition, "reference": reference}
                if pool_note:
                    r["note"] = "; ".join(x for x in (r["note"], "classified in different pools") if x)
            out += sorted(_adjust(rows), key=lambda r: (not np.isfinite(r["p"]), r["p"] if np.isfinite(r["p"]) else 0))
    return pd.DataFrame(out, columns=COLUMNS), notes, conditions, reference


def _design(ref: list[str], cond: list[str]) -> tuple[str, dict[str, list[str]], str]:
    """'paired' (replicates with both conditions), 'unpaired' or '' (no test), the replicates used per side, and a
    note on the design."""
    both = sorted(set(ref) & set(cond))
    if len(both) >= 2:
        left = sorted(set(ref) ^ set(cond))
        extra = f"; without both conditions, left out: {', '.join(left)}" if left else ""
        few = "; with two replicates the interval is wide (1 degree of freedom)" if len(both) == 2 else ""
        return "paired", {"ref": both, "cond": both}, (f"paired over {len(both)} replicates ({', '.join(both)})"
                                                       f"{extra}{few}")
    if len(ref) >= 2 and len(cond) >= 2:
        return "unpaired", {"ref": ref, "cond": cond}, (
            f"unpaired: {len(ref)} and {len(cond)} replicates (Welch); give both conditions of a donor the same "
            "replicate to pair them")
    return "", {"ref": ref, "cond": cond}, "no test: two or more replicates per condition are needed"


ValueFn = Callable[[Path, pd.DataFrame, "Project"], dict[tuple[str, str], dict[str, Any]]]


def _compare_measure(project: Project, measure: str, values: ValueFn, a: pd.DataFrame, b: pd.DataFrame, design: str,
                     used: dict[str, list[str]]) -> list[dict[str, Any]]:
    """One row per item (cell type, partner) of the measure."""
    folder = Path(project.folder)
    per_unit = {("ref", r): values(folder, a[a["replicate"] == r], project) for r in used["ref"]}
    per_unit |= {("cond", r): values(folder, b[b["replicate"] == r], project) for r in used["cond"]}
    items = sorted({k for v in per_unit.values() for k in v})
    rows = []
    for item in items:
        got = {u: v.get(item) for u, v in per_unit.items()}
        short = [f"{r} ({'reference' if side == 'ref' else 'condition'})" for (side, r), g in got.items()
                 if g is None or not g["ok"]]
        ref = {r: got[("ref", r)] for r in used["ref"]}
        cond = {r: got[("cond", r)] for r in used["cond"]}
        row = {"measure": measure, "cell type": item[0], "partner": item[1], "design": design,
               "replicates": len(used["ref"]) if design == "paired" else f"{len(used['ref'])} + {len(used['cond'])}",
               "per replicate": _per_replicate(design, ref, cond), "note": "", "effect": np.nan, "95% low": np.nan,
               "95% high": np.nan, "p": np.nan, "q": np.nan}
        checks = sorted({c for g in got.values() if g for c in g.get("checks", ())} - {"", "none", OK_CHECK})
        notes = [("too few contacts" if measure == SUBSETS else "too few interfaces") + " in " + ", ".join(short)
                 ] if short else []
        if checks:
            notes.append("marker check: " + ", ".join(checks))
        row["note"] = "; ".join(notes)
        if design and not short:
            x = {r: g["value"] for r, g in ref.items() if g}
            y = {r: g["value"] for r, g in cond.items() if g}
            est, low, high, p = (_paired(x, y) if design == "paired" else _welch(x, y))
            if measure == SUBSETS:  # the ratio of the two observed / expected
                est, low, high = np.exp(est), np.exp(low), np.exp(high)
            row |= {"effect": est, "95% low": low, "95% high": high, "p": p}
        rows.append(row)
    return rows


def _subsets(folder: Path, units: pd.DataFrame, project: Project) -> dict[tuple[str, str], dict[str, Any]]:
    """Per pair: the contacts and expected contacts summed over the datasets of one replicate and condition; value
    log((O + 0.5) / (E + 0.5)); ok when E >= MIN_EXPECTED. Each immune–immune pair once (S–T and T–S are the same
    contacts)."""
    parts = [pd.read_csv(f, sep="\t") for f in (folder / u / "Subset preference.tsv" for u in units["folder"])
             if f.exists()]
    if not parts:
        return {}
    t = pd.concat(parts, ignore_index=True)
    t = t[(t["neighbour type"] == "Tumor") | (t["cell type"] <= t["neighbour type"])]
    s = t.groupby(["cell type", "neighbour type"])[["contacts", "expected (subtypes shuffled)"]].sum()
    out = {}
    for key, o, e in zip(s.index, s["contacts"].astype(float), s["expected (subtypes shuffled)"].astype(float),
                         strict=True):
        out[(str(key[0]), str(key[1]))] = {"value": float(np.log((o + 0.5) / (e + 0.5))),
                                           "show": o / e if e > 0 else np.nan, "ok": e >= MIN_EXPECTED}
    return out


def _membranes(folder: Path, units: pd.DataFrame, project: Project) -> dict[tuple[str, str], dict[str, Any]]:
    """Per membrane pair of the cell of interest with the other cell type: the median SD above mean of the datasets
    of one replicate and condition (weighted by their cells with an interface); ok with MIN_INTERFACES of them."""
    parts = [pd.read_csv(f, sep="\t") for f in (folder / u / "Membrane pairs.tsv" for u in units["folder"])
             if f.exists()]
    if not parts:
        return {}
    t = pd.concat(parts, ignore_index=True)
    tumor = project.settings.channels.cell_of_interest == CellOfInterest.TUMOR_CELL
    mine = t["membrane of"].astype(str).str.startswith("Tumor") == tumor
    theirs = t["interface with"].astype(str).str.startswith("Tumor") != tumor
    t = t[mine & theirs & t["median stdDevs above mean"].notna()]
    out = {}
    for (m, w), g in t.groupby(["membrane of", "interface with"]):
        n = g["cells with an interface"].to_numpy(dtype=float)
        v = g["median stdDevs above mean"].to_numpy(dtype=float)
        value = float(np.sum(n * v) / n.sum()) if n.sum() else np.nan
        out[(str(m), str(w))] = {"value": value, "show": value, "ok": n.sum() >= MIN_INTERFACES,
                                 "checks": set(g["marker check"].fillna("").astype(str))}
    return out


def _paired(x: dict[str, float], y: dict[str, float]) -> tuple[float, float, float, float]:
    d = np.array([y[r] - x[r] for r in sorted(x)])
    n = len(d)
    mean, se = float(d.mean()), float(d.std(ddof=1) / np.sqrt(n))
    if not se > 0:
        return mean, np.nan, np.nan, np.nan
    half = stats.t.ppf(0.975, n - 1) * se
    return mean, mean - half, mean + half, float(2 * stats.t.sf(abs(mean / se), n - 1))


def _welch(x: dict[str, float], y: dict[str, float]) -> tuple[float, float, float, float]:
    a, b = np.array(list(x.values())), np.array(list(y.values()))
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    mean, se = float(b.mean() - a.mean()), float(np.sqrt(va + vb))
    if not se > 0:
        return mean, np.nan, np.nan, np.nan
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    half = stats.t.ppf(0.975, df) * se
    return mean, mean - half, mean + half, float(2 * stats.t.sf(abs(mean / se), df))


def _per_replicate(design: str, ref: dict[str, dict[str, Any] | None], cond: dict[str, dict[str, Any] | None]) -> str:
    def f(g: dict[str, Any] | None) -> str:
        return "–" if g is None or not np.isfinite(g["show"]) else f"{g['show']:.2f}"
    if design == "paired":
        return "; ".join(f"{r} {f(ref[r])} → {f(cond[r])}" for r in sorted(ref))
    return ("reference " + ", ".join(f"{r} {f(g)}" for r, g in sorted(ref.items())) + "; condition "
            + ", ".join(f"{r} {f(g)}" for r, g in sorted(cond.items())))


def _adjust(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Benjamini-Hochberg q values over the rows with a p value."""
    tested = [r for r in rows if np.isfinite(r["p"])]
    q = bh([r["p"] for r in tested])
    for r, v in zip(tested, q, strict=True):
        r["q"] = float(v)
    return rows


def bh(p: list[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values (q)."""
    a = np.asarray(p, dtype=float)
    n = len(a)
    if not n:
        return a
    order = np.argsort(a)
    q = np.minimum.accumulate((a[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1)
    return out


def _pool_note(project: Project, units: pd.DataFrame) -> str:
    """A note when the compared datasets were not classified together and the subtype thresholds are automatic."""
    c = project.settings.classification
    markers = project.settings.channels.functional or []
    if markers and all(m in c.thresholds for m in markers):
        return ""
    if c.pool == Pool.DATASET and len(units) > 1:
        return ("Every dataset was classified on its own: a subtype can differ between conditions by its threshold "
                "alone; set the thresholds by hand or classify per experiment")
    experiments = sorted(set(units["experiment"]))
    if len(experiments) > 1:
        return (f"The datasets were classified in different pools ({', '.join(experiments)}): a subtype can differ "
                "between conditions by its threshold alone; give them one experiment label or set the thresholds by "
                "hand")
    return ""


def write_comparison(project: Project, table: pd.DataFrame) -> pd.DataFrame:
    """condition comparison.tsv in the project folder (removed when no conditions can be compared)."""
    result, notes, _, _ = compare(project, table)
    path = Path(project.folder) / COMPARISON_FILE
    for n in notes:
        log.info("Comparison: %s", n)
    if len(result):
        save_table(result, path)
    elif path.exists():
        path.unlink()
    return result
