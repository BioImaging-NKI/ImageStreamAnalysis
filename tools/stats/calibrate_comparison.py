"""Calibration of the comparison of conditions (comparison.py; manual: analysis/comparing-conditions.md).

Simulated experiments with a known truth, analysed with the tool's own tests (comparison._paired, _welch, bh):

1. the false-positive rate of the paired t-test across replicates, and of Welch's test, when the conditions do not
   differ (with differences between donors and between samples);
2. a test of one donor's two samples against each other (event sampling only), which the tool does not offer: why;
3. Benjamini-Hochberg over all subset pairs: the chance of any false finding;
4. the power for a real, small change of preference;
5. the membrane enrichment: medians of skewed per-cell scores.

The events are made as in the synthetic test data: one or two tumor cells and up to three immune cells of five
subtypes, placed touching at random, except that a cell of subtype A is placed next to a tumor cell with probability
p_A (the preference) and Gamma cells next to an immune cell with probability 0.8. Only the contacts are simulated,
not the images: this checks the statistics, not the segmentation.

    python tools/stats/calibrate_comparison.py            (about 3 minutes)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from imagestream_analysis.comparison import _paired, _welch, bh  # noqa: E402

COMPOSITIONS = [((1, 1), 0.44), ((2, 1), 0.12), ((1, 2), 0.15), ((1, 0), 0.06), ((0, 1), 0.06), ((2, 0), 0.04),
                ((0, 2), 0.05), ((2, 2), 0.04), ((1, 3), 0.04)]  # (tumor cells, immune cells) per event
KINDS = {"A": 0.28, "A+B": 0.06, "B": 0.22, "Gamma": 0.15, "neg": 0.26}
SUBTYPES = sorted(KINDS)
GRID = np.round(np.arange(0.20, 1.001, 0.05), 2)  # p_A of the event libraries
P0, TAU = 0.4, 0.1  # the preference of the control, and its SD between donors
ALPHA_TUMOR = (SUBTYPES.index("A"), 0)


def place(rng, kinds, p_a):
    cells = []
    for kind in kinds:
        r = rng.uniform(11, 14) if kind == "tumor" else rng.uniform(6, 8)
        for _ in range(300):
            if not cells:
                cy = cx = 0.0
            else:
                pool = cells
                if kind in ("A", "A+B") and rng.random() < p_a:
                    pool = [c for c in cells if c["kind"] == "tumor"] or cells
                elif kind == "Gamma" and rng.random() < 0.8:
                    pool = [c for c in cells if c["kind"] != "tumor"] or cells
                anchor = pool[rng.integers(len(pool))]
                a = rng.uniform(0, 2 * np.pi)
                d = anchor["r"] + r - 0.5
                cy, cx = anchor["cy"] + d * np.sin(a), anchor["cx"] + d * np.cos(a)
            if all(np.hypot(cy - c["cy"], cx - c["cx"]) >= c["r"] + r - 0.6 for c in cells):
                cells.append({"kind": kind, "cy": cy, "cx": cx, "r": r})
                break
    return cells


def library(p_a, n, seed):
    """Per event: observed and expected subset contacts (k subtypes x [tumor, k subtypes]), as subset_preference."""
    rng = np.random.default_rng(seed)
    kw = np.array(list(KINDS.values()))
    cw = np.array([w for _, w in COMPOSITIONS])
    k = len(SUBTYPES)
    obs, exp = np.zeros((n, k, k + 1)), np.zeros((n, k, k + 1))
    for e in range(n):
        nt, ni = COMPOSITIONS[rng.choice(len(COMPOSITIONS), p=cw / cw.sum())][0]
        cells = place(rng, ["tumor"] * nt + list(rng.choice(list(KINDS), ni, p=kw / kw.sum())), p_a)
        sub = np.array([-1 if c["kind"] == "tumor" else SUBTYPES.index(c["kind"]) for c in cells])
        m = int((sub >= 0).sum())
        if m < 2:
            continue  # no choice of immune cells
        c = np.bincount(sub[sub >= 0], minlength=k).astype(float)
        it = ii = 0
        for i in range(len(cells)):
            for j in range(i + 1, len(cells)):
                a, b = cells[i], cells[j]
                if np.hypot(a["cy"] - b["cy"], a["cx"] - b["cx"]) > a["r"] + b["r"] + 0.5:
                    continue
                si, sj = sub[i], sub[j]
                if si >= 0 and sj >= 0:
                    obs[e, si, sj + 1] += 1
                    obs[e, sj, si + 1] += 1
                    ii += 1
                elif max(si, sj) >= 0:
                    obs[e, max(si, sj), 0] += 1
                    it += 1
        exp[e, :, 0] = it * c / m
        exp[e, :, 1:] = ii * 2 * (np.outer(c, c) - np.diag(c)) / (m * (m - 1))
    return obs, exp


def snap(p):
    return float(GRID[np.argmin(abs(GRID - np.clip(p, GRID[0], GRID[-1])))])


class Sim:
    def __init__(self, n=20000):
        print(f"making event libraries ({len(GRID)} preferences x {n} events)...", flush=True)
        self.lib = {float(p): library(float(p), n, 1000 + i) for i, p in enumerate(GRID)}
        o, e = self.lib[P0]
        t, f = self.lib[0.95]
        print(f"  true A-tumor preference: {o[..., 0, 0].sum() / e[..., 0, 0].sum():.3f} (p_A {P0}), "
              f"{t[..., 0, 0].sum() / f[..., 0, 0].sum():.3f} (p_A 0.95)")

    def sample(self, rng, p, n_events):
        """Summed observed and expected contacts of one sample, and the events (for the per-sample error)."""
        o, e = self.lib[snap(p)]
        idx = rng.integers(len(o), size=n_events)
        return o[idx], e[idx]

    def donors(self, rng, n, delta=0.0, events=300, tech=0.0, unpaired=False):
        """Per donor the log ratios (all pairs) of control and treated; tech: sample-to-sample SD of the log ratio."""
        ref, cond, samples = [], [], []
        for _ in range(n):
            p = P0 + rng.normal(0, TAU)
            o0, e0 = self.sample(rng, p, events)
            o1, e1 = self.sample(rng, (P0 + rng.normal(0, TAU) if unpaired else p) + delta, events)
            y0 = np.log((o0.sum(0) + .5) / (e0.sum(0) + .5)) + rng.normal(0, tech, o0.shape[1:])
            y1 = np.log((o1.sum(0) + .5) / (e1.sum(0) + .5)) + rng.normal(0, tech, o1.shape[1:])
            ref.append(y0)
            cond.append(y1)
            samples.append(((o0, e0), (o1, e1), (y0, y1)))
        ok = np.minimum(np.min([e.sum(0) for (_, e), _, _ in samples], axis=0),
                        np.min([e.sum(0) for _, (_, e), _ in samples], axis=0)) >= 5
        return np.array(ref), np.array(cond), ok, samples


def p_of(test, ref, cond, pair):
    x = {str(i): float(v[pair]) for i, v in enumerate(ref)}
    y = {str(i) if test is _paired else f"b{i}": float(v[pair]) for i, v in enumerate(cond)}
    return test(x, y)[3]


def jackknife(o, e):
    """The event-sampling standard error of log((O + .5) / (E + .5)) (leave one event out)."""
    big_o, big_e = o.sum(), e.sum()
    loo = np.log((big_o - o + .5) / (big_e - e + .5))
    return np.sqrt((len(o) - 1) / len(o) * np.sum((loo - loo.mean()) ** 2))


def main():
    rng = np.random.default_rng(7)
    sim = Sim()
    n_sim = 2000

    print("\n1. No difference between the conditions: false-positive rate at p < 0.05 (A with tumor)")
    print("   donors | paired t | Welch (other donors) | paired t, sample-to-sample SD 0.05")
    for n in (2, 3, 4, 6):
        a = np.mean([p_of(_paired, *sim.donors(rng, n)[:2], ALPHA_TUMOR) < .05 for _ in range(n_sim)])
        b = np.mean([p_of(_welch, *sim.donors(rng, n, unpaired=True)[:2], ALPHA_TUMOR) < .05 for _ in range(n_sim)])
        c = np.mean([p_of(_paired, *sim.donors(rng, n, tech=0.05)[:2], ALPHA_TUMOR) < .05 for _ in range(n_sim)])
        print(f"   {n:6d} | {a:8.3f} | {b:20.3f} | {c:.3f}")

    print("\n2. One donor's two samples against each other (z = difference / event-sampling SE), no difference")
    print("   sample-to-sample SD of the log ratio | 300 events | 3000 events")
    for tech in (0.0, 0.03, 0.05, 0.1):
        rates = []
        for events in (300, 3000):
            hits = 0
            for _ in range(n_sim // 4):
                (o0, e0), (o1, e1), (y0, y1) = sim.donors(rng, 1, events=events, tech=tech)[3][0]
                i, j = ALPHA_TUMOR
                se = np.hypot(jackknife(o0[:, i, j], e0[:, i, j]), jackknife(o1[:, i, j], e1[:, i, j]))
                hits += abs(y1[i, j] - y0[i, j]) / se > 1.96
            rates.append(hits / (n_sim // 4))
        print(f"   {tech:36.2f} | {rates[0]:10.3f} | {rates[1]:.3f}")

    print("\n3. All subset pairs at once (no difference): chance of at least one finding")
    print("   donors, events | pairs tested | any p < 0.05 | any q < 0.05")
    for n, events in ((3, 300), (4, 1000), (6, 300)):
        any_p = any_q = tested = 0
        for _ in range(n_sim // 2):
            ref, cond, ok, _ = sim.donors(rng, n, events=events)
            ps = [p_of(_paired, ref, cond, pair) for pair in zip(*np.nonzero(ok), strict=True)
                  if pair[1] == 0 or pair[1] - 1 >= pair[0]]  # each immune pair once
            tested += len(ps)
            any_p += min(ps, default=1) < .05
            any_q += min(bh(ps), default=1) < .05
        print(f"   {n}, {events:5d} | {tested / (n_sim // 2):12.1f} | {any_p / (n_sim // 2):12.3f} | "
              f"{any_q / (n_sim // 2):.3f}")

    print(f"\n4. Power (p < 0.05): A placed next to a tumor cell with probability {P0} -> 0.95 after treatment")
    print("   events per sample | 2 donors | 3 | 4 | 6 | 10")
    for events in (300, 1000, 3000):
        row = [np.mean([p_of(_paired, *sim.donors(rng, n, delta=0.55, events=events)[:2], ALPHA_TUMOR) < .05
                        for _ in range(n_sim // 2)]) for n in (2, 3, 4, 6, 10)]
        print(f"   {events:17d} | " + " | ".join(f"{v:.2f}" for v in row))

    print("\n5. Membrane enrichment, no difference: paired t on per-donor medians of skewed per-cell scores")
    print("   cells with an interface | 2 donors | 3 | 4 | 6")

    def cells(n, shift):
        x = rng.normal(0.3, 1.0, n)
        tail = rng.random(n) < 0.25
        x[tail] += rng.exponential(3.0, tail.sum())
        return x + shift
    for lo, hi in ((10, 30), (30, 200)):
        row = []
        for n in (2, 3, 4, 6):
            hits = 0
            for _ in range(n_sim):
                u = rng.normal(0, 0.5, n)
                x = {str(i): float(np.median(cells(rng.integers(lo, hi), u[i]))) for i in range(n)}
                y = {str(i): float(np.median(cells(rng.integers(lo, hi), u[i]))) for i in range(n)}
                hits += _paired(x, y)[3] < .05
            row.append(hits / n_sim)
        print(f"   {lo:3d}-{hi:<19d} | " + " | ".join(f"{v:.3f}" for v in row))


if __name__ == "__main__":
    main()
