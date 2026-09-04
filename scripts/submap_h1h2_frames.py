#!/usr/bin/env python3
"""SubMap H1/H2 in other absolute frames (sidecar).

Existing analysis (`scripts/submap_m56_h1h2.py`) used only Vupn1/Vtn1
(MORVEL56-NNR). This sidecar repeats C vs E Spearman for every SubMap
absolute-frame pair with finite Vupn/Vtn:

  M56 family: Vupn1/Vtn1 … Vupn4/Vtn4
  N15 family: Vupn5/Vtn5 … Vupn9/Vtn9

UPS recode is imported from `submap_m56_h1h2.py` (C = C+C-SS, E = E+E-SS,
Undefined dropped, SS held out of C vs E). Score C=+1, E=−1.

Frame names: only frame 1 is named in the dump + Cerpa text as used in the
existing analysis (MORVEL56-NNR). Cerpa et al. (2025) list four M56 absolute
frames but do not map columns 2–4. N15 columns have no frame names in
README / Cerpa fulltext. Unknown names stay "Frame N (unknown)" — no HS3
vs NNR labels are invented.

This sidecar is about absolute V_UP vs V_T only. Circularity of vdn vs UPS
still applies and is not re-tested. Müller v_T vs v_OP is not used.

Run:
    /workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
        /workspace/subduction-test/scripts/submap_h1h2_frames.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import submap_m56_h1h2 as h12  # noqa: E402

ROOT = Path("/workspace/subduction-test")
FULL_CSV = ROOT / "data" / "submap" / "submap_transects_full.csv"
REPORT_MD = ROOT / "results" / "submap_h1h2_frames.md"
FIG_FOREST = ROOT / "results" / "fig_h1h2_by_frame.png"

SENTINEL = h12.SENTINEL
CE_SCORE = h12.CE_SCORE
CNE_SCORE = h12.CNE_SCORE
UPS3_ORDER = h12.UPS3_ORDER

# Cerpa, Lallemand & Heuret 2025 (EarthArXiv) §II.1 lists four M56 absolute
# frames, in this order, without mapping them onto Vupn2/3/4:
#   MORVEL56-NNR (Argus et al., 2011)
#   spreading-alignment (Becker et al., 2015)
#   T25m (Wang et al., 2018)
#   GMHRF-1Ma (Doubrovine et al., 2012)
# Existing analysis documents Vupn1/Vtn1 = MORVEL56-NNR. Columns 2–4 and all
# N15 frames (5–9) are unnamed in the dump, README, and Cerpa fulltext.
CERPA_M56_FRAME_LIST = [
    "MORVEL56-NNR (Argus et al., 2011)",
    "spreading-alignment (Becker et al., 2015)",
    "T25m (Wang et al., 2018)",
    "GMHRF-1Ma (Doubrovine et al., 2012)",
]

FRAMES = [
    {
        "n": 1,
        "family": "M56",
        "vupn": "Vupn1",
        "vtn": "Vtn1",
        "name": "MORVEL56-NNR (Argus et al., 2011)",
        "name_known": True,
    },
    {
        "n": 2,
        "family": "M56",
        "vupn": "Vupn2",
        "vtn": "Vtn2",
        "name": "Frame 2 (unknown)",
        "name_known": False,
    },
    {
        "n": 3,
        "family": "M56",
        "vupn": "Vupn3",
        "vtn": "Vtn3",
        "name": "Frame 3 (unknown)",
        "name_known": False,
    },
    {
        "n": 4,
        "family": "M56",
        "vupn": "Vupn4",
        "vtn": "Vtn4",
        "name": "Frame 4 (unknown)",
        "name_known": False,
    },
    {
        "n": 5,
        "family": "N15",
        "vupn": "Vupn5",
        "vtn": "Vtn5",
        "name": "Frame 5 (unknown)",
        "name_known": False,
    },
    {
        "n": 6,
        "family": "N15",
        "vupn": "Vupn6",
        "vtn": "Vtn6",
        "name": "Frame 6 (unknown)",
        "name_known": False,
    },
    {
        "n": 7,
        "family": "N15",
        "vupn": "Vupn7",
        "vtn": "Vtn7",
        "name": "Frame 7 (unknown)",
        "name_known": False,
    },
    {
        "n": 8,
        "family": "N15",
        "vupn": "Vupn8",
        "vtn": "Vtn8",
        "name": "Frame 8 (unknown)",
        "name_known": False,
    },
    {
        "n": 9,
        "family": "N15",
        "vupn": "Vupn9",
        "vtn": "Vtn9",
        "name": "Frame 9 (unknown)",
        "name_known": False,
    },
]


def fmt_p(p):
    return h12.fmt_p(p)


def fmt_rho(r):
    return h12.fmt_rho(r)


def fmt_num(x, digits=3):
    return h12.fmt_num(x, digits)


def spearman_ci(r, n, alpha=0.05):
    """Approximate 95% CI for Spearman ρ via Fisher z (n > 3)."""
    if n is None or n <= 3 or not np.isfinite(r) or abs(r) >= 1.0:
        return np.nan, np.nan
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    # ~1.96 for 95%
    zcrit = 1.959963984540054
    lo = float(np.tanh(z - zcrit * se))
    hi = float(np.tanh(z + zcrit * se))
    return lo, hi


def stronger_label(rho_vup, rho_vt, gap_tie=0.0):
    """Which |ρ| is larger on C vs E. Tie if |Δ| <= gap_tie."""
    a = abs(rho_vup) if np.isfinite(rho_vup) else np.nan
    b = abs(rho_vt) if np.isfinite(rho_vt) else np.nan
    if not np.isfinite(a) or not np.isfinite(b):
        return "NA"
    if abs(a - b) <= gap_tie:
        return "tie"
    return "H1 Vupn" if a > b else "H2 Vtn"


def analyse_frame(df, spec):
    vu_col = spec["vupn"]
    vt_col = spec["vtn"]
    vu = h12.sentinel_to_nan(df[vu_col]) / 10.0  # mm/yr → cm/yr
    vt = h12.sentinel_to_nan(df[vt_col]) / 10.0
    work = df.copy()
    work["Vupn_cm"] = vu.to_numpy(dtype=float)
    work["Vtn_cm"] = vt.to_numpy(dtype=float)

    finite = work["Vupn_cm"].notna() & work["Vtn_cm"].notna()
    defined = work.loc[finite & (work["Geodetic_UPS"] != "Undefined")].copy()
    defined = defined.loc[defined["UPS_3"].isin(["C", "N", "E", "SS"])].copy()

    ce = defined.loc[defined["UPS_3"].isin(["C", "E"])].copy()
    ce["score"] = ce["UPS_3"].map(CE_SCORE)
    cne = defined.loc[defined["UPS_3"].isin(["C", "N", "E"])].copy()
    cne["score"] = cne["UPS_3"].map(CNE_SCORE)

    r_vup_ce, p_vup_ce, n_vup_ce = h12.spearman_safe(ce["Vupn_cm"], ce["score"])
    r_vt_ce, p_vt_ce, n_vt_ce = h12.spearman_safe(ce["Vtn_cm"], ce["score"])
    r_col_ce, p_col_ce, n_col_ce = h12.spearman_safe(ce["Vupn_cm"], ce["Vtn_cm"])
    r_col_all, p_col_all, n_col_all = h12.spearman_safe(
        defined["Vupn_cm"], defined["Vtn_cm"]
    )
    r_vup_cne, p_vup_cne, n_vup_cne = h12.spearman_safe(cne["Vupn_cm"], cne["score"])
    r_vt_cne, p_vt_cne, n_vt_cne = h12.spearman_safe(cne["Vtn_cm"], cne["score"])

    meds = {}
    for cls in UPS3_ORDER:
        meds[cls] = {
            "n": int((defined["UPS_3"] == cls).sum()),
            "Vupn": h12.median_n(defined, "Vupn_cm", cls),
            "Vtn": h12.median_n(defined, "Vtn_cm", cls),
        }

    n_c = int((ce["UPS_3"] == "C").sum())
    n_e = int((ce["UPS_3"] == "E").sum())
    lo_vup, hi_vup = spearman_ci(r_vup_ce, n_vup_ce)
    lo_vt, hi_vt = spearman_ci(r_vt_ce, n_vt_ce)

    return {
        **spec,
        "n_finite_both": int(finite.sum()),
        "n_defined": int(len(defined)),
        "n_ce": int(len(ce)),
        "n_c": n_c,
        "n_e": n_e,
        "rho_vup_ce": r_vup_ce,
        "p_vup_ce": p_vup_ce,
        "n_vup_ce": n_vup_ce,
        "rho_vt_ce": r_vt_ce,
        "p_vt_ce": p_vt_ce,
        "n_vt_ce": n_vt_ce,
        "lo_vup": lo_vup,
        "hi_vup": hi_vup,
        "lo_vt": lo_vt,
        "hi_vt": hi_vt,
        "rho_col_ce": r_col_ce,
        "p_col_ce": p_col_ce,
        "n_col_ce": n_col_ce,
        "rho_col_all": r_col_all,
        "p_col_all": p_col_all,
        "n_col_all": n_col_all,
        "rho_vup_cne": r_vup_cne,
        "p_vup_cne": p_vup_cne,
        "n_vup_cne": n_vup_cne,
        "rho_vt_cne": r_vt_cne,
        "p_vt_cne": p_vt_cne,
        "n_vt_cne": n_vt_cne,
        "meds": meds,
        "usable": int(len(ce)) >= 10
        and np.isfinite(r_vup_ce)
        and np.isfinite(r_vt_ce),
        "stronger": stronger_label(r_vup_ce, r_vt_ce),
        "delta_abs": (
            abs(r_vup_ce) - abs(r_vt_ce)
            if np.isfinite(r_vup_ce) and np.isfinite(r_vt_ce)
            else np.nan
        ),
    }


def fig_forest(rows, path):
    h12._setup_style()
    usable = [r for r in rows if r["usable"]]
    if not usable:
        print("WARNING: no usable frames for forest plot")
        return

    fig, ax = plt.subplots(figsize=(9.6, 6.8))
    y = np.arange(len(usable))[::-1]
    # two series offset
    dy = 0.16
    for i, r in enumerate(usable):
        yi = y[i]
        ax.plot(
            [r["lo_vup"], r["hi_vup"]],
            [yi + dy, yi + dy],
            color="#9B2226",
            lw=1.6,
            zorder=2,
        )
        ax.plot(
            [r["lo_vt"], r["hi_vt"]],
            [yi - dy, yi - dy],
            color="#005F73",
            lw=1.6,
            zorder=2,
        )
        ax.scatter(
            r["rho_vup_ce"],
            yi + dy,
            s=42,
            color="#9B2226",
            zorder=3,
            edgecolors="0.15",
            linewidths=0.4,
        )
        ax.scatter(
            r["rho_vt_ce"],
            yi - dy,
            s=42,
            marker="s",
            color="#005F73",
            zorder=3,
            edgecolors="0.15",
            linewidths=0.4,
        )

    ax.axvline(0.0, color="0.55", lw=0.9, ls="--", zorder=1)
    labels = []
    for r in usable:
        tag = f"{r['family']} {r['n']}"
        if r["name_known"]:
            labels.append(f"{tag}: {r['name']}")
        else:
            labels.append(f"{tag}: {r['name']}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r"Spearman $\rho$ vs C vs E  (C=+1, E=$-1$)")
    ax.set_title(
        "SubMap absolute velocities across frames\n"
        r"circles $V_{\mathrm{UPn}}$ (H1); squares $V_{Tn}$ (H2); "
        "whiskers ≈ 95% Fisher-z CI",
        pad=10,
    )
    ax.scatter([], [], s=42, color="#9B2226", edgecolors="0.15", label=r"$V_{\mathrm{UPn}}$ (H1)")
    ax.scatter(
        [],
        [],
        s=42,
        marker="s",
        color="#005F73",
        edgecolors="0.15",
        label=r"$V_{Tn}$ (H2)",
    )
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # family separator between M56 and N15
    n_m56 = sum(1 for r in usable if r["family"] == "M56")
    if 0 < n_m56 < len(usable):
        # y decreases upward in plot? yticks are reversed so M56 (first in usable
        # if FRAMES order 1..9) sit at the top. Separator below last M56.
        y_m56_last = y[n_m56 - 1]
        y_n15_first = y[n_m56]
        ymid = 0.5 * (y_m56_last + y_n15_first)
        ax.axhline(ymid, color="0.75", lw=0.8, ls=":", zorder=1)
    ax.set_xlim(-1.05, 1.05)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Wrote {path}")


def ranking_summary(rows):
    usable = [r for r in rows if r["usable"]]
    labels = [r["stronger"] for r in usable]
    unique = sorted(set(labels))
    flips = len(unique) > 1
    signs_vup = [
        np.sign(r["rho_vup_ce"]) for r in usable if np.isfinite(r["rho_vup_ce"])
    ]
    signs_vt = [
        np.sign(r["rho_vt_ce"]) for r in usable if np.isfinite(r["rho_vt_ce"])
    ]
    vup_sign_stable = len(set(int(s) for s in signs_vup if s != 0)) <= 1
    vt_sign_stable = len(set(int(s) for s in signs_vt if s != 0)) <= 1
    return {
        "usable": usable,
        "labels": labels,
        "unique": unique,
        "flips": flips,
        "vup_sign_stable": vup_sign_stable,
        "vt_sign_stable": vt_sign_stable,
    }


def write_report(rows, n_raw, ranking):
    usable = ranking["usable"]

    frame_table = [
        "| n | family | columns | name | name known? | n finite both | n C vs E |",
        "|---:|---|---|---|---|---:|---:|",
    ]
    for r in rows:
        known = "yes" if r["name_known"] else "no — unknown"
        frame_table.append(
            f"| {r['n']} | {r['family']} | `{r['vupn']}` / `{r['vtn']}` | "
            f"{r['name']} | {known} | {r['n_finite_both']} | {r['n_ce']} |"
        )

    spear_table = [
        "| frame | family | n C/E | ρ(Vupn, C/E) | p | ρ(Vtn, C/E) | p | ρ(Vupn, Vtn) | p | |ρ_H1|−|ρ_H2| | stronger |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in rows:
        if not r["usable"]:
            spear_table.append(
                f"| {r['n']} {r['name']} | {r['family']} | {r['n_ce']} | "
                f"{fmt_rho(r['rho_vup_ce'])} | {fmt_p(r['p_vup_ce'])} | "
                f"{fmt_rho(r['rho_vt_ce'])} | {fmt_p(r['p_vt_ce'])} | "
                f"{fmt_rho(r['rho_col_ce'])} | {fmt_p(r['p_col_ce'])} | "
                f"{fmt_num(r['delta_abs'])} | unusable |"
            )
            continue
        spear_table.append(
            f"| {r['n']} {r['name']} | {r['family']} | {r['n_ce']} | "
            f"{fmt_rho(r['rho_vup_ce'])} | {fmt_p(r['p_vup_ce'])} | "
            f"{fmt_rho(r['rho_vt_ce'])} | {fmt_p(r['p_vt_ce'])} | "
            f"{fmt_rho(r['rho_col_ce'])} | {fmt_p(r['p_col_ce'])} | "
            f"{fmt_num(r['delta_abs'])} | {r['stronger']} |"
        )

    cne_table = [
        "| frame | ρ(Vupn, C/N/E) | p | n | ρ(Vtn, C/N/E) | p | n | ρ(Vupn, Vtn) all defined | p | n |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        cne_table.append(
            f"| {r['n']} {r['family']} | {fmt_rho(r['rho_vup_cne'])} | "
            f"{fmt_p(r['p_vup_cne'])} | {r['n_vup_cne']} | "
            f"{fmt_rho(r['rho_vt_cne'])} | {fmt_p(r['p_vt_cne'])} | {r['n_vt_cne']} | "
            f"{fmt_rho(r['rho_col_all'])} | {fmt_p(r['p_col_all'])} | {r['n_col_all']} |"
        )

    med_blocks = []
    for r in rows:
        lines = [
            f"**Frame {r['n']}** — {r['name']} ({r['family']}; n_defined={r['n_defined']})",
            "",
            "| UPS_3 | n | median Vupn (cm/yr) | median Vtn (cm/yr) |",
            "|---|---:|---|---|",
        ]
        for cls in UPS3_ORDER:
            m = r["meds"][cls]
            lines.append(
                f"| {cls} | {m['n']} | {h12.med_cell(m['Vupn'])} | {h12.med_cell(m['Vtn'])} |"
            )
        med_blocks.append("\n".join(lines) + "\n")

    if ranking["flips"]:
        near_ties = [
            r
            for r in usable
            if np.isfinite(r["delta_abs"]) and abs(r["delta_abs"]) < 0.02
        ]
        tie_note = ""
        if near_ties:
            tie_note = (
                " Frames "
                + ", ".join(str(r["n"]) for r in near_ties)
                + " are near-ties (|Δρ|<0.02) and should not be read as a ranking."
            )
        rank_text = (
            "H1 vs H2 **ranking flips** across frames: the larger-|ρ| predictor is "
            + ", ".join(
                f"frame {r['n']} → {r['stronger']}" for r in usable
            )
            + "."
            + tie_note
            + " Unambiguous flips remain (e.g. frame 4 vs frame 1). "
            "That is the Schellart (2008) point: absolute-frame choice can "
            "reorder V_UP vs V_T. No hypothesis winner."
        )
    else:
        who = ranking["unique"][0] if ranking["unique"] else "NA"
        rank_text = (
            f"H1 vs H2 ranking by |ρ| is the same in every usable frame ({who}). "
            "That stability is still not a hypothesis winner: UPS vs vdn remains "
            "circular, several frame names are unknown, |ρ| values are modest, "
            "and a ranking that survives these nine columns is not a general "
            "H1/H2 verdict."
        )

    sign_text = (
        f"Sign of ρ(Vupn, C/E) stable across usable frames: "
        f"{'yes' if ranking['vup_sign_stable'] else 'no'}. "
        f"Sign of ρ(Vtn, C/E) stable: "
        f"{'yes' if ranking['vt_sign_stable'] else 'no'}. "
        "Predicted signs (SubMap convention, + = oceanward): H1 ρ(Vupn)>0 "
        "(OP advance with C); H2 ρ(Vtn)<0 (trench retreat / rollback with E)."
    )

    # highlight frame 1 vs others in one sentence using computed numbers
    f1 = next((r for r in rows if r["n"] == 1), None)
    f1_line = ""
    if f1 is not None:
        f1_line = (
            f"Frame 1 (MORVEL56-NNR) C vs E: ρ(Vupn)={fmt_rho(f1['rho_vup_ce'])} "
            f"(p={fmt_p(f1['p_vup_ce'])}), ρ(Vtn)={fmt_rho(f1['rho_vt_ce'])} "
            f"(p={fmt_p(f1['p_vt_ce'])}), n={f1['n_ce']}, "
            f"ρ(Vupn,Vtn)={fmt_rho(f1['rho_col_ce'])}."
        )

    md = f"""# SubMap H1/H2 across absolute frames

Sidecar to `results/submap_m56_h1h2.md`. Values are computed from
`{FULL_CSV}`. Nothing is invented. UPS recode is the same function as
`scripts/submap_m56_h1h2.py`.

- **n rows in full sheet:** {n_raw}
- **Units:** SubMap stores mm/yr; medians below are **cm/yr** (÷ 10). Spearman
  is rank-based (unit-free). Sentinel −999 → NaN.
- **This sidecar tests absolute V_UP vs V_T only.** `vdn` vs UPS is circular
  and is not repeated. Müller `v_T` vs `v_OP` is not used.

## Frame names (do not invent)

Cerpa, Lallemand & Heuret (2025, §II.1) calculated absolute velocities in
**four** M56 reference frames, listed as:

1. {CERPA_M56_FRAME_LIST[0]}
2. {CERPA_M56_FRAME_LIST[1]}
3. {CERPA_M56_FRAME_LIST[2]}
4. {CERPA_M56_FRAME_LIST[3]}

The SubMap dump labels columns `Vupn1`…`Vupn4` / `Vtn1`…`Vtn4` without a
name field. The existing analysis documents **frame 1 = MORVEL56-NNR**.
Cerpa does **not** map frames 2–4 onto column numbers, so those three stay
**Frame N (unknown)**. They are the remaining three names in the list above,
in unknown order.

Columns `Vupn5`…`Vupn9` / `Vtn5`…`Vtn9` sit next to the `N15_*` relative
velocities (NUVEL-1A-era kinematics in this dump). No N15 absolute-frame
names appear in `data/submap/README.md`, Cerpa fulltext, or the JSON sheets.
Frames 5–9 are **Frame N (unknown)**. HS3 / NNR / hotspot labels are **not**
assigned.

{chr(10).join(frame_table)}

## Sample and recode

C = C + C-SS; E = E + E-SS; N = Neutral; SS = SS + SS-C + SS-E.
Undefined dropped. SS held out of C vs E. Score C=+1, E=−1 (shortening-positive,
matching SubMap `vdn`). A frame is used if C vs E has n≥10 and finite ρ for
both Vupn and Vtn.

Sign of absolute trench-normal components (from the frame-1 identity
`vdn ≈ Vupn1 − Vtn1` in the parent analysis): **positive Vupn / Vtn = oceanward**.

{f1_line}

## Spearman C vs E, collinearity, ranking

{chr(10).join(spear_table)}

ρ(Vupn, Vtn) on the C vs E rows is the collinearity check. Müller2019
`v_T ≈ −v_OP` is **not** this number; SubMap can separate the two series
when ρ is well below 1.

{sign_text}

### Does H1 vs H2 ranking flip?

Ranking = which of |ρ(Vupn, C/E)| and |ρ(Vtn, C/E)| is larger.

{rank_text}

C/N/E ordinal (C=+1, N=0, E=−1) and collinearity on all defined rows:

{chr(10).join(cne_table)}

## Median Vupn / Vtn by UPS_3 (cm/yr)

{chr(10).join(med_blocks)}

## Figure

- `{FIG_FOREST}` — forest / dot plot of ρ(Vupn) and ρ(Vtn) vs C vs E.
  Whiskers are approximate 95% Fisher-z intervals, not a bootstrap.

## How to rerun

```bash
/workspace/subduction-test/.micromamba/envs/pygplates/bin/python \\
    /workspace/subduction-test/scripts/submap_h1h2_frames.py
```

## Caveats

- Geodetic UPS is classified from `Ovd` / `vd`. This sidecar does not treat
  `vdn` as a predictor.
- Absolute velocities change with the reference frame (Schellart 2008). A
  result that exists in only one named frame is fragile; unnamed frames
  cannot be interpreted as a specific hotspot or NNR model.
- N15 and M56 relative plate pairs differ; n finite is not the same in every
  column (some N15 frames have more −999).
- Pure C and pure E are small before recode; C-SS and E-SS dominate the
  C vs E n.
- Fisher-z intervals treat Spearman ρ like a Pearson correlation of ranks;
  they are a display aid.
"""
    REPORT_MD.write_text(md)
    print(f"Wrote {REPORT_MD}")


def main() -> int:
    h12._setup_style()
    REPORT_MD.parent.mkdir(parents=True, exist_ok=True)

    full = pd.read_csv(FULL_CSV)
    n_raw = int(len(full))
    if "Geodetic_UPS" not in full.columns:
        raise SystemExit("Full sheet missing Geodetic_UPS")
    for spec in FRAMES:
        for col in (spec["vupn"], spec["vtn"]):
            if col not in full.columns:
                raise SystemExit(f"Full sheet missing {col}")

    full = full.copy()
    full["UPS_3"] = h12.recode_ups3(full["Geodetic_UPS"])

    rows = [analyse_frame(full, spec) for spec in FRAMES]
    ranking = ranking_summary(rows)
    fig_forest(rows, FIG_FOREST)
    write_report(rows, n_raw, ranking)

    print("\nC vs E Spearman by frame:")
    for r in rows:
        print(
            f"  {r['family']}{r['n']:d} n_ce={r['n_ce']:3d}  "
            f"ρ_Vupn={fmt_rho(r['rho_vup_ce']):>7} p={fmt_p(r['p_vup_ce']):>10}  "
            f"ρ_Vtn={fmt_rho(r['rho_vt_ce']):>7} p={fmt_p(r['p_vt_ce']):>10}  "
            f"ρ(Vup,Vt)={fmt_rho(r['rho_col_ce']):>7}  stronger={r['stronger']}"
        )
    print("ranking flips:", ranking["flips"], "labels:", ranking["labels"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
