# -*- coding: utf-8 -*-
"""Mimari preprint figurleri — data/mimari_hesap.json -> figs/*.pdf+png."""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle

HERE = Path(__file__).resolve().parent
D = json.loads((HERE.parent / "results" / "mimari_hesap.json").read_text(encoding="utf-8"))
FIGS = HERE / "figs"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 9, "axes.titlesize": 9,
    "legend.fontsize": 7.0, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "figure.dpi": 150, "savefig.bbox": "tight"})

C1 = "#1f77b4"; C2 = "#d62728"; C3 = "#2ca02c"; C4 = "#9467bd"
HBAR = 0.65821195


def kaydet(fig, ad):
    fig.savefig(FIGS / f"{ad}.pdf")
    fig.savefig(FIGS / f"{ad}.png", dpi=300)
    plt.close(fig)
    print("yazildi:", ad)


# ---------------------------------------------------------------- fig 1
def fig1_mimari():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.2),
                               gridspec_kw={"width_ratios": [1.15, 1.0]})
    # (a) d=3 cell: 3x3 data dual-rails + 8 ancilla, checkerboard
    a.set_xlim(-1.9, 6.1); a.set_ylim(-1.9, 5.9)
    a.set_aspect("equal"); a.axis("off")
    a.set_title("(a) $d=3$ XZZX cell: 26 traps")
    # rotated-code checkerboard: 4 bulk + 4 boundary ancillas
    anc_xy = [(1, 1), (3, 1), (1, 3), (3, 3), (2, 5.0), (2, -1.0),
              (-1.0, 2), (5.0, 2)]
    for (x, y) in anc_xy:
        a.add_patch(Circle((x, y), 0.22, fc="#f6d7d7", ec=C2, lw=1.1))
        for (dx, dy) in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
            xx, yy = x + dx, y + dy
            if -0.1 <= xx <= 4.1 and -0.1 <= yy <= 4.1:
                a.plot([x, xx], [y, yy], color="0.8", lw=0.6, zorder=0)
    for i in range(3):
        for j in range(3):
            x, y = 2 * j, 2 * (2 - i)
            a.add_patch(FancyBboxPatch((x - 0.52, y - 0.28), 1.04, 0.56,
                                       boxstyle="round,pad=0.03",
                                       fc="#dbe9f6", ec=C1, lw=1.0))
            a.add_patch(Circle((x - 0.26, y), 0.15, fc=C1, ec="none"))
            a.add_patch(Circle((x + 0.26, y), 0.15, fc="white", ec=C1, lw=1.1))
    a.plot([], [], "s", ms=9, mfc="#dbe9f6", mec=C1, label="data dual-rail "
           "($|L\\rangle$/$|R\\rangle$)")
    a.plot([], [], "o", ms=8, mfc="#f6d7d7", mec=C2, label="ancilla (single trap)")
    a.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=1,
             fontsize=6.6, frameon=False)
    a.text(2, -1.65, "tunnel-contact row (ns reset)  |  crossbar gates",
           ha="center", fontsize=6.4, color="0.35")
    a.text(2, 5.6, "moiré lattice = registration grid "
           "($\\mathrm{d}a_M/\\mathrm{d}\\theta = 0$)",
           ha="center", fontsize=6.4, color="0.35")

    # (b) layer/clock diagram
    b.axis("off"); b.set_xlim(0, 10); b.set_ylim(0, 10)
    b.set_title("(b) protection stack and the two clocks")
    katman = [
        ("L0  passive: dark shelf + dual-rail DFS + quadrupole park", "#eef5ec", C3),
        ("L1  autonomous loop (4--10 ns, measurement-free):\n"
         "      syndrome $\\to$ coherent Z-fix $\\to$ contact dump $\\to$ fresh IX",
         "#fdf3e3", "#e09a2f"),
        ("L2  biased-erasure XZZX, $d=3$", "#f6e3e3", C2),
        ("L3  slow clock: flash/click every $K\\geq 25$ rounds\n"
         "      (updates decoder frame only)", "#e8ecf6", C1),
    ]
    y = 9.0
    for et, fc, ec in katman:
        n = et.count("\n") + 1
        h = 1.0 + 0.85 * (n - 1)
        b.add_patch(FancyBboxPatch((0.3, y - h), 9.4, h,
                                   boxstyle="round,pad=0.08", fc=fc, ec=ec,
                                   lw=1.1))
        b.text(0.65, y - h / 2, et, va="center", fontsize=6.8)
        y -= h + 0.55
    b.annotate("", xy=(9.75, 1.2), xytext=(9.75, 8.9),
               arrowprops=dict(arrowstyle="->", lw=0.8, color="0.4"))
    b.text(9.55, 5.0, "entropy out", rotation=90, fontsize=6.2, color="0.4",
           va="center")
    fig.tight_layout()
    kaydet(fig, "fig1_mimari")


# ---------------------------------------------------------------- fig 2
def fig2_tunel():
    h1 = D["H1_cift_tuzak"]
    V0 = np.array(h1["V0_meV"]); A = np.array(h1["a_nm"])
    T = np.array(h1["t_ueV"])    # [V0, a]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9))
    pc = a.pcolormesh(A, V0, np.log10(T), shading="nearest", cmap="viridis",
                      vmin=-7, vmax=2.2)
    cb = fig.colorbar(pc, ax=a)
    cb.set_label("$\\log_{10} t$ ($\\mu$eV)", fontsize=8)
    cs = a.contour(A, V0, np.log10(T), levels=[math.log10(0.05),
                                               math.log10(0.3)],
                   colors=["w", "w"], linewidths=1.2,
                   linestyles=["-", "--"])
    a.clabel(cs, fmt={math.log10(0.05): "t=0.05 $\\mu$eV",
                      math.log10(0.3): "0.3 $\\mu$eV"}, fontsize=6.5)
    a.axhline(200, color="w", ls=":", lw=1.0)
    a.text(10.3, 207, "$V_0=200$ meV accessibility line", color="w",
           fontsize=6.2)
    a.plot([20], [150], "w*", ms=10)
    a.text(19.0, 157, "v2 OFF", color="w", fontsize=6.8, ha="right")
    a.set_xlabel("trap separation $a$ (nm)")
    a.set_ylabel("trap depth $V_0$ (meV)")
    a.set_title("(a) OFF-state tunneling $t(V_0,a)$")

    xk = h1["X_kapisi"]
    dV = -np.array(xk["dV_meV"]); tX = np.array(xk["tX_ps"], float)
    gap = np.array(xk["E_gap_meV"])
    b.semilogy(dV, tX, "o-", color=C1, ms=4, lw=1.2)
    b.axhspan(20, 200, color=C3, alpha=0.15)
    b.text(12, 52, "target X-gate band\n20--200 ps", fontsize=6.5, color=C3)
    b.set_xlabel("barrier-gate dip $|\\Delta V_c|$ (meV)")
    b.set_ylabel("rail-SWAP time $t_X$ (ps)", color=C1)
    b.set_ylim(0.3, 1e7)
    b2 = b.twinx()
    b2.plot(dV, gap, "s--", color=C2, ms=3.5, lw=1.0)
    b2.set_ylabel("orbital gap $E_2-E_1$ (meV)", color=C2, fontsize=8)
    b2.set_ylim(0, 23)
    b.set_title("(b) barrier-controlled X gate")
    fig.tight_layout()
    kaydet(fig, "fig2_tunel")


# ---------------------------------------------------------------- fig 3
def fig3_ekran():
    h2 = D["H2_ekran"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9))
    stil = {"bare": ("k", "-", "unscreened"),
            "dm=20": (C1, "--", "$d_m=20$ nm"),
            "dm=10": (C3, "-.", "$d_m=10$ nm"),
            "dm=5": (C2, ":", "$d_m=5$ nm")}
    for k, (clr, ls, et) in stil.items():
        r = np.array(h2["V_r"][k]["r_nm"], float)
        V = np.array(h2["V_r"][k]["V_ueV"], float)
        a.loglog(r, V, color=clr, ls=ls, lw=1.2, label=et)
    a.axvline(20, color="0.7", lw=0.7); a.axvline(40, color="0.7", lw=0.7)
    a.text(20, 42, "$a$", fontsize=7.5, ha="center", color="0.4")
    a.text(40, 42, "$R$", fontsize=7.5, ha="center", color="0.4")
    a.set_xticks([10, 20, 30, 40, 60, 80])
    a.set_xticklabels(["10", "20", "30", "40", "60", "80"])
    a.tick_params(axis="x", which="minor", labelbottom=False)
    a.set_xlabel("lateral separation $r$ (nm)")
    a.set_ylabel("$V_{zz}(r)$ ($\\mu$eV)")
    a.set_title("(a) single-gate image correction to $V(r)$")
    a.legend(); a.grid(alpha=0.25, which="both")

    konum = {"a=10": 0, "a=20": 1, "a=35": 2}
    w = 0.2
    for i, (k, (clr, ls, et)) in enumerate(stil.items()):
        key = k
        xs, ys = [], []
        for ak, pos in konum.items():
            ys.append(h2["J_eff"][key][ak]["J_ueV"])
            xs.append(pos + (i - 1.5) * w)
        b.bar(xs, ys, width=w * 0.92, color=clr, alpha=0.85, label=et)
    b.axhline(4.0, color="k", ls="--", lw=0.8)
    b.text(1.62, 4.6, "v2 requirement $J\\geq 4\\,\\mu$eV", fontsize=6.5)
    b.set_yscale("log")
    b.set_xticks(list(konum.values()))
    b.set_xticklabels(["$a=10$ nm\n($R=20$)", "$a=20$ nm\n($R=40$)",
                       "$a=35$ nm\n($R=70$)"], fontsize=7.5)
    b.set_ylabel("$J_{\\rm eff}$ ($\\mu$eV)")
    b.set_title("(b) conditional coupling with gate screening")
    b.legend(fontsize=6.2)
    b.grid(alpha=0.25, axis="y", which="both")
    fig.tight_layout()
    kaydet(fig, "fig3_ekran")


# ---------------------------------------------------------------- fig 4
def fig4_cz():
    h4 = D["H4a_zz_cz"]["satirlar"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9))
    tc_stil = {"1 ns": (C2, "o"), "10 ns": (C3, "s"), "1 us": (C1, "^")}
    for tc, (clr, mk) in tc_stil.items():
        xs = [r["sigma_eff_ueV"] for r in h4 if r["eko"] and r["tau_c"] == tc]
        ys = [r["hata"] for r in h4 if r["eko"] and r["tau_c"] == tc]
        o = np.argsort(xs)
        a.loglog(np.array(xs)[o], np.array(ys)[o], marker=mk, ms=4, lw=1.1,
                 color=clr, label="echo, $\\tau_c=$ " + tc)
    xs = [r["sigma_eff_ueV"] for r in h4 if (not r["eko"])
          and r["tau_c"] == "1 us"]
    ys = [r["hata"] for r in h4 if (not r["eko"]) and r["tau_c"] == "1 us"]
    o = np.argsort(xs)
    a.loglog(np.array(xs)[o], np.array(ys)[o], "k--", marker="x", ms=4,
             lw=1.0, label="no echo (quasi-static)")
    a.axhline(1e-3, color="0.4", ls=":", lw=0.9)
    a.text(0.42, 1.35e-3, "$10^{-3}$ target", fontsize=6.5, color="0.3")
    # sigma_E bands
    for sE, lo, hi, yy in [("$\\sigma_E=3.7$", 0.4, 1.2, 3.3e-6),
                           ("10", 1.0, 3.3, 1.3e-6),
                           ("47 $\\mu$eV", 4.7, 15.7, 3.3e-6)]:
        a.annotate("", xy=(lo, yy), xytext=(hi, yy),
                   arrowprops=dict(arrowstyle="<->", lw=0.7, color="0.5"))
        a.text(math.sqrt(lo * hi), yy * 1.5, sE, fontsize=6.0, ha="center",
               color="0.4")
    a.set_xlabel("$\\sigma_{\\rm eff}$ ($\\mu$eV)  [$=\\sigma_E\\,"
                 "\\Delta x/\\lambda_c$]")
    a.set_ylabel("CZ error $\\varepsilon$")
    a.set_ylim(5e-7, 1.2)
    a.set_title("(a) echoed polarize-and-phase CZ (MC)")
    a.legend(fontsize=6.0, loc="upper left")
    a.grid(alpha=0.25, which="both")

    hb = D["H4b_blokaj_cz"]["satirlar"]
    for V, clr, mk in ((59.2, C1, "o"), (473.0, C2, "s")):
        xs = [max(r["sigma_ueV"], 0.2) for r in hb if r["V_ueV"] == V]
        ys = [r["hata"] for r in hb if r["V_ueV"] == V]
        es = [r.get("hata_std", 0) for r in hb if r["V_ueV"] == V]
        o = np.argsort(xs)
        b.errorbar(np.array(xs)[o], np.array(ys)[o], yerr=np.array(es)[o],
                   marker=mk, ms=4, lw=1.1, color=clr, capsize=2,
                   label=f"$V={V:g}\\,\\mu$eV ($r={10 if V<100 else 5}$ nm)")
    b.set_xscale("log"); b.set_yscale("log")
    b.axhline(1e-2, color="0.4", ls=":", lw=0.9)
    b.text(0.25, 1.3e-2, "percent class", fontsize=6.5, color="0.3")
    b.set_xlabel("residual site detuning $\\sigma$ ($\\mu$eV)\n"
                 "(0.2 $\\approx$ trimmed floor; 47 = untrimmed @ "
                 "$0.1^\\circ$)", fontsize=7.5)
    b.set_ylabel("blockade-CZ error")
    b.set_title("(b) Jaksch blockade variant (3-level Lindblad)")
    b.legend(fontsize=6.2, loc="lower right")
    b.grid(alpha=0.25, which="both")
    fig.tight_layout()
    kaydet(fig, "fig4_cz")


# ---------------------------------------------------------------- fig 5
def fig5_silinme():
    sil = D["silinme"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.8))
    for Tet, clr, mk in (("1.5K", C1, "o"), ("10K", C2, "s")):
        xs = [p["eta"] for p in sil["pencere"] if p["T"] == Tet]
        ys = [p["t_pencere_ns"] for p in sil["pencere"] if p["T"] == Tet]
        a.loglog(xs, ys, marker=mk, ms=4.5, lw=1.2, color=clr,
                 label=f"$T=${Tet} ($\\beta=${sil['beta'][Tet]})")
    a.set_xticks([0.05, 0.1, 0.2])
    a.set_xticklabels(["0.05", "0.1", "0.2"])
    a.tick_params(axis="x", which="minor", labelbottom=False)
    a.set_xlabel("collection efficiency $\\eta$")
    a.set_ylabel("flash window $t_w$ (ns)")
    a.set_title("(a) occupancy-flash window\n"
                "($T_1'=2.2$ ns, $\\eta_{\\rm SNSPD}=0.8$, err$<10^{-2}$)",
                fontsize=8)
    a.legend(); a.grid(alpha=0.25, which="both")

    raf = np.geomspace(0.3, 30, 200)   # us
    for tw, clr, ls in ((129, C3, "-"), (257, C1, "--"), (514, C2, ":")):
        b.loglog(raf, tw * 1e-3 / (raf * 25) * 100, color=clr, ls=ls, lw=1.2,
                 label=f"$t_w={tw}$ ns")
    b.axhline(0.5, color="k", ls="--", lw=0.8)
    b.text(0.35, 0.57, "0.5%/round budget", fontsize=6.5)
    b.axvspan(1, 10, color="0.85", alpha=0.5)
    b.text(1.45, 0.028, "shelf design\nband 1--10 $\\mu$s", fontsize=6.3,
           ha="center", color="0.35")
    b.set_xlabel("dark-shelf lifetime $T_{\\rm shelf}$ ($\\mu$s)")
    b.set_ylabel("amortized flash loss (%/round), $K=25$")
    b.set_title("(b) slow-clock amortization", fontsize=8.5)
    b.legend(fontsize=6.3)
    b.grid(alpha=0.25, which="both")
    fig.tight_layout()
    kaydet(fig, "fig5_silinme")


if __name__ == "__main__":
    fig1_mimari(); fig2_tunel(); fig3_ekran(); fig4_cz(); fig5_silinme()
    print("tamam")
