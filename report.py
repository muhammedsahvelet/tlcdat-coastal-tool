"""
report.py — one-page branded PDF reports for TLCDAT (matplotlib only, no extra deps).
Each pdf_*(...) returns PDF bytes ready for Streamlit's download_button.
"""
import io, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import coastal_core as cc

PUR = "#1F3864"; EDGE = "#2E6CA4"


# ----------------------------------------------------------------- page frame
def _kv(fig, x, y, title, items, dy=0.020, colw=0.20):
    fig.text(x, y, title, fontsize=11, fontweight="bold", color=PUR)
    yy = y - 0.026
    for k, v in items.items():
        fig.text(x, yy, str(k), fontsize=8.3, color="0.30")
        fig.text(x + colw, yy, str(v), fontsize=8.3, fontweight="bold", color="black")
        yy -= dy
    return yy


def build_pdf(module_title, inputs, results, draw, note):
    fig = plt.figure(figsize=(8.5, 11))
    fig.patch.set_facecolor("white")
    # header
    fig.text(0.06, 0.955, "TLCDAT", fontsize=22, fontweight="bold", color=PUR)
    fig.text(0.225, 0.9575, "Coastal Diffraction Analysis Tool", fontsize=9.5, color=EDGE)
    fig.text(0.06, 0.936, module_title, fontsize=13, fontweight="bold", color="black")
    fig.text(0.94, 0.958, "Troy Lab · Purdue University", ha="right", fontsize=8, color="0.35")
    fig.text(0.94, 0.944, "Lyles School of Civil & Construction Eng.", ha="right", fontsize=7.5, color="0.45")
    fig.text(0.94, 0.930, datetime.date.today().isoformat(), ha="right", fontsize=8, color="0.45")
    fig.add_artist(Line2D([0.06, 0.94], [0.905, 0.905], transform=fig.transFigure,
                          color=PUR, lw=1.6))
    # inputs (left) + results (right)
    _kv(fig, 0.06, 0.875, "Inputs", inputs)
    _kv(fig, 0.55, 0.875, "Results", results)
    # main plot region (callback adds axes)
    draw(fig)
    # footer
    fig.add_artist(Line2D([0.06, 0.94], [0.075, 0.075], transform=fig.transFigure,
                          color="0.7", lw=0.8))
    fig.text(0.06, 0.055, "Method: " + note, fontsize=7.3, color="0.25", wrap=True)
    fig.text(0.06, 0.022, "Generated with TLCDAT — M. N. Sahvelet and C. D. Troy, "
                          "Troy Lab, Purdue University.", fontsize=7, color=EDGE)
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    plt.close(fig)
    return buf.getvalue()


def _heat(fig, rect, xs, ys, Z, title, marks=None, beam=None, vmax=1.2):
    ax = fig.add_axes(rect)
    cf = ax.contourf(ys, xs, Z, levels=np.linspace(0, vmax, 25), cmap="jet", extend="max")
    cs = ax.contour(ys, xs, Z, levels=[.2, .4, .5, .6, .8, 1.0], colors="k", linewidths=0.3)
    ax.clabel(cs, fontsize=5, fmt="%.1f")
    if marks:
        for m in marks: ax.axvline(m, color="w", ls=":", lw=1.0)
    if beam is not None:
        ax.plot([0, beam * xs.max()], [0, xs.max()], "w--", lw=1.1)
    ax.set_title(title, fontsize=9); ax.set_xlabel("lateral y (m)", fontsize=8)
    ax.set_ylabel("distance into basin x (m)", fontsize=8)
    ax.tick_params(labelsize=7)
    fig.colorbar(cf, ax=ax, label="K = H/H₀", fraction=0.046, pad=0.02)


# ----------------------------------------------------------------- modules
def pdf_gap(T, d, B, theta0, beta_k, xmax, ymax):
    L = cc.wavelength(T, d); s = np.sin(np.radians(theta0))
    ymax_plot = max(ymax, xmax * abs(s) + B)
    xs, ys, Z, _ = cc.gap_map(T, d, B, xmax, ymax_plot, nx=110, ny=150, L=L,
                              theta0_deg=theta0, beta_k=beta_k)
    Ra, Ta = cc.barrier_RT(theta0, beta_k)
    inp = {"Wave period T": f"{T} s", "Water depth d": f"{d} m", "Gap width B": f"{B} m",
           "Angle θ₀": f"{theta0}°", "Permeability β/k": f"{beta_k:.2f}",
           "Basin x·max": f"{xmax} m", "Half-width y": f"{ymax} m"}
    res = {"Wavelength L": f"{L:.1f} m", "B/L": f"{B/L:.2f}",
           "Barrier |Ta|": f"{Ta:.2f}", "Barrier |Ra|": f"{Ra:.2f}",
           "K, beam @ x=2L": f"{cc.gap_K(2*L, 2*L*s, B, L, theta0, beta_k):.2f}",
           "K, centre @ x=5L": f"{cc.gap_K(5*L, 5*L*s, B, L, theta0, beta_k):.2f}"}

    def draw(fig):
        _heat(fig, [0.10, 0.40, 0.78, 0.28], xs, ys, Z,
              f"Gap K   (B/L={B/L:.1f}, θ₀={theta0}°, β/k={beta_k:.2f})",
              marks=[-B/2, B/2], beam=(s if theta0 != 0 else None))
        ax = fig.add_axes([0.10, 0.13, 0.80, 0.15])
        kc = [cc.gap_K(x, x*s, B, L, theta0, beta_k) for x in xs]
        ax.plot(xs, kc, color=PUR, lw=1.8); ax.grid(alpha=.3); ax.set_ylim(0, 1.4)
        ax.set_xlabel("x into basin (m)", fontsize=8); ax.set_ylabel("K (beam axis)", fontsize=8)
        ax.tick_params(labelsize=7); ax.set_title("K along the beam axis", fontsize=9)

    note = ("Fresnel/Kirchhoff single-slit diffraction (Penney & Price 1952; SPM Fig 2-44..2-52); "
            "oblique = paraxial beam-steering; permeable arms via Bowen & McIver (2002) β/k, "
            "K = √(|Ta|² + |Ra|²·K_solid²).")
    return build_pdf("Entrance Gap Diffraction Report", inp, res, draw, note)


def pdf_breakwater(T, d, B, ang, refl, xmax, ymax):
    L = cc.wavelength(T, d)
    xs, ys, Z = cc.detached_map(T, d, B, xmax, ymax, ang, refl, nx=70, ny=90)
    inp = {"Wave period T": f"{T} s", "Water depth d": f"{d} m", "Breakwater B": f"{B} m",
           "Obliquity θ": f"{ang}°", "Reflection": f"{refl:.2f}",
           "Basin x·max": f"{xmax} m", "Half-width y": f"{ymax} m"}
    res = {"Wavelength L": f"{L:.1f} m", "B/L": f"{B/L:.2f}",
           "K min": f"{Z.min():.2f}", "K max": f"{Z.max():.2f}",
           "K behind centre": f"{cc.end_diff_wiegel(T, d, B, xmax*0.5, ang, B/2, refl):.2f}"}

    def draw(fig):
        _heat(fig, [0.10, 0.40, 0.78, 0.28], xs, ys, Z,
              f"Detached breakwater K   (θ={ang}°)", marks=[-B/2, B/2])
        ax = fig.add_axes([0.10, 0.13, 0.80, 0.15])
        kc = [cc.end_diff_wiegel(T, d, B, x, ang, B/2, refl) for x in xs]
        ax.plot(xs, kc, color=PUR, lw=1.8); ax.grid(alpha=.3)
        ax.set_xlabel("x behind structure (m)", fontsize=8); ax.set_ylabel("K (centre)", fontsize=8)
        ax.tick_params(labelsize=7); ax.set_title("Centreline behind breakwater", fontsize=9)

    note = ("Wiegel (1962) K' diffraction table with coherent two-tip end-diffraction "
            "superposition (SPM guidance).")
    return build_pdf("Detached Breakwater Diffraction Report", inp, res, draw, note)


def pdf_dock(T, d, dock, barr, sep, ang, beam, draft, refl):
    L = cc.wavelength(T, d); Kt = cc.kt_kbc(T, d, beam, draft)
    pos = np.linspace(0, dock, 25)
    ed = [cc.end_diff_wiegel(T, d, barr, sep, ang, p, refl) for p in pos]
    comb = [np.sqrt(e*e + Kt*Kt) for e in ed]
    inp = {"Wave period T": f"{T} s", "Water depth d": f"{d} m", "Dock length": f"{dock} m",
           "Barrier length": f"{barr} m", "Barrier–dock": f"{sep} m", "Obliquity θ": f"{ang}°",
           "Beam width": f"{beam} m", "Draft": f"{draft} m"}
    res = {"Wavelength L": f"{L:.1f} m", "Transmission Kt": f"{Kt:.3f}",
           "Barrier B/L": f"{barr/L:.2f}", "Combined max": f"{max(comb):.2f}",
           "Combined min": f"{min(comb):.2f}"}

    def draw(fig):
        ax = fig.add_axes([0.12, 0.16, 0.78, 0.50])
        ax.plot(pos/L, ed, color=EDGE, marker="o", ms=3, label="End diffraction")
        ax.axhline(Kt, color="crimson", ls="--", label=f"Transmission Kt={Kt:.2f}")
        ax.plot(pos/L, comb, color="green", marker="s", ms=3, label="Combined")
        ax.set_xlabel("position along dock (×L)", fontsize=9); ax.set_ylabel("coefficient", fontsize=9)
        ax.set_title("End diffraction + transmission along the dock", fontsize=10)
        ax.legend(fontsize=8); ax.grid(alpha=.3); ax.set_ylim(0, max(1.1, max(comb)*1.1))
        ax.tick_params(labelsize=8)

    note = ("Kriebel/Cox transmission (ASCE Manual 50, eq 2-32) combined with Wiegel end "
            "diffraction: Combined = √(diffraction² + transmission²).")
    return build_pdf("Dock Response Report", inp, res, draw, note)


def pdf_refraction(T, d, a0, H0):
    r = cc.refraction_shoaling(T, d, a0, H0)
    inp = {"Wave period T": f"{T} s", "Local depth d": f"{d} m",
           "Deep angle α₀": f"{a0}°", "Deep height H₀": f"{H0} m"}
    res = {"Wavelength L": f"{r['L']:.1f} m", "Local angle α": f"{r['alpha_deg']:.1f}°",
           "Refraction Kr": f"{r['Kr']:.3f}", "Shoaling Ks": f"{r['Ks']:.3f}",
           "Wave height H": f"{r['H']:.3f} m"}

    def draw(fig):
        ax = fig.add_axes([0.12, 0.24, 0.76, 0.42])
        ds = np.linspace(max(0.3, d/8), max(d*2, 20), 60)
        ax.plot(ds, [cc.shoaling_Ks(T, dd)[0] for dd in ds], color=EDGE, lw=2, label="Ks")
        ax.axvline(d, color="crimson", ls="--", lw=1); ax.plot([d], [r["Ks"]], "ro")
        ax.set_xlabel("water depth d (m)", fontsize=9); ax.set_ylabel("shoaling Ks", fontsize=9)
        ax.set_title(f"Shoaling curve (T={T}s) — your case marked", fontsize=10)
        ax.grid(alpha=.3); ax.legend(fontsize=8); ax.tick_params(labelsize=8)

    note = ("Snell's law for straight, parallel offshore contours (sin α = (L/L₀) sin α₀); "
            "shoaling Ks from SPM eq. 2-44; H = H₀·Kr·Ks.")
    return build_pdf("Refraction + Shoaling Report", inp, res, draw, note)


def pdf_wavelength(T, d):
    methods = cc.wavelength_methods(T, d)
    inp = {"Wave period T": f"{T} s", "Water depth d": f"{d} m",
           "Deep-water L₀": f"{cc.deep_L0(T):.2f} m"}
    res = {k: f"{v[0]:.2f} m  ({v[1]:+.2f}%)" for k, v in methods.items()}

    def draw(fig):
        ax = fig.add_axes([0.12, 0.24, 0.76, 0.42])
        ds = np.linspace(1, max(40, d*1.5), 60)
        for nm, f in [("Fenton", cc.wl_fenton), ("Guo", cc.wl_guo),
                      ("You", cc.wl_you), ("Hunt", cc.wl_hunt)]:
            ax.plot(ds, [100*(f(T, dd)-cc.wavelength(T, dd))/cc.wavelength(T, dd) for dd in ds], label=nm)
        ax.axvline(d, color="crimson", ls="--", lw=1)
        ax.set_xlabel("water depth d (m)", fontsize=9); ax.set_ylabel("error vs exact (%)", fontsize=9)
        ax.set_title(f"Wavelength approximation error (T={T}s)", fontsize=10)
        ax.grid(alpha=.3); ax.legend(fontsize=8); ax.tick_params(labelsize=8)

    note = ("Exact linear dispersion L = L₀·tanh(2πd/L) by iteration; Hunt (1979), "
            "Fenton & McKee (1990), Guo (2002), You (2008) explicit approximations.")
    return build_pdf("Wavelength / Dispersion Report", inp, res, draw, note)
