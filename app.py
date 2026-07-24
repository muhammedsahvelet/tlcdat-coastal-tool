"""
TLCDAT — Troy Lab Coastal Diffraction Analysis Tool (web edition)
A free, no-code coastal-engineering calculator.  Streamlit front-end over
coastal_core.py.  Run locally with:   streamlit run app.py
Deploy free on Streamlit Community Cloud (see README).

Purdue University · Lyles School of Civil & Construction Engineering 
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st
import coastal_core as cc
import report

@st.cache_data(show_spinner=False)
def make_pdf(kind, *args):
    """Cached PDF generation — only rebuilds when inputs change."""
    return getattr(report, "pdf_" + kind)(*args)

st.set_page_config(page_title="TLCDAT — Coastal Diffraction Tool",
                   page_icon="🌊", layout="wide")

# ----------------------------------------------------------------- styling
PURDUE = "#1F3864"; EDGE = "#2E6CA4"
st.markdown(f"""
<style>
.block-container {{padding-top: 1.4rem;}}
h1, h2, h3 {{color: {PURDUE};}}
.stTabs [data-baseweb="tab-list"] {{gap: 4px;}}
.stTabs [data-baseweb="tab"] {{padding: 8px 14px;}}
</style>""", unsafe_allow_html=True)

st.title("🌊 TLCDAT — Coastal Diffraction Analysis Tool")
st.caption("Troy Lab · Purdue University (Lyles School of Civil & Construction Engineering) · "
           "free wave-transformation calculators for coastal design.")

def heatmap(xs, ys, Z, title, xlabel="lateral y (m)", ylabel="into basin x (m)",
            marks=None, vmax=1.2, beam_slope=None):
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    cf = ax.contourf(ys, xs, Z, levels=np.linspace(0, vmax, 25), cmap="jet", extend="max")
    cs = ax.contour(ys, xs, Z, levels=[.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],
                    colors="k", linewidths=0.35)
    ax.clabel(cs, fontsize=6, fmt="%.1f")
    if marks:
        for m in marks: ax.axvline(m, color="w", ls=":", lw=1.1)
    if beam_slope is not None:                      # geometric beam axis (oblique)
        ax.plot([0, beam_slope*xs.max()], [0, xs.max()], color="w", ls="--", lw=1.2)
    ax.set_title(title, fontsize=10); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    fig.colorbar(cf, ax=ax, label="K = H/H₀", fraction=0.046, pad=0.02)
    fig.tight_layout()
    return fig

tabs = st.tabs(["📏 Wavelength", "↩️ Refraction + Shoaling",
                "🧱 Breakwater diffraction", "🚪 Entrance gap", "⚓ Dock transmission",
                "ℹ️ Method & notes"])

# ============================================================= 1. WAVELENGTH
with tabs[0]:
    st.subheader("Wavelength / dispersion")
    c1, c2 = st.columns([1, 2])
    with c1:
        T = st.number_input("Wave period T (s)", 0.5, 60.0, 6.0, 0.1, key="wl_T")
        d = st.number_input("Water depth d (m)", 0.1, 500.0, 3.0, 0.5, key="wl_d")
    with c2:
        res = cc.wavelength_methods(T, d)
        rows = [{"Method": k, "L (m)": round(v[0], 3), "error vs exact (%)": round(v[1], 3)}
                for k, v in res.items()]
        st.dataframe(rows, width="stretch", hide_index=True)
        st.metric("Deep-water wavelength L₀ (m)", f"{cc.deep_L0(T):.2f}")
        st.info("The **Iteration** value is the exact linear-theory solution "
                "(L = L₀·tanh(2πd/L)). Hunt, Fenton, Guo and You are explicit approximations — "
                "all within ~1.7 %. Use any of them where a closed form is convenient.")
        st.download_button("📄 Download PDF report", data=make_pdf("wavelength", T, d),
                           file_name="TLCDAT_wavelength_report.pdf", mime="application/pdf",
                           key="dl_wl")

# ============================================================= 2. REFRACTION + SHOALING
with tabs[1]:
    st.subheader("Refraction + shoaling (straight, parallel contours — Snell's law)")
    c1, c2 = st.columns([1, 2])
    with c1:
        T = st.number_input("Wave period T (s)", 0.5, 60.0, 8.0, 0.1, key="rf_T")
        d = st.number_input("Local water depth d (m)", 0.1, 500.0, 10.0, 0.5, key="rf_d")
        a0 = st.number_input("Deep-water wave angle to shore-normal α₀ (deg)",
                             0.0, 89.0, 20.0, 1.0, key="rf_a")
        H0 = st.number_input("Deep-water wave height H₀ (m)", 0.0, 30.0, 1.0, 0.1, key="rf_H")
    r = cc.refraction_shoaling(T, d, a0, H0)
    with c2:
        m = st.columns(4)
        m[0].metric("Wavelength L (m)", f"{r['L']:.2f}")
        m[1].metric("Refraction Kr", f"{r['Kr']:.3f}")
        m[2].metric("Shoaling Ks", f"{r['Ks']:.3f}")
        m[3].metric("Wave height H (m)", f"{r['H']:.3f}")
        st.write(f"Local wave angle to shore-normal **α = {r['alpha_deg']:.1f}°** "
                 f"(refracted from {a0:.0f}°). H = H₀·Kr·Ks.")
        # simple refraction sketch: incoming ray (deep) vs refracted ray (local)
        fig, ax = plt.subplots(figsize=(6.5, 3.0))
        ax.axhline(0, color="0.6", lw=6, alpha=.4)             # shoreline
        ax.annotate("", xy=(0, 0), xytext=(-np.sin(np.radians(a0))*1.2, np.cos(np.radians(a0))*1.2),
                    arrowprops=dict(arrowstyle="->", color=EDGE, lw=2))
        ax.annotate("", xy=(np.sin(np.radians(r["alpha_deg"]))*1.2, np.cos(np.radians(r["alpha_deg"]))*1.2),
                    xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="crimson", lw=2))
        ax.plot([0, 0], [0, 1.35], color="0.5", ls=":", lw=1)  # shore normal
        ax.text(-np.sin(np.radians(a0))*1.25, np.cos(np.radians(a0))*1.25,
                f"deep α₀={a0:.0f}°", color=EDGE, fontsize=8, ha="right")
        ax.text(np.sin(np.radians(r["alpha_deg"]))*1.25, np.cos(np.radians(r["alpha_deg"]))*1.25,
                f"local α={r['alpha_deg']:.0f}°", color="crimson", fontsize=8)
        ax.text(0.03, 1.30, "shore-normal", color="0.4", fontsize=7)
        ax.text(1.0, -0.16, "shoreline", color="0.4", fontsize=8, ha="right")
        ax.set_xlim(-1.4, 1.4); ax.set_ylim(-0.3, 1.5); ax.axis("off")
        ax.set_title("Wave crest turns toward the shore-normal as depth decreases", fontsize=9)
        st.pyplot(fig)
        st.download_button("📄 Download PDF report", data=make_pdf("refraction", T, d, a0, H0),
                           file_name="TLCDAT_refraction_report.pdf", mime="application/pdf",
                           key="dl_rf")

# ============================================================= 3. BREAKWATER DIFFRACTION
with tabs[2]:
    st.subheader("Detached / semi-infinite breakwater diffraction (Wiegel + end superposition)")
    c1, c2 = st.columns([1, 2])
    with c1:
        T = st.number_input("Wave period T (s)", 0.5, 60.0, 4.0, 0.1, key="bw_T")
        d = st.number_input("Water depth d (m)", 0.1, 500.0, 30.0, 0.5, key="bw_d")
        B = st.number_input("Breakwater length B (m)", 1.0, 5000.0, 30.0, 1.0, key="bw_B")
        ang = st.number_input("Wave obliquity θ (deg, 90 = normal)", 15.0, 180.0, 90.0, 5.0, key="bw_a")
        refl = st.number_input("Reflection factor (0.7–1.0)", 0.0, 1.0, 1.0, 0.05, key="bw_r")
        xmax = st.number_input("Basin extent x (m)", 5.0, 5000.0, 60.0, 5.0, key="bw_x")
        ymax = st.number_input("Half-width y (m)", 5.0, 5000.0, 45.0, 5.0, key="bw_y")
    with c2:
        xs, ys, Z = cc.detached_map(T, d, B, xmax, ymax, ang, refl, nx=70, ny=90)
        st.pyplot(heatmap(xs, ys, Z, f"Detached breakwater K   (T={T}s, d={d}m, B={B}m, θ={ang:.0f}°)",
                          marks=[-B/2, B/2]))
        kc = [cc.end_diff_wiegel(T, d, B, x, ang, B/2, refl) for x in xs]
        fig, ax = plt.subplots(figsize=(7.4, 2.4))
        ax.plot(xs, kc, color=PURDUE, lw=2); ax.grid(alpha=.3)
        ax.set_xlabel("distance behind structure x (m)"); ax.set_ylabel("K (centre)")
        ax.set_title("Centreline behind breakwater", fontsize=9)
        st.pyplot(fig)
        st.caption("Two breakwater tips, coherent (phase-aware) superposition of Wiegel end "
                   "diffraction — this is Jack's original method, reproduced macro-free.")
        st.download_button("📄 Download PDF report",
                           data=make_pdf("breakwater", T, d, B, ang, refl, xmax, ymax),
                           file_name="TLCDAT_breakwater_report.pdf", mime="application/pdf",
                           key="dl_bw")

# ============================================================= 4. GAP
with tabs[3]:
    st.subheader("Entrance gap diffraction (Fresnel single-slit — SPM Fig 2-44…2-52)")
    c1, c2 = st.columns([1, 2])
    with c1:
        T = st.number_input("Wave period T (s)", 0.5, 60.0, 6.0, 0.1, key="gp_T")
        d = st.number_input("Water depth d (m)", 0.1, 500.0, 8.0, 0.5, key="gp_d")
        B = st.number_input("Gap width B (m)", 1.0, 5000.0, 120.0, 5.0, key="gp_B")
        theta0 = st.slider("Wave angle from normal θ₀ (deg, 0 = head-on)",
                           -60, 60, 0, 5, key="gp_th")
        beta_k = st.slider("Breakwater permeability β/k  (0 = solid wall)",
                           0.0, 1.5, 0.0, 0.05, key="gp_bk")
        xmax = st.number_input("Basin extent x (m)", 10.0, 8000.0, 300.0, 10.0, key="gp_x")
        ymax = st.number_input("Half-width y (m)", 10.0, 8000.0, 250.0, 10.0, key="gp_y")
        L = cc.wavelength(T, d)
        st.metric("Wavelength L (m)", f"{L:.1f}")
        st.metric("Gap width B/L", f"{B/L:.2f}")
        if beta_k > 0:
            Ra, Ta = cc.barrier_RT(theta0, beta_k)
            mm = st.columns(2)
            mm[0].metric("Barrier transmission |Ta|", f"{Ta:.2f}")
            mm[1].metric("Barrier reflection |Ra|", f"{Ra:.2f}")
        if B/L < 1:      st.warning("Narrow gap (B/L<1): radial spreading — estimate only.")
        elif B/L < 5:    st.success("Slit-diffraction range (1 ≤ B/L < 5).")
        else:            st.info("Wide gap (B/L ≥ 5): edges act almost independently.")
        if abs(theta0) > 40:
            st.warning("Large obliquity: paraxial steering is approximate beyond ~40°.")
        st.markdown("**Point check**")
        px = st.number_input("x point (m)", 0.1, 8000.0, 100.0, 5.0, key="gp_px")
        py = st.number_input("y point (m)", -8000.0, 8000.0, 0.0, 5.0, key="gp_py")
        st.metric("K at point", f"{cc.gap_K(px, py, B, L, theta0, beta_k):.3f}")
    with c2:
        s = np.sin(np.radians(theta0))
        ymax_plot = max(ymax, xmax*abs(s) + B)          # keep the steered beam in frame
        xs, ys, Z, L = cc.gap_map(T, d, B, xmax, ymax_plot, nx=110, ny=170, L=L,
                                  theta0_deg=theta0, beta_k=beta_k)
        ttl = (f"Gap K   (T={T}s, d={d}m, B={B}m, B/L={B/L:.1f}, "
               f"θ₀={theta0}°, β/k={beta_k:.2f})")
        st.pyplot(heatmap(xs, ys, Z, ttl, marks=[-B/2, B/2], vmax=1.2,
                          beam_slope=(s if theta0 != 0 else None)))
        # K along the geometric beam axis  y = x*sin(theta0)
        kc = [cc.gap_K(x, x*s, B, L, theta0, beta_k) for x in xs]
        fig, ax = plt.subplots(figsize=(7.4, 2.4))
        ax.plot(xs, kc, color=PURDUE, lw=2); ax.grid(alpha=.3); ax.set_ylim(0, 1.4)
        ax.set_xlabel("distance into basin x (m)"); ax.set_ylabel("K (beam axis)")
        ax.set_title("K along the beam axis" + ("" if theta0 == 0 else f" (steered {theta0}°)"),
                     fontsize=9)
        st.pyplot(fig)
        if theta0 != 0:
            st.caption("White dashed line = geometric beam axis (y = x·sin θ₀). "
                       "The transmitted beam steers toward the incoming wave direction.")
        if beta_k != 0:
            st.caption("Permeable arms (Bowen & McIver 2002, β/k): waves also leak straight "
                       "through the structure, so the sheltered lee no longer goes to zero — "
                       "deep-shadow K rises toward the barrier transmission |Ta|.")
        st.download_button("📄 Download PDF report",
                           data=make_pdf("gap", T, d, B, theta0, beta_k, xmax, ymax),
                           file_name="TLCDAT_gap_report.pdf", mime="application/pdf",
                           key="dl_gp")

# ============================================================= 5. DOCK
with tabs[4]:
    st.subheader("Dock response: end diffraction + barrier transmission")
    c1, c2 = st.columns([1, 2])
    with c1:
        T = st.number_input("Wave period T (s)", 0.5, 60.0, 3.0, 0.1, key="dk_T")
        d = st.number_input("Water depth d (m)", 0.1, 500.0, 3.0, 0.5, key="dk_d")
        dock = st.number_input("Dock length (m)", 1.0, 2000.0, 18.0, 1.0, key="dk_dock")
        barr = st.number_input("Barrier length (m)", 1.0, 2000.0, 18.0, 1.0, key="dk_barr")
        sep = st.number_input("Barrier–dock distance (m)", 0.1, 2000.0, 6.0, 0.5, key="dk_sep")
        ang = st.number_input("Obliquity θ (deg)", 15.0, 180.0, 45.0, 5.0, key="dk_a")
        beam = st.number_input("Barrier beam width (m)", 0.1, 200.0, 3.0, 0.1, key="dk_beam")
        draft = st.number_input("Barrier draft (m, =d for full depth)", 0.1, 500.0, 2.0, 0.1, key="dk_draft")
        refl = st.number_input("Reflection factor", 0.0, 1.0, 1.0, 0.05, key="dk_r")
    with c2:
        Kt = cc.kt_kbc(T, d, beam, draft)
        pos = np.linspace(0, dock, 25)
        ed = [cc.end_diff_wiegel(T, d, barr, sep, ang, p, refl) for p in pos]
        comb = [np.sqrt(e**2 + Kt**2) for e in ed]  # end diffraction + transmission energy
        L = cc.wavelength(T, d)
        fig, ax = plt.subplots(figsize=(7.6, 4.2))
        ax.plot(pos/L, ed, label="End diffraction", color=EDGE, marker="o", ms=3)
        ax.axhline(Kt, color="crimson", ls="--", label=f"Transmission Kt={Kt:.2f}")
        ax.plot(pos/L, comb, label="Combined", color="green", marker="s", ms=3)
        ax.set_xlabel("position along dock  (×L)"); ax.set_ylabel("coefficient")
        ax.set_title("End diffraction + transmission along the dock", fontsize=10)
        ax.legend(); ax.grid(alpha=.3); ax.set_ylim(0, max(1.1, max(comb)*1.1))
        st.pyplot(fig)
        cc1 = st.columns(3)
        cc1[0].metric("Wavelength L (m)", f"{L:.1f}")
        cc1[1].metric("Transmission Kt", f"{Kt:.3f}")
        cc1[2].metric("Barrier B/L", f"{barr/L:.2f}")
        st.caption("Combined = √(diffraction² + transmission²): energy reaching the dock past the "
                   "barrier by wrapping around its ends and by passing through it.")
        st.download_button("📄 Download PDF report",
                           data=make_pdf("dock", T, d, dock, barr, sep, ang, beam, draft, refl),
                           file_name="TLCDAT_dock_report.pdf", mime="application/pdf",
                           key="dl_dk")

# ============================================================= 6. NOTES
with tabs[5]:
    st.subheader("Method & references")
    st.markdown("""
**What this is.** A free Troy-Lab coastal spreadsheet, plus a
new continuous **entrance-gap** calculator. Everything runs from one clean Python core  no license.

**Methods**
- *Wavelength* — exact linear dispersion by iteration; Hunt (1979), Fenton & McKee (1990),
  Guo (2002) and You (2008) explicit approximations for comparison.
- *Refraction + shoaling* — Snell's law for straight, parallel contours; shoaling Ks from SPM eq. 2-44.
- *Breakwater diffraction* — Wiegel (1962) K′ table with coherent two-tip end-diffraction
  superposition (the original workbook's method).
- *Entrance gap* — Fresnel / Kirchhoff single-slit diffraction, the optical analogy of
  Penney & Price (1952) that underlies the SPM gap diagrams. Continuous in gap width,
  with **oblique incidence** (paraxial beam steering) and **permeable / partially-reflecting
  breakwater arms** via the Bowen & McIver (2002) β/k reflection–transmission (Ra, Ta):
  K = √(|Ta|² + |Ra|²·K_solid²).
- *Dock transmission* — Kriebel/Cox transmission (ASCE Manual 50, eq. 2-32).

**Assumptions** — constant depth, monochromatic linear waves, diffraction only
(no breaking/friction). Gap: normal **and oblique** incidence (paraxial beam-steering by x·sin θ₀,
most accurate for angles ≲40°) and **solid or permeable** breakwater arms (β/k; real β is
lossless, |Ra|²+|Ta|²=1). Breakwater/dock modules assume rigid arms.

**References** — Penney & Price (1952); Wiegel (1962); SPM (1984) §2; Kriebel & Cox;
Bowen & McIver (2002); Abramowitz & Stegun (1964) §7.3.
""")
    st.caption("Prepared by M. N. Sahvelet & C. D. Troy (Troy Lab, Purdue) ")
