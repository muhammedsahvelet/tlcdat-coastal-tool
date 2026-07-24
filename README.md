# TLCDAT — Troy Lab Coastal Diffraction Analysis Tool (web edition)

A free, coastal-engineering diffractio calculator.


**Five modules**
1. **Wavelength / dispersion** — exact iteration + Hunt, Fenton, Guo, You approximations.
2. **Refraction + shoaling** — Snell's law (straight, parallel contours) + Ks, Kr.
3. **Breakwater diffraction** — Wiegel K′ table + coherent two-tip end superposition.
4. **Entrance gap** — Fresnel / Kirchhoff single-slit diffraction (SPM Fig 2-44…2-52), continuous in gap width, **normal & oblique incidence** (the beam steers with the incoming wave angle).
5. **Dock transmission** — Kriebel/Cox transmission + combined dock response.

Every module has a one-click **PDF report** button (branded header, inputs, results and plots).

```
coastal_core.py   the engineering engine (pure Python, importable, unit-tested)
app.py            the Streamlit user interface (5 tabs)
report.py         one-page branded PDF report generator (matplotlib, no extra deps)
requirements.txt  dependencies
run_windows.bat   one-click launch on Windows
run_mac.command   one-click launch on macOS
```

---

## 1. Run it on your own computer (2 minutes)

You need **Python 3.9+** (get it from python.org — tick *"Add Python to PATH"* on Windows).

**Windows:** double-click **`run_windows.bat`**.
**macOS:** double-click **`run_mac.command`** (first time: right-click → Open).

Or from a terminal in this folder:
```bash
pip install -r requirements.txt
streamlit run app.py
```
Your browser opens at `http://localhost:8501`. The engineer only ever sees the web page —
they never touch code.

---

## 2. Put it online for the whole office (free, recommended)

So colleagues just click a link and nothing is installed on their machines:

1. Create a free **GitHub** account and upload this folder to a new repository.
2. Go to **share.streamlit.io** (Streamlit Community Cloud), sign in with GitHub.
3. **New app** → pick the repo → main file `app.py` → **Deploy**.
4. You get a public URL like `https://your-tool.streamlit.app`. Share it with Jack and the team.

Cost: **$0**. Updates: push to GitHub and the site refreshes automatically.

---

## 3. Offline / desktop `.exe` (same code, no internet)

The same `app.py` can become a **downloadable desktop app that runs with no server and no
internet** using **stlite** (Streamlit compiled to WebAssembly):

- **Single offline web page:** wrap `app.py` with `@stlite/mountable` → one `index.html` that
  runs entirely in the browser. Host it free on **GitHub Pages** (no server at all).
- **Windows `.exe` / macOS `.app`:** use **`@stlite/desktop`** (Electron wrapper) to package the
  same page into an installer engineers double-click. No Python required on their machine.

(Traditional **PyInstaller** also works via a small launcher that runs `streamlit run app.py`
and opens the browser, but stlite gives the cleaner, license-free, no-install result.)

---

## Method & references
- Penney & Price (1952) — gap / optical-analogy diffraction
- Wiegel (1962) — semi-infinite breakwater K′ table
- Shore Protection Manual (1984) §2 — refraction, shoaling, gap diagrams (Fig 2-44…2-52)
- Kriebel & Cox — barrier transmission (ASCE Manual 50, eq 2-32)
- Bowen & McIver (2002) — permeable-gap large-gap approximation
- Abramowitz & Stegun (1964) §7.3 — Fresnel integrals

Prepared by **M. N. Sahvelet & C. D. Troy** (Troy Lab, Purdue University, Lyles School of Civil &
Construction Engineering) 
