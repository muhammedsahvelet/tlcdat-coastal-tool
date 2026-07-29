"""
cad_export.py — CAD-compatible export for TLCDAT diffraction fields.

Turns a diffraction-coefficient field K(x,y) into engineering deliverables that
drop straight into AutoCAD / Civil 3D and GIS:

  * DXF  (AutoCAD R12 ASCII) — labelled K contour polylines on colour-coded
    layers, plus the breakwater / gap geometry, at real-world scale (metres).
    Import into CAD and overlay on the site plan; the gap centre is at (0,0)
    and the structure line lies along the X-axis (Y = 0).
  * CSV / XYZ grid — one row per grid point (x, y, K) for GIS or surface import.

Pure Python, no third-party dependency (uses only matplotlib, already required,
to trace the contour paths).  Authors: M. N. Sahvelet & C. D. Troy, Troy Lab, Purdue.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# K level -> AutoCAD Colour Index (blue = sheltered, red = amplified)
def _aci(level):
    if level <= 0.2: return 5   # blue
    if level <= 0.4: return 4   # cyan
    if level <= 0.6: return 3   # green
    if level <= 0.8: return 2   # yellow
    return 1                    # red

DEFAULT_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


# --------------------------------------------------------------------- DXF core
class _DXF:
    def __init__(self):
        self.e = []            # entity group-code list
        self.layers = {}       # name -> color
    def _g(self, code, val):   # append one group code/value pair
        self.e.append(str(code)); self.e.append(str(val))
    def layer(self, name, color):
        self.layers[name] = color
    def line(self, layer, x1, y1, x2, y2, color=None):
        self._g(0, "LINE"); self._g(8, layer)
        if color is not None: self._g(62, color)
        self._g(10, f"{x1:.4f}"); self._g(20, f"{y1:.4f}"); self._g(30, "0.0")
        self._g(11, f"{x2:.4f}"); self._g(21, f"{y2:.4f}"); self._g(31, "0.0")
    def text(self, layer, x, y, h, s, color=None):
        self._g(0, "TEXT"); self._g(8, layer)
        if color is not None: self._g(62, color)
        self._g(10, f"{x:.4f}"); self._g(20, f"{y:.4f}"); self._g(30, "0.0")
        self._g(40, f"{h:.4f}"); self._g(1, str(s))
    def polyline(self, layer, pts, color=None):
        self._g(0, "POLYLINE"); self._g(8, layer)
        if color is not None: self._g(62, color)
        self._g(66, 1); self._g(70, 0)
        self._g(10, "0.0"); self._g(20, "0.0"); self._g(30, "0.0")
        for (x, y) in pts:
            self._g(0, "VERTEX"); self._g(8, layer)
            self._g(10, f"{x:.4f}"); self._g(20, f"{y:.4f}"); self._g(30, "0.0")
        self._g(0, "SEQEND"); self._g(8, layer)
    def build(self):
        out = ["0", "SECTION", "2", "HEADER",
               "9", "$ACADVER", "1", "AC1009",
               "9", "$INSUNITS", "70", "6",          # 6 = metres
               "0", "ENDSEC",
               "0", "SECTION", "2", "TABLES",
               "0", "TABLE", "2", "LAYER", "70", str(len(self.layers))]
        for name, col in self.layers.items():
            out += ["0", "LAYER", "2", name, "70", "0", "62", str(col), "6", "CONTINUOUS"]
        out += ["0", "ENDTAB", "0", "ENDSEC",
                "0", "SECTION", "2", "ENTITIES"]
        out += self.e
        out += ["0", "ENDSEC", "0", "EOF"]
        return ("\r\n".join(out) + "\r\n").encode("ascii", "replace")


def kd_dxf(xs, ys, Z, B, title="TLCDAT diffraction K field", levels=None, kind="gap"):
    """Build a DXF (bytes) from a K field.
    xs = into-basin distance (m), ys = lateral offset (m), Z[len(xs)][len(ys)].
    kind = "gap" (solid arms outside +/-B/2, open gap in the middle) or
           "breakwater" (a solid detached breakwater of length B centred at origin).
    Plan orientation: CAD-X = lateral y, CAD-Y = into-basin x; origin at the gap /
    breakwater centre, structure line along Y = 0."""
    if levels is None: levels = DEFAULT_LEVELS
    Z = np.asarray(Z); xs = np.asarray(xs); ys = np.asarray(ys)
    txt = max(0.02*max(xs.max(), ys.max()-ys.min()), 0.5)   # text height ~ 2% of domain

    dxf = _DXF()
    for lv in levels: dxf.layer(f"Kd_{lv:.2f}", _aci(lv))
    dxf.layer("BREAKWATER", 8); dxf.layer("GAP", 1); dxf.layer("NOTES", 7)

    # contour polylines (CAD-X = ys, CAD-Y = xs)
    fig = plt.figure()
    cs = plt.contour(ys, xs, Z, levels=levels)
    for i, lv in enumerate(cs.levels):
        layer = f"Kd_{lv:.2f}"; col = _aci(lv)
        segs = cs.allsegs[i] if i < len(cs.allsegs) else []
        for seg in segs:
            if len(seg) < 2: continue
            dxf.polyline(layer, seg, col)
            mid = seg[len(seg)//2]
            dxf.text(layer, mid[0], mid[1], txt*0.8, f"{lv:.1f}", col)
    plt.close(fig)

    # structure along Y = 0
    y0, y1 = float(ys.min()), float(ys.max())
    tick = txt*1.5
    if kind == "breakwater":
        dxf.line("BREAKWATER", -B/2, 0.0, B/2, 0.0, 8)          # solid detached breakwater
        dxf.line("GAP", -B/2, 0.0, -B/2, -tick, 1)
        dxf.line("GAP",  B/2, 0.0,  B/2, -tick, 1)
        dxf.text("GAP", 0.0, -tick*1.6, txt, f"BREAKWATER  L={B:g} m", 1)
    else:
        dxf.line("BREAKWATER", y0, 0.0, -B/2, 0.0, 8)           # solid arms, open gap
        dxf.line("BREAKWATER", B/2, 0.0, y1, 0.0, 8)
        dxf.line("GAP", -B/2, 0.0, -B/2, -tick, 1)
        dxf.line("GAP",  B/2, 0.0,  B/2, -tick, 1)
        dxf.text("GAP", 0.0, -tick*1.6, txt, f"GAP  B={B:g} m", 1)
    # domain frame
    dxf.line("NOTES", y0, 0.0, y0, float(xs.max()), 7)
    dxf.line("NOTES", y1, 0.0, y1, float(xs.max()), 7)
    dxf.line("NOTES", y0, float(xs.max()), y1, float(xs.max()), 7)
    # title + orientation note
    dxf.text("NOTES", y0, float(xs.max())+txt*1.6, txt*1.2, title, 7)
    dxf.text("NOTES", y0, float(xs.max())+txt*0.2, txt*0.8,
             "X=lateral(m)  Y=into basin(m)  origin=gap centre  units=m  K=H/H0", 7)
    return dxf.build()


def kd_csv(xs, ys, Z):
    """Flatten the K field to CSV bytes: x_into_basin_m, y_lateral_m, K."""
    xs = np.asarray(xs); ys = np.asarray(ys); Z = np.asarray(Z)
    lines = ["x_into_basin_m,y_lateral_m,K"]
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            lines.append(f"{x:.3f},{y:.3f},{Z[i][j]:.4f}")
    return ("\n".join(lines) + "\n").encode("ascii")
