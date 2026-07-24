"""
coastal_core.py  —  Troy Lab Coastal Diffraction Analysis Tool (TLCDAT)
Pure-Python engineering core.  Ported faithfully from


Modules:
  1. Dispersion / wavelength   (Newton + Hunt, Fenton, Guo, You approximations)
  2. Refraction + shoaling      (Snell's law, Ks, Kr  — USACE Vol. 2 §2)
  3. Wiegel diffraction         (K' table + End_Diff superposition = detached breakwater)
  4. Dock transmission          (Kriebel/Cox Kt + combined dock response)
  5. Entrance gap diffraction   (Fresnel / Kirchhoff single-slit — SPM Fig 2-44..2-52)

Authors: M. N. Sahvelet & C. D. Troy Ph.D. (Troy Lab, Purdue University).
"""
import numpy as np

G = 9.81

# ============================================================ 1. DISPERSION
def deep_L0(T):
    return G * T**2 / (2*np.pi)

def wavelength(T, d, tol=1e-9, itmax=200):
    """Exact linear-theory wavelength by fixed-point iteration L = L0*tanh(2*pi*d/L).
    Matches the VBA Td_L() result to machine precision."""
    L0 = deep_L0(T); L = L0
    for _ in range(itmax):
        Ln = L0 * np.tanh(2*np.pi*d/L)
        if abs(Ln - L) < tol:
            return Ln
        L = Ln
    return L

def _kd_from_L(L, d):
    return 2*np.pi*d/L

def wl_hunt(T, d):
    """Hunt (1979) 9th-order approximation for (kd)^2."""
    y = (2*np.pi/T)**2 * d / G                      # sigma^2 d / g
    dd = [0.6666666666, 0.3555555555, 0.1608465608,
          0.0632098765, 0.0217540484, 0.0065407983]
    denom = 1.0 + sum(dn*y**(n+1) for n, dn in enumerate(dd))
    kd = np.sqrt(y*y + y/denom)
    return 2*np.pi*d/kd

def wl_fenton(T, d):
    """Fenton & McKee (1990): kd = (s2d/g)*[coth((s*sqrt(d/g))^1.5)]^(2/3)."""
    s = 2*np.pi/T
    a = (s*np.sqrt(d/G))**1.5
    kd = (s*s*d/G) * (1.0/np.tanh(a))**(2.0/3.0)
    return 2*np.pi*d/kd

def wl_guo(T, d):
    """Guo (2002): kd = (s2d/g)*(1 - exp(-(s*sqrt(d/g))^2.5))^(-2/5)."""
    s = 2*np.pi/T
    x = s*np.sqrt(d/G)
    kd = (s*s*d/G) * (1.0 - np.exp(-(x**2.5)))**(-2.0/5.0)
    return 2*np.pi*d/kd

def wl_you(T, d):
    """You (2008): L = L0*tanh(xi0), xi0=(k0 d)^0.5*(1 + k0d/6 + (k0d)^2/30)."""
    L0 = deep_L0(T); k0 = 2*np.pi/L0
    k0d = k0*d
    xi0 = np.sqrt(k0d) * (1 + k0d/6 + k0d**2/30)
    return L0*np.tanh(xi0)

def wavelength_methods(T, d):
    """Return every method + % error vs the exact iteration."""
    exact = wavelength(T, d)
    out = {"Iteration (exact)": (exact, 0.0)}
    for name, f in [("Hunt (1979)", wl_hunt), ("Fenton (1990)", wl_fenton),
                    ("Guo (2002)", wl_guo), ("You (2008)", wl_you)]:
        try:
            L = f(T, d); err = 100*(L-exact)/exact
        except Exception:
            L, err = float("nan"), float("nan")
        out[name] = (L, err)
    return out

# ============================================================ 2. REFRACTION + SHOALING
def shoaling_Ks(T, d):
    """Shoaling coefficient Ks = H/H0'  (SPM eq 2-44)."""
    L = wavelength(T, d); kd = 2*np.pi*d/L
    n = 0.5*(1 + (2*kd)/np.sinh(2*kd))          # 4*pi*d/L = 2*kd
    Ks = np.sqrt(1.0/np.tanh(kd) * 1.0/(2*n))
    return Ks, n, L

def refraction_shoaling(T, d, alpha0_deg, H0=1.0):
    """Straight, parallel-contour refraction (Snell) + shoaling.
    alpha0_deg = deep-water wave angle measured from the shore normal (deg).
    Returns dict with L, local wave angle, Kr, Ks, H."""
    L0 = deep_L0(T); L = wavelength(T, d)
    a0 = np.radians(alpha0_deg)
    ratio = L/L0                                  # C/C0 = L/L0
    sin_a = np.clip(ratio*np.sin(a0), -1, 1)
    a = np.arcsin(sin_a)
    Kr = np.sqrt(max(np.cos(a0), 1e-9)/max(np.cos(a), 1e-9))
    Ks, n, _ = shoaling_Ks(T, d)
    return {"L0": L0, "L": L, "alpha0_deg": alpha0_deg,
            "alpha_deg": np.degrees(a), "Kr": Kr, "Ks": Ks, "n": n,
            "H": H0*Kr*Ks}

# ============================================================ 3. WIEGEL DIFFRACTION
_BETA = np.array([0,15,30,45,60,75,90,105,120,135,150,165,180.0])
_RL   = np.array([0,0.5,1,2,5,10.0])
_THETA= np.array([15,30,45,60,75,90,105,120,135,150,165,180.0])
# Wiegel K' table lifted verbatim from the VBA (theta -> 6 rows[r/L] x 13 cols[beta])
_W = {
15:[[0.8]*13,[0.49,0.79,0.83,0.9,0.97,1.01,1.03,1.02,1.01,0.99,0.99,1,1],[0.38,0.73,0.83,0.95,1.04,1.04,0.99,0.98,1.01,1.01,1,1,1],[0.21,0.68,0.86,1.05,1.03,0.97,1.02,0.99,1,1,1,1,1],[0.13,0.63,0.99,1.04,1.03,1.02,0.99,0.99,1,1.01,1,1,1],[0.35,0.58,1.1,1.05,0.98,0.99,1.01,1,1,1,1,1,1]],
30:[[0.7]*13,[0.61,0.63,0.68,0.76,0.87,0.97,1.03,1.05,1.03,1.01,0.99,0.95,1],[0.5,0.53,0.63,0.78,0.95,1.06,1.05,0.98,0.98,1.01,1.01,0.97,1],[0.4,0.44,0.59,0.84,1.07,1.03,0.96,1.02,0.98,1.01,0.99,0.95,1],[0.27,0.32,0.55,1,1.04,1.04,1.02,0.99,0.99,1,1.01,0.97,1],[0.2,0.24,0.54,1.12,1.06,0.97,0.99,1.01,1,1,1,0.98,1]],
45:[[0.65]*13,[0.49,0.5,0.55,0.63,0.73,0.85,0.96,1.04,1.06,1.04,1,0.99,1],[0.38,0.4,0.47,0.59,0.76,0.95,1.07,1.06,0.98,0.97,1.01,1.01,1],[0.29,0.31,0.39,0.56,0.83,1.08,1.04,0.96,1.03,0.98,1.01,1,1],[0.18,0.2,0.29,0.54,1.01,1.04,1.05,1.03,1,0.99,1.01,1,1],[0.13,0.15,0.22,0.53,1.13,1.07,0.96,0.98,1.02,0.99,1,1,1]],
60:[[0.6]*13,[0.4,0.41,0.45,0.52,0.6,0.72,0.85,1.13,1.04,1.06,1.03,1.01,1],[0.31,0.32,0.36,0.44,0.57,0.75,0.96,1.08,1.06,0.98,0.98,1.01,1],[0.22,0.23,0.28,0.37,0.55,0.83,1.08,1.04,0.96,1.03,0.98,1.01,1],[0.14,0.15,0.18,0.28,0.53,1.01,1.04,1.05,1.03,0.99,0.99,1,1],[0.1,0.11,0.13,0.21,0.52,1.14,1.07,0.96,0.98,1.01,1,1,1]],
75:[[0.55]*13,[0.34,0.35,0.38,0.42,0.5,0.59,0.71,0.85,0.97,1.04,1.05,1.02,1],[0.25,0.26,0.29,0.34,0.43,0.56,0.75,0.95,1.02,1.06,0.98,0.98,1],[0.18,0.19,0.22,0.26,0.36,0.54,0.83,1.09,1.04,0.96,1.03,0.99,1],[0.12,0.12,0.13,0.17,0.27,0.52,1.01,1.04,1.05,1.03,0.99,0.99,1],[0.08,0.08,0.1,0.13,0.2,0.52,1.14,1.07,0.96,0.98,1.01,1,1]],
90:[[0.52]*13,[0.31,0.31,0.33,0.36,0.41,0.49,0.59,0.71,0.85,0.96,1.03,1.03,1],[0.22,0.23,0.24,0.28,0.33,0.42,0.56,0.75,0.96,1.07,1.05,0.99,1],[0.16,0.16,0.18,0.2,0.26,0.35,0.54,0.69,1.08,1.04,0.96,1.02,1],[0.1,0.1,0.11,0.13,0.16,0.27,0.53,1.01,1.04,1.05,1.02,0.99,1],[0.07,0.07,0.08,0.09,0.13,0.2,0.52,1.14,1.07,0.96,0.99,1.01,1]],
105:[[0.55]*13,[0.28,0.28,0.29,0.32,0.35,0.41,0.49,0.59,0.72,0.85,0.97,1.01,1],[0.2,0.2,0.24,0.23,0.27,0.33,0.42,0.56,0.75,0.95,1.06,1.04,1],[0.14,0.14,0.13,0.17,0.2,0.25,0.35,0.54,0.83,1.08,1.03,0.97,1],[0.09,0.09,0.1,0.11,0.13,0.17,0.27,0.52,1.02,1.04,1.04,1.02,1],[0.07,0.06,0.08,0.08,0.09,0.12,0.2,0.52,1.14,1.07,0.97,0.99,1]],
120:[[0.6]*13,[0.25,0.26,0.27,0.28,0.31,0.35,0.41,0.5,0.6,0.73,0.87,0.97,1],[0.18,0.19,0.19,0.21,0.23,0.27,0.33,0.43,0.57,0.76,0.95,1.04,1],[0.13,0.13,0.14,0.14,0.17,0.2,0.26,0.16,0.55,0.83,1.07,1.03,1],[0.08,0.08,0.08,0.09,0.11,0.13,0.16,0.27,0.53,1.01,1.04,1.03,1],[0.06,0.06,0.06,0.07,0.07,0.09,0.13,0.2,0.52,1.13,1.06,0.98,1]],
135:[[0.6]*13,[0.24,0.24,0.25,0.26,0.28,0.32,0.36,0.42,0.52,0.63,0.76,0.9,1],[0.18,0.17,0.18,0.19,0.21,0.23,0.28,0.34,0.44,0.59,0.78,0.95,1],[0.12,0.12,0.13,0.14,0.14,0.17,0.2,0.26,0.37,0.56,0.84,1.05,1],[0.08,0.07,0.08,0.08,0.09,0.11,0.13,0.17,0.28,0.54,1,1.04,1],[0.05,0.06,0.06,0.06,0.07,0.08,0.09,0.13,0.21,0.53,1.12,1.05,1]],
150:[[0.5]*13,[0.23,0.23,0.24,0.25,0.27,0.29,0.33,0.38,0.45,0.55,0.68,0.83,1],[0.16,0.17,0.17,0.18,0.19,0.22,0.24,0.29,0.36,0.47,0.63,0.83,1],[0.12,0.12,0.12,0.13,0.14,0.15,0.18,0.22,0.28,0.39,0.59,0.86,1],[0.07,0.07,0.08,0.08,0.08,0.1,0.11,0.13,0.18,0.29,0.55,0.99,1],[0.05,0.05,0.05,0.06,0.06,0.07,0.08,0.1,0.13,0.22,0.54,1.1,1]],
165:[[0.5]*13,[0.23,0.23,0.23,0.24,0.26,0.28,0.31,0.35,0.41,0.5,0.63,0.79,1],[0.16,0.16,0.17,0.17,0.19,0.2,0.23,0.26,0.32,0.4,0.53,0.73,1],[0.11,0.11,0.12,0.12,0.13,0.14,0.16,0.19,0.23,0.31,0.44,0.68,1],[0.07,0.07,0.07,0.07,0.08,0.09,0.1,0.12,0.15,0.2,0.32,0.63,1],[0.05,0.05,0.05,0.06,0.06,0.06,0.07,0.08,0.11,0.11,0.21,0.58,1]],
180:[[0.5]*13,[0.2,0.25,0.23,0.24,0.25,0.28,0.31,0.34,0.4,0.49,0.61,0.78,1],[0.1,0.17,0.16,0.18,0.18,0.23,0.22,0.25,0.31,0.38,0.5,0.7,1],[0.02,0.09,0.12,0.12,0.13,0.18,0.16,0.18,0.22,0.29,0.4,0.6,1],[0.02,0.06,0.07,0.07,0.07,0.08,0.1,0.12,0.14,0.18,0.27,0.46,1],[0.1,0.05,0.05,0.04,0.06,0.07,0.07,0.08,0.1,0.13,0.2,0.36,1]],
}
_WARR = {k: np.array(v) for k, v in _W.items()}

def _interp(xs, ys, x0):
    """Linear interp with end-segment extrapolation (matches VBA Interpolate_arr)."""
    xs = np.asarray(xs, float); ys = np.asarray(ys, float); n = len(xs)
    if x0 <= xs[0]: a, b = 0, 1
    elif x0 > xs[-1]: a, b = n-2, n-1
    else:
        b = int(np.searchsorted(xs, x0)); a = b-1
    if xs[b] == xs[a]: return ys[a]
    return ys[a] + (ys[b]-ys[a])*(x0-xs[a])/(xs[b]-xs[a])

def kd_wiegel(theta, rL, beta):
    """Wiegel K' for one semi-infinite breakwater tip (trilinear on the table)."""
    kb = np.empty(len(_THETA))
    for ti, th in enumerate(_THETA):
        rows = _WARR[int(th)]
        over_beta = [ _interp(_BETA, rows[ri], beta) for ri in range(6) ]
        kb[ti] = _interp(_RL, over_beta, rL)
    return _interp(_THETA, kb, theta)

def end_diff_wiegel(T, d, attn, dist, angle, loc, reflec=1.0):
    """Detached-breakwater end diffraction: coherent superposition of the two
    breakwater tips. Faithful port of VBA End_Diff_Wiegel (r/L handled explicitly)."""
    L = wavelength(T, d)
    r_l = np.hypot(loc, dist)
    beta_l = 0.0 if r_l==0 else (angle if dist>=r_l else np.degrees(np.arcsin(min(dist/r_l,1))))
    RL = kd_wiegel(angle, r_l/L, beta_l)
    ph_l = (r_l/L - np.floor(r_l/L))*360
    r_r = np.hypot(attn-loc, dist)
    beta_r = 0.0 if r_r==0 else ((180-angle) if dist>=r_r else np.degrees(np.arcsin(min(dist/r_r,1))))
    RR = kd_wiegel(180-angle, r_r/L, beta_r)
    ph_r = (r_r/L - np.floor(r_r/L))*360
    ph = abs(ph_l - ph_r)
    return np.sqrt(RL*RL + RR*RR + 2*RL*RR*np.cos(np.radians(ph)))*reflec

# ============================================================ 4. DOCK TRANSMISSION
def kt_kbc(T, d, beam, draft):
    """Kriebel/Cox transmission coefficient (ASCE Manual 50, eq 2-32)."""
    L = wavelength(T, d)
    x1 = 4*np.pi*(d-draft)/L; x2 = 4*np.pi*d/L
    p = (x1 + np.sinh(x1))/(x2 + np.sinh(x2))
    Kt = 2*p/(1+p)
    y = 2*np.pi*beam/L
    Cb = (2*np.sqrt(1+y*y))/(2+y*y)
    return Cb*Kt

# ============================================================ 5. GAP (FRESNEL SLIT)
def fresnel_CS(v):
    """Fresnel cosine/sine integrals via Abramowitz & Stegun 7.3.32/7.3.33."""
    s = np.sign(v); x = np.abs(v)
    a = 0.5*np.pi*x*x
    f = (1+0.926*x)/(2+1.792*x+3.104*x*x)
    g = 1.0/(2+4.142*x+3.492*x*x+6.670*x**3)
    C = 0.5 + f*np.sin(a) - g*np.cos(a)
    S = 0.5 - f*np.cos(a) - g*np.sin(a)
    return s*C, s*S

def barrier_RT(theta0_deg, beta_k):
    """Bowen & McIver (2002) infinite permeable-barrier reflection |Ra| and
    transmission |Ta|, for permeability parameter beta_k = beta/k.
        Ra = -i k cos(theta0) / (2*beta - i k cos(theta0)),   Ta = 1 - Ra
    beta_k = 0  -> solid impermeable wall (|Ra|=1, |Ta|=0)
    beta_k -> inf -> fully transparent (|Ra|=0, |Ta|=1)
    For real beta_k the barrier is lossless: |Ra|^2 + |Ta|^2 = 1.
    theta0_deg is measured from the gap normal (0 = head-on)."""
    c = abs(np.cos(np.radians(theta0_deg)))
    denom = np.sqrt(4.0*beta_k*beta_k + c*c)
    if denom == 0:
        return 1.0, 0.0
    return c/denom, 2.0*beta_k/denom              # |Ra|, |Ta|

def gap_K(x, y, B, L, theta0_deg=0.0, beta_k=0.0):
    """Diffraction coefficient K behind a breakwater gap of width B.
    x = distance into basin (>0), y = lateral from gap centre.
    theta0_deg = wave angle from the gap normal (0 = head-on; +ve steers the
    transmitted beam toward +y). Fresnel single-slit with paraxial beam steering.
    beta_k = beta/k permeability of the breakwater arms (Bowen & McIver 2002);
    0 = solid/fully-reflecting (reduces exactly to the classic slit solution),
    larger = more wave energy leaks through the structure into the lee."""
    x = max(x, 1e-6); k = np.sqrt(2.0/(L*x))
    yc = y - x*np.sin(np.radians(theta0_deg))          # beam-steering shift
    v1 = (yc + B/2)*k; v2 = (yc - B/2)*k
    C1, S1 = fresnel_CS(v1); C2, S2 = fresnel_CS(v2)
    reC = C1-C2; imC = S1-S2
    ur = (reC+imC)/2; ui = (imC-reC)/2
    Kg = np.hypot(ur, ui)                              # solid-barrier gap diffraction
    if beta_k <= 0:
        return Kg
    Ra, Ta = barrier_RT(theta0_deg, beta_k)
    # lee field = uniform barrier transmission + diffraction scaled by reflection
    return np.hypot(Ta, Ra*Kg)

def gap_map(T, d, B, xmax, ymax, nx=120, ny=160, L=None, theta0_deg=0.0, beta_k=0.0):
    if L is None: L = wavelength(T, d)
    xs = np.linspace(xmax/40, xmax, nx); ys = np.linspace(-ymax, ymax, ny)
    Z = np.array([[gap_K(x, y, B, L, theta0_deg, beta_k) for y in ys] for x in xs])
    return xs, ys, Z, L

def detached_map(T, d, B, xmax, ymax, angle=90, reflec=1.0, nx=90, ny=120):
    """Detached-breakwater K map using the Wiegel superposition (his method)."""
    xs = np.linspace(0.5, xmax, nx); ys = np.linspace(-ymax, ymax, ny)
    Z = np.array([[end_diff_wiegel(T, d, B, x, angle, y+B/2, reflec) for y in ys] for x in xs])
    return xs, ys, Z
