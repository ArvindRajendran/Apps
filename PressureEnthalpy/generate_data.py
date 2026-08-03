#!/usr/bin/env python3
"""Generate embedded property data for PressureEnthalpy.html.

Uses CoolProp (reference EOS: Span-Wagner CO2, IAPWS-95 water,
Buecker-Wagner ethane; default = NIST WebBook reference states).
Injects JSON between /*__D__*/ and /*__E__*/ markers in the HTML.

Run:  python3 generate_data.py
"""
import base64
import json
import math
import re
import struct

import CoolProp
from CoolProp.CoolProp import PropsSI

ATM = 101325.0
NP, NH = 120, 180          # crosshair grid resolution (log10 P_atm x h)
NCRV = 130                 # samples per iso-curve

FLUIDS = [
    dict(key="CO2", cp="CO2", label="CO₂ (carbon dioxide)",
         Pmin=1.0*ATM, Pmax=300.0*ATM, Tmin=223.15, Tmax=473.15,   # -50..200 C
         Tlevels=list(range(-40, 201, 20)), sStep=0.25),
    dict(key="H2O", cp="Water", label="H₂O (water / steam)",
         Pmin=0.01*ATM, Pmax=300.0*ATM, Tmin=274.15, Tmax=973.15,  # 1..700 C
         Tlevels=list(range(50, 701, 50)), sStep=1.0),
    dict(key="C2H6", cp="Ethane", label="C₂H₆ (ethane)",
         Pmin=0.1*ATM, Pmax=100.0*ATM, Tmin=93.15, Tmax=473.15,    # -180..200 C
         Tlevels=list(range(-180, 201, 20)), sStep=0.5),
]
RHO_CAND = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50,
            100, 200, 300, 500, 700, 900, 1100]

NAN = float("nan")


def logspace(a, b, n):
    la, lb = math.log10(a), math.log10(b)
    return [10 ** (la + (lb - la) * i / (n - 1)) for i in range(n)]


def linspace(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def q16(vals):
    """Quantize list (row-major) to uint16 base64 with min/max scaling."""
    fin = [v for v in vals if v == v]
    lo, hi = (min(fin), max(fin)) if fin else (0.0, 1.0)
    if hi <= lo:
        hi = lo + 1.0
    sc = 65534.0 / (hi - lo)
    words = [65535 if v != v else max(0, min(65534, int(round((v - lo) * sc))))
             for v in vals]
    raw = struct.pack("<%dH" % len(words), *words)
    return dict(min=lo, max=hi, b64=base64.b64encode(raw).decode())


def gen_fluid(cfg):
    cp = cfg["cp"]
    AS = CoolProp.AbstractState("HEOS", cp)
    Tc = PropsSI("Tcrit", cp)
    Pc = PropsSI("pcrit", cp)
    Ttr = PropsSI("Ttriple", cp)
    Ptr = PropsSI("ptriple", cp)
    Pmin, Pmax, Tmin, Tmax = cfg["Pmin"], cfg["Pmax"], cfg["Tmin"], cfg["Tmax"]
    la = lambda P: math.log10(P / ATM)

    def flash(pair, a, b):
        try:
            AS.update(pair, a, b)
            T = AS.T()
            if T < Tmin - 10 or T > Tmax + 10:
                return None
            # exclude pseudo-liquid/solid states below the triple pressure
            if AS.p() < Ptr and AS.rhomass() > 100:
                return None
            return (T - 273.15, AS.smass() / 1000.0, AS.rhomass(),
                    AS.hmass() / 1000.0)
        except Exception:
            return None

    # ---- saturation dome (dense near critical) ----
    dome = dict(la=[], hl=[], hv=[], ts=[], sl=[], sv=[])
    Tlo = max(Ttr, Tmin)
    n = 150
    for i in range(n):
        u = i / (n - 1)
        T = Tc - (Tc - Tlo) * (1 - u) ** 2
        if i == n - 1:
            T = Tc - 1e-4 * (Tc - Tlo)
        try:
            P = PropsSI("P", "T", T, "Q", 0, cp)
            hl = PropsSI("H", "T", T, "Q", 0, cp) / 1000.0
            hv = PropsSI("H", "T", T, "Q", 1, cp) / 1000.0
            sl = PropsSI("S", "T", T, "Q", 0, cp) / 1000.0
            sv = PropsSI("S", "T", T, "Q", 1, cp) / 1000.0
        except Exception:
            continue
        if P < Pmin or P > Pmax:
            continue
        dome["la"].append(round(la(P), 5))
        dome["hl"].append(round(hl, 2))
        dome["hv"].append(round(hv, 2))
        dome["ts"].append(round(T - 273.15, 3))
        dome["sl"].append(round(sl, 4))
        dome["sv"].append(round(sv, 4))

    # ---- h range from corners + dome ----
    cand = []
    for P in (Pmin, Pmax):
        for T in (Tmin, Tmax):
            try:
                cand.append(PropsSI("H", "P", P, "T", T, cp) / 1000.0)
            except Exception:
                pass
    cand += dome["hl"] + dome["hv"]
    h0, h1 = min(cand), max(cand)
    pad = 0.03 * (h1 - h0)
    h0, h1 = h0 - pad, h1 + pad

    # ---- crosshair grid: T, s, rho on (la, h) ----
    Ps = logspace(Pmin, Pmax, NP)
    Hs = linspace(h0, h1, NH)
    gT, gS, gR = [], [], []
    for P in Ps:
        for h in Hs:
            r = flash(CoolProp.HmassP_INPUTS, h * 1000.0, P)
            if r is None:
                gT.append(NAN); gS.append(NAN); gR.append(NAN)
            else:
                gT.append(r[0]); gS.append(r[1])
                gR.append(math.log10(r[2]) if r[2] > 0 else NAN)  # log10(rho)

    # ---- iso-curves ----
    def curve_pts(fn):
        """fn(P) -> h or None; sample over log P, None = break."""
        hs, las = [], []
        for P in logspace(Pmin, Pmax, NCRV):
            h = fn(P)
            hs.append(None if h is None else round(h, 2))
            las.append(round(la(P), 5))
        # trim leading/trailing breaks
        return dict(h=hs, la=las)

    curves = dict(T=[], s=[], r=[])

    # isotherms: exact dome split
    for TdegC in cfg["Tlevels"]:
        T = TdegC + 273.15
        if T < Tmin - 1e-9 or T > Tmax + 1e-9:
            continue
        hs, las = [], []
        try:
            Psat = PropsSI("P", "T", T, "Q", 0, cp) if T < Tc else None
        except Exception:
            Psat = None
        if Psat is not None and Pmin < Psat < Pmax:
            for P in logspace(Pmin, Psat, 60)[:-1]:      # vapor branch
                r = flash(CoolProp.PT_INPUTS, P, T)
                hs.append(None if r is None else round(r[3], 2))
                las.append(round(la(P), 5))
            try:                                          # dome jump
                hv = PropsSI("H", "T", T, "Q", 1, cp) / 1000.0
                hl = PropsSI("H", "T", T, "Q", 0, cp) / 1000.0
                hs += [round(hv, 2), round(hl, 2)]
                las += [round(la(Psat), 5)] * 2
            except Exception:
                pass
            liq_ok = T > Ttr                              # no liquid below triple
            if liq_ok:
                for P in logspace(Psat, Pmax, 60)[1:]:    # liquid branch
                    r = flash(CoolProp.PT_INPUTS, P, T)
                    hs.append(None if r is None else round(r[3], 2))
                    las.append(round(la(P), 5))
        else:
            for P in logspace(Pmin, Pmax, NCRV):
                r = flash(CoolProp.PT_INPUTS, P, T)
                hs.append(None if r is None else round(r[3], 2))
                las.append(round(la(P), 5))
        if any(v is not None for v in hs):
            curves["T"].append(dict(v=TdegC, h=hs, la=las))

    # isentropes (PS flash handles the dome continuously)
    sfin = [v for v in gS if v == v]
    slo, shi = sorted(sfin)[int(0.005 * len(sfin))], sorted(sfin)[int(0.995 * len(sfin))]
    step = cfg["sStep"]
    lev = math.ceil(slo / step) * step
    while lev <= shi:
        def f(P, s=lev):
            r = flash(CoolProp.PSmass_INPUTS, P, s * 1000.0)
            return None if r is None else r[3]
        c = curve_pts(f)
        if any(v is not None for v in c["h"]):
            curves["s"].append(dict(v=round(lev, 3), **c))
        lev += step

    # isochores (DmassP flash, two-phase included)
    rfin = [v for v in gR if v == v]
    rlo, rhi = 10 ** min(rfin), 10 ** max(rfin)
    for rho in RHO_CAND:
        if rho < rlo * 1.1 or rho > rhi * 0.95:
            continue
        def f(P, d=rho):
            r = flash(CoolProp.DmassP_INPUTS, d, P)
            return None if r is None else r[3]
        c = curve_pts(f)
        if any(v is not None for v in c["h"]):
            curves["r"].append(dict(v=rho, **c))

    hC = PropsSI("H", "T", Tc - 1e-3, "Q", 1, cp) / 1000.0

    return dict(
        label=cfg["label"],
        la0=round(la(Pmin), 5), la1=round(la(Pmax), 5), nP=NP,
        h0=round(h0, 2), h1=round(h1, 2), nH=NH,
        laC=round(la(Pc), 5), hC=round(hC, 2), TcC=round(Tc - 273.15, 2),
        laTr=round(la(Ptr), 5),
        gT=q16(gT), gS=q16(gS), gR=q16(gR),
        dome=dome, curves=curves,
    )


def main():
    out = {}
    for cfg in FLUIDS:
        print("generating", cfg["key"], "...")
        out[cfg["key"]] = gen_fluid(cfg)
    js = json.dumps(out, separators=(",", ":"), allow_nan=False)
    path = "PressureEnthalpy.html"
    src = open(path, encoding="utf-8").read()
    src, n = re.subn(r"/\*__D__\*/.*?/\*__E__\*/",
                     lambda m: "/*__D__*/" + js + "/*__E__*/", src, flags=re.S)
    assert n == 1, "marker not found"
    open(path, "w", encoding="utf-8").write(src)
    print("injected %.1f kB of data into %s" % (len(js) / 1024, path))


if __name__ == "__main__":
    main()
