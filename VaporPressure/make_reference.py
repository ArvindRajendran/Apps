"""
make_reference.py - reference values for validating the JavaScript engine.

For every component and every correlation in vp_data.json, evaluate Psat with
the *Python* implementations in `chemicals` at 9 points spanning the published
validity range, and write reference.json.  check_engine.js then runs the same
cases through the JS engine shipped in the app and compares.

    python make_reference.py            # writes reference.json
"""

import json

from chemicals.dippr import EQ101
from chemicals.vapor_pressure import (Ambrose_Walton, Antoine, Lee_Kesler,
                                      TRC_Antoine_extended, Wagner,
                                      Wagner_original)

E = 2.718281828459045

data = json.load(open("vp_data.json"))
cases = []

for c in data["components"]:
    Tc, Pc, w = c["Tc"], c["Pc"], c["w"]
    methods = dict(c["m"])
    if "sub" in c:
        methods["sub"] = c["sub"]          # solid-vapour line, same ln convention
    for meth, p in methods.items():
        if meth in ("wagner36",):
            A, B, C, D, mTc, mPc, Tmin, Tmax = p
            f = lambda T: Wagner_original(T, mTc, mPc, A, B, C, D)
        elif meth in ("wagner25", "ppds"):
            A, B, C, D, mTc, mPc, Tmin, Tmax = p
            f = lambda T: Wagner(T, mTc, mPc, A, B, C, D)
        elif meth == "dippr":
            C1, C2, C3, C4, C5, Tmin, Tmax = p
            f = lambda T: EQ101(T, C1, C2, C3, C4, C5)
        elif meth == "antoine":
            A, B, C, Tmin, Tmax = p
            f = lambda T: Antoine(T, A, B, C)
        elif meth in ("landolt", "sub"):
            A, B, C, Tmin, Tmax = p
            f = lambda T: Antoine(T, A, B, C, base=E)
        elif meth == "antoineext":
            A, B, C, eTc, to, n, Ee, F, Tmin, Tmax = p
            f = lambda T: TRC_Antoine_extended(T, eTc, to, A, B, C, n, Ee, F)
        else:
            continue
        if Tmin is None or Tmax is None or Tmax <= Tmin:
            continue
        for i in range(9):
            T = Tmin + (Tmax - Tmin) * i / 8.0
            # A few published ranges run past the critical point, where there is
            # no saturation state.  The app refuses to evaluate there (several
            # correlations would silently clamp Tr to 1 and return Pc), so the
            # reference must not claim a value either.
            if meth != "sub" and Tc and T > Tc:
                continue
            try:
                P = f(T)
            except Exception:
                continue
            if P is None or P != P or P <= 0.0:
                continue
            cases.append([c["cas"], meth, round(T, 6), P])

    # corresponding-states methods, evaluated over the liquid range
    if Tc and Pc and w is not None:
        Tlo = max(c["Tm"] or 0.4 * Tc, 0.4 * Tc)
        Tlo = min(Tlo, 0.95 * Tc)          # never let the grid run past Tc
        for i in range(9):
            T = round(Tlo + (Tc - Tlo) * i / 8.0, 6)
            for name, fn in (("leekesler", Lee_Kesler),
                             ("ambrosewalton", Ambrose_Walton)):
                # Lee-Kesler overflows for strongly negative omega at low Tr
                try:
                    P = fn(T, Tc, Pc, w)
                except (OverflowError, ValueError, ZeroDivisionError):
                    continue
                if P == P and P not in (float("inf"), float("-inf")) and P > 0.0:
                    cases.append([c["cas"], name, T, P])

json.dump(cases, open("reference.json", "w"))
print(f"{len(cases)} reference cases written")
from collections import Counter
for k, v in Counter(x[1] for x in cases).most_common():
    print(f"  {k:16s} {v:6d}")
