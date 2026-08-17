"""
generate_data.py - build the embedded component database for VaporPressure.html

Source: the `chemicals` package (MIT, Caleb Bell), which redistributes the
published coefficient tables listed below.  Run:

    pip install chemicals
    python generate_data.py > vp_data.json

The curated component set is the union of the six "engineering" coefficient
tables (Poling/Reid, Perry's, VDI, McGarry).  A compound is in the set if at
least one of those tables carries it -- i.e. it was considered common enough
to appear in a standard handbook.  The much larger Landolt-Boernstein Antoine
table (6346 compounds) is used only to add a method to components already in
the set, never to add new components.

CONVENTIONS (the main trap in this data):
  antoine     log10(P/Pa) = A - B/(T+C)          T in K
  landolt     ln(P/Pa)    = A - B/(T+C)          T in K   <-- natural log!
  antoineext  TRC extended form, log10(P/Pa)
  wagner36    Wagner 3,6 form (McGarry 1983)
  wagner25    Wagner 2.5,5 form (Poling)
  ppds        Wagner 2.5,5 form with VDI PPDS coefficients
  dippr       DIPPR-101: ln(P/Pa) = C1 + C2/T + C3 ln T + C4 T^C5
Feeding Landolt coefficients through the log10 form gives answers wrong by
seven orders of magnitude, so every method carries its own evaluator.
"""

import json
import math
import re
import sys

from chemicals import vapor_pressure as vp
from chemicals.acentric import omega as omega_f
from chemicals.critical import Pc as Pc_f
from chemicals.critical import Tc as Tc_f
from chemicals.identifiers import search_chemical
from chemicals.phase_change import Tb as Tb_f
from chemicals.phase_change import Tm as Tm_f

vp.load_vapor_pressure_dfs()

# Tables that define membership in the curated set, in descending order of
# the accuracy generally attributed to them.
PREMIUM = [
    "Psat_data_WagnerMcGarry",
    "Psat_data_WagnerPoling",
    "Psat_data_VDI_PPDS_3",
    "Psat_data_Perrys2_8",
    "Psat_data_AntoineExtended",
    "Psat_data_AntoinePoling",
]


def num(x, nd=8):
    """Round a value for compact JSON; NaN/inf -> None."""
    if x is None:
        return None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    if math.isnan(x) or math.isinf(x):
        return None
    if x == 0.0:
        return 0.0
    # keep nd significant figures, not decimal places
    mag = math.floor(math.log10(abs(x)))
    r = round(x, max(0, nd - 1 - mag))
    return r


def disp_hint(hit, cas):
    """Best available name at the point the sanity checks run."""
    return hit.common_name if hit and hit.common_name else cas


def safe(fn, cas):
    try:
        v = fn(cas)
        return None if v is None else float(v)
    except Exception:
        return None


# Synonyms worth keeping for search: structural prefixes a chemical engineer
# actually types (n-hexane, isobutane, o-xylene) and refrigerant designations
# (R134a, HFC-134a).  Everything else in the synonym lists is brand names and
# transliterations, which only make the picker noisier.
ALIAS_OK = re.compile(
    r"^(n-|iso|sec-|tert-|neo|cyclo|o-|m-|p-|ortho-|meta-|para-)"
    r"|^(r|hfc|hcfc|cfc)[- ]?\d"
    r"|^refrigerant\s+r",
    re.I,
)


# Nomenclature prefixes that stay lower-case however the name is displayed.
PREFIX = re.compile(r"^(n|o|m|p|d|l|r|s|z|e|sec|tert|cis|trans|alpha|beta|"
                    r"gamma|ortho|meta|para|sym|asym)-", re.I)


def aliases(hit, table_names):
    """Search aliases: the handbook names plus useful structural synonyms."""
    out, seen = [], {hit.common_name.lower()}
    for nm in table_names:
        k = nm.lower()
        if k and k not in seen:
            seen.add(k)
            out.append(nm)
    n = 0
    for s in hit.synonyms:
        k = s.lower()
        if k in seen or len(s) > 32 or not ALIAS_OK.match(s):
            continue
        seen.add(k)
        out.append(s)
        n += 1
        if n >= 5:
            break
    return out


def main():
    members = set()
    for s in PREMIUM:
        members |= set(getattr(vp, s).index)

    landolt = vp.Psat_data_Landolt_Antoine
    sublim = vp.Psat_data_Landolt_Antoine  # placeholder, sublimation added below
    psub = vp.Psub_data_Landolt_Antoine

    comps = []
    dropped = []
    bad_critical = []
    bad_melting = []
    for cas in sorted(members):
        try:
            hit = search_chemical(cas)
        except Exception:
            hit = None
        if hit is None or not hit.common_name:
            dropped.append(cas)
            continue

        Tc = safe(Tc_f, cas)
        Pc = safe(Pc_f, cas)
        w = safe(omega_f, cas)
        Tb = safe(Tb_f, cas)
        Tm = safe(Tm_f, cas)

        # Sanity-check the critical data against the phase-change data.  Two
        # records in this set are internally impossible: phenanthrene is listed
        # with Tc = 0.869 K (a factor-1000 slip -- its own Tb is 611 K), and
        # vinyl acetylene with Tm = 476 K above its Tc of 455 K.  A boiling
        # point at or above the critical temperature condemns Tc itself, so
        # drop the critical set; otherwise the melting point is the outlier.
        if Tc and Tb and Tb >= Tc:
            bad_critical.append((cas, disp_hint(hit, cas), Tc, Tb))
            Tc = Pc = w = None
        elif Tc and Tm and Tm >= Tc:
            bad_melting.append((cas, disp_hint(hit, cas), Tm, Tc))
            Tm = None

        m = {}
        tnames = []          # names this compound carries in the handbook tables
        ppds_name = None     # VDI's names are the best-formatted of the lot

        # --- Wagner 3,6 (McGarry 1983) -----------------------------------
        if cas in vp.Psat_data_WagnerMcGarry.index:
            r = vp.Psat_data_WagnerMcGarry.loc[cas]
            tnames.append(str(r.Name).strip())
            m["wagner36"] = [num(r.A), num(r.B), num(r.C), num(r.D),
                             num(r.Tc), num(r.Pc), num(r.Tmin), num(r.Tc)]

        # --- Wagner 2.5,5 (Poling) ---------------------------------------
        if cas in vp.Psat_data_WagnerPoling.index:
            r = vp.Psat_data_WagnerPoling.loc[cas]
            tnames.append(str(r.Name).strip())
            m["wagner25"] = [num(r.A), num(r.B), num(r.C), num(r.D),
                             num(r.Tc), num(r.Pc), num(r.Tmin), num(r.Tmax)]

        # --- VDI PPDS (also Wagner 2.5,5) --------------------------------
        if cas in vp.Psat_data_VDI_PPDS_3.index:
            r = vp.Psat_data_VDI_PPDS_3.loc[cas]
            ppds_name = str(r.Chemical).strip()
            tnames.append(ppds_name)
            # no explicit range published; melting point to critical point
            m["ppds"] = [num(r.A), num(r.B), num(r.C), num(r.D),
                         num(r.Tc), num(r.Pc), num(r.Tm), num(r.Tc)]

        # --- DIPPR-101 (Perry's 8th ed. Table 2-8) -----------------------
        if cas in vp.Psat_data_Perrys2_8.index:
            r = vp.Psat_data_Perrys2_8.loc[cas]
            tnames.append(str(r.Chemical).strip())
            m["dippr"] = [num(r.C1), num(r.C2), num(r.C3), num(r.C4),
                          num(r.C5), num(r.Tmin), num(r.Tmax)]

        # --- TRC extended Antoine ----------------------------------------
        if cas in vp.Psat_data_AntoineExtended.index:
            r = vp.Psat_data_AntoineExtended.loc[cas]
            tnames.append(str(r.Chemical).strip())
            m["antoineext"] = [num(r.A), num(r.B), num(r.C), num(r.Tc),
                               num(r.to), num(r.n), num(r.E), num(r.F),
                               num(r.Tmin), num(r.Tmax)]

        # --- Antoine, log10, Pa (Poling/Reid) ----------------------------
        if cas in vp.Psat_data_AntoinePoling.index:
            r = vp.Psat_data_AntoinePoling.loc[cas]
            tnames.append(str(r.Chemical).strip())
            m["antoine"] = [num(r.A), num(r.B), num(r.C),
                            num(r.Tmin), num(r.Tmax)]

        # --- Antoine, ln, Pa (Landolt-Boernstein) ------------------------
        if cas in landolt.index:
            r = landolt.loc[cas]
            m["landolt"] = [num(r.A), num(r.B), num(r.C),
                            num(r.Tmin), num(r.Tmax)]

        # --- sublimation (Landolt), optional solid branch ----------------
        sub = None
        if cas in psub.index:
            r = psub.loc[cas]
            sub = [num(r.A), num(r.B), num(r.C), num(r.Tmin), num(r.Tmax)]

        if not m:
            dropped.append(cas)
            continue

        # Display name: VDI's names are cleanest.  Failing that, if the
        # database's common name never appears in any handbook table it is
        # usually a pharmaceutical or trade name (CAS 811-97-2 is "norflurane",
        # not something anyone would look up), so take the most descriptive
        # handbook name instead.
        lower_tnames = {t.lower() for t in tnames}
        if ppds_name:
            disp = ppds_name
        elif tnames and hit.common_name.lower() not in lower_tnames:
            disp = max(tnames, key=len)
        else:
            disp = hit.common_name
        # Capitalize for display, but never touch a nomenclature prefix:
        # "n-hexane" must not become "N-Hexane", and cis-/tert-/o- likewise.
        if disp.islower() and not PREFIX.match(disp):
            disp = disp[:1].upper() + disp[1:]
        al = [a for a in aliases(hit, tnames) if a.lower() != disp.lower()]
        if hit.common_name.lower() != disp.lower():
            al.insert(0, hit.common_name)

        c = {
            "cas": cas,
            "n": disp,
            "f": hit.formula,
            "mw": num(hit.MW, 6),
            "Tc": num(Tc), "Pc": num(Pc), "w": num(w),
            "Tb": num(Tb), "Tm": num(Tm),
            "m": m,
        }
        if al:
            c["a"] = al
        if sub:
            c["sub"] = sub
        comps.append(c)

    comps.sort(key=lambda c: c["n"].lower())

    out = {
        "meta": {
            "generated_by": "generate_data.py via the chemicals package (MIT)",
            "n_components": len(comps),
            "methods": {
                "wagner36": "Wagner (3,6) - McGarry, Ind. Eng. Chem. Process Des. Dev. 22 (1983) 313",
                "wagner25": "Wagner (2.5,5) - Poling, Prausnitz & O'Connell, Properties of Gases and Liquids, 5th ed., App. A",
                "ppds": "VDI PPDS (Wagner 2.5,5 form) - VDI Heat Atlas, 2nd ed., 2010",
                "dippr": "DIPPR-101 - Perry's Chemical Engineers' Handbook, 8th ed., Table 2-8",
                "antoineext": "TRC extended Antoine - Poling et al., 5th ed., App. A",
                "antoine": "Antoine (log10, Pa) - Poling et al., 5th ed., App. A",
                "landolt": "Antoine (ln, Pa) - Landolt-Boernstein, New Series IV/20",
                "sub": "Sublimation Antoine (ln, Pa) - Landolt-Boernstein IV/20",
            },
        },
        "components": comps,
    }
    json.dump(out, sys.stdout, separators=(",", ":"))
    print(file=sys.stderr)
    print(f"components: {len(comps)}   dropped (no name or no coefficients): {len(dropped)}",
          file=sys.stderr)
    from collections import Counter
    cnt = Counter()
    for c in comps:
        for k in c["m"]:
            cnt[k] += 1
        if "sub" in c:
            cnt["sub"] += 1
    for k, v in cnt.most_common():
        print(f"  {k:12s} {v:5d}", file=sys.stderr)
    for cas, nm, a, b in bad_critical:
        print(f"  ! {nm} ({cas}): Tb={b} >= Tc={a} — critical data dropped", file=sys.stderr)
    for cas, nm, a, b in bad_melting:
        print(f"  ! {nm} ({cas}): Tm={a} >= Tc={b} — melting point dropped", file=sys.stderr)
    nmulti = sum(1 for c in comps if len(c["m"]) >= 2)
    print(f"  components with >=2 methods: {nmulti}", file=sys.stderr)
    ncs = sum(1 for c in comps if c["Tc"] and c["Pc"] and c["w"] is not None)
    print(f"  components usable by Lee-Kesler / Ambrose-Walton: {ncs}", file=sys.stderr)


if __name__ == "__main__":
    main()
