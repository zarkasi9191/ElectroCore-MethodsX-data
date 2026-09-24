"""Recompute every number of the article that follows from the data, and compare with the printed value.

Run:  python scripts/reproduce_numbers.py      (writes results/reproduced_numbers.csv)
Numbers that require the ElectroCore software itself (the automatic ranking, the identifiability
penalties and the sensitivity analysis) are read from the exports / results files instead.
"""
import csv, math
import numpy as np
from scipy import stats
from ecore_io import ROOT, load_fit_export, load_measurement, sheet_dict, table
from ecore_kk import kk_test, window

checks = []
def check(item, article, value, rel_tol=0.02, where=""):
    ok = math.isclose(value, article, rel_tol=rel_tol, abs_tol=1e-12)
    checks.append([item, where, article, value, "PASS" if ok else "FAIL"])

# ---------- Kramers-Kronig validation (Eqs. 1-3), software-equivalent unweighted fit
fA, ZA, ZfA = load_fit_export("A"); fB, ZB, ZfB = load_fit_export("B")
check("chi2_KK dataset A, M=20", 2.7e-5, kk_test(fA, ZA)["chi2_kk"], 0.03, "Sec. 4.1, 10; Fig. 6c")
check("chi2_KK dataset B, M=20", 7.29e-5, kk_test(fB, ZB)["chi2_kk"], 0.01, "MV 1, 7; Fig. 6f")
for M, a, b in [(10, 2.5e-5, 1.6e-4), (30, 2.1e-5, 6.4e-5)]:
    check(f"chi2_KK dataset A, M={M}", a, kk_test(fA, ZA, M)["chi2_kk"], 0.05, "Sec. 4.2")
    check(f"chi2_KK dataset B, M={M}", b, kk_test(fB, ZB, M)["chi2_kk"], 0.05, "Sec. 4.2")
check("chi2_KK dataset B without two highest-frequency points", 4.4e-5, kk_test(fB[:-2], ZB[:-2])["chi2_kk"], 0.03, "MV 1")
rL = kk_test(fA, ZA, L="free")
check("dataset A with sign-free L: chi2_KK", 8.4e-6, rL["chi2_kk"], 0.03, "Sec. 4.1")
check("dataset A with sign-free L: L (H)", -57e-6, rL["L"], 0.03, "Sec. 4.1")
check("dataset B, |Z|^-1-weighted, L>=0", 7.0e-5, kk_test(fB, ZB, weighting="modulus", L="positive")["chi2_kk"], 0.03, "MV 7")

# ferri/ferrocyanide spectra (MV 6, Fig. 8)
art = {10: (4.1e-2, 7.6e-5, None), 15: (5.0e-2, 4.4e-4, 1.7e-5), 20: (2.3e-1, 1.4e-4, 6.1e-5)}
for c, (full, le10k, trimmed) in art.items():
    f, Z = load_measurement(c)
    check(f"{c} mM, 10 Hz-100 kHz", full, kk_test(f, Z)["chi2_kk"], 0.03, "MV 6")
    check(f"{c} mM, 10 Hz-10 kHz", le10k, kk_test(*window(f, Z, 10, 1e4)[:2])["chi2_kk"], 0.03, "MV 6")
    if trimmed:
        check(f"{c} mM, 30 Hz-10 kHz", trimmed, kk_test(*window(f, Z, 30, 1e4)[:2])["chi2_kk"], 0.03, "MV 6")
wfull = [kk_test(*load_measurement(c), weighting="modulus")["chi2_kk"] for c in (10, 15, 20)]
wval = [kk_test(*window(*load_measurement(10), 10, 1e4)[:2], weighting="modulus")["chi2_kk"]] + \
       [kk_test(*window(*load_measurement(c), 30, 1e4)[:2], weighting="modulus")["chi2_kk"] for c in (15, 20)]
check("weighted check, full window, lowest", 2.7e-2, min(wfull), 0.05, "MV 6")
check("weighted check, full window, highest", 1.8e-1, max(wfull), 0.05, "MV 6")
check("weighted check, validated windows, lowest", 1.6e-5, min(wval), 0.05, "MV 6")
check("weighted check, validated windows, highest", 5.2e-5, max(wval), 0.05, "MV 6")

# ---------- information criteria from the exported weighted chi-square (Eqs. 6-9, 14)
mA = sheet_dict("A", "Advanced_Metrics"); N = int(mA["N_Points"]); nobs = 2 * N
head, rows = table("A", "Auto_Fit_Ranking"); rk = {r[1].split(".")[0]: dict(zip(head, r)) for r in rows}
chi2w = {"8": mA["Chi_Squared_Weighted"], "7": rk["7"]["Chi2_red"] * (nobs - 3 - 1)}
aic = lambda c, P: nobs * math.log(c / nobs) + 2 * P
bic = lambda c, P: nobs * math.log(c / nobs) + P * math.log(nobs)
check("AIC M7 (dataset A)", -959.998, aic(chi2w["7"], 3), 1e-5, "Table 8; MV 7")
check("AIC M8 (dataset A)", -976.146, aic(chi2w["8"], 4), 1e-5, "Table 8; MV 7")
check("AICc M8 (dataset A)", -975.697, aic(chi2w["8"], 4) + 2 * 4 * 5 / (nobs - 4 - 1), 1e-5, "Eq. 8")
check("BIC M8 (dataset A)", -965.973, bic(chi2w["8"], 4), 1e-5, "Table 8")
check("Delta AIC M8-M7, Eq. (14)", -16.15, aic(chi2w["8"], 4) - aic(chi2w["7"], 3), 1e-3, "Eq. 14")
check("Delta BIC M8-M7", -13.61, bic(chi2w["8"], 4) - bic(chi2w["7"], 3), 1e-3, "Sec. 10.2")
pA = {r[0]: r for r in table("A", "Fitting_Parameters")[1]}
n, se = pA["CPE_n"][1], pA["CPE_n"][3]; nu = nobs - 4 - 1
t = (1 - n) / se
check("t statistic on n", 4.47, t, 2e-3, "Sec. 10.2")
check("t-test p (t distribution, nu=89)", 2.3e-5, 2 * stats.t.sf(t, nu), 0.03, "Sec. 10.2")
Fv = (chi2w["7"] - chi2w["8"]) / (chi2w["8"] / nu)
check("F statistic (1, 89)", 19.0, Fv, 5e-3, "Sec. 10.2")
check("F-test p", 3.6e-5, stats.f.sf(Fv, 1, nu), 0.03, "Sec. 10.2")
Rs, Rct, Q = pA["Rs"][1], pA["Rct"][1], pA["CPE_Q"][1]
Rp = Rs * Rct / (Rs + Rct); Ceff = Q ** (1 / n) * Rp ** ((1 - n) / n)
check("C_eff (Brug), Eq. (13), pF", 31715, Ceff * 1e12, 1e-3, "Table 7; MV 7")
check("R_p (ohm)", 494.7, Rp, 1e-3, "Sec. 10.1")
check("R_s deviation from 1 kOhm (%)", -2.3, 100 * (Rs / 1000 - 1), 0.02, "Sec. 10.1")
check("R_1 deviation from 1 kOhm (%)", 0.2, 100 * (Rct / 1000 - 1), 0.1, "Sec. 10.1")

# ---------- dataset B parameter recovery (Table 9)
pB = {r[0]: r for r in table("B", "Fitting_Parameters")[1]}
mB = sheet_dict("B", "Advanced_Metrics")
check("AIC M13 (dataset B)", -1021.41, mB["AIC"], 1e-5, "Table 10")
# export names: Rct/Cdl form the 470 nF arc (R1, C1 of Table 9), R2/C2 the 33 nF arc; capacitances exported in pF
for lab, key, nom, scale in [("R_s", "Rs", 1000, 1), ("R_1", "Rct", 1000, 1), ("C_1 (nF)", "Cdl", 470, 1e-3),
                             ("R_2", "R2", 1000, 1), ("C_2 (nF)", "C2", 33, 1e-3)]:
    v = pB[key][1] * scale
    art_dev = {"R_s": -0.17, "R_1": -1.72, "C_1 (nF)": 2.00, "R_2": -2.36, "C_2 (nF)": 1.70}[lab]
    check(f"dataset B {lab}: deviation from nominal (%)", art_dev, 100 * (v / nom - 1), 0.02, "Table 9")

# ---------- residual diagnostics (MV 3), recomputed from the exported fit
dw = lambda e: float(np.sum(np.diff(e) ** 2) / np.sum(e ** 2))
DW, KU, QQ, SWP = [], [], [], []
for f, Z, Zf in (load_fit_export("A"), load_fit_export("B")):
    raw = Z - Zf                                   # Durbin-Watson: raw residual components, ordered by f
    DW += [dw(raw.real), dw(raw.imag)]
    for e in (raw.real / np.abs(Z), raw.imag / np.abs(Z)):   # shape statistics: normalized residuals
        KU.append(float(stats.kurtosis(e, bias=False))); SWP.append(float(stats.shapiro(e).pvalue))
        osm, osr = stats.probplot(e, dist="norm")[0]; QQ.append(float(np.corrcoef(osm, osr)[0, 1] ** 2))
check("Durbin-Watson, lowest component", 0.09, min(DW), 0.06, "MV 3")
check("Durbin-Watson, highest component", 0.17, max(DW), 0.03, "MV 3")
check("excess kurtosis, lowest", 3.4, min(KU), 0.03, "MV 3")
check("excess kurtosis, highest", 11.7, max(KU), 0.03, "MV 3")
check("Q-Q R^2, lowest", 0.55, min(QQ), 0.03, "MV 3")
checks.append(["Shapiro-Wilk p, largest (article: p <= 1e-4)", "MV 3", 1e-4, max(SWP), "PASS" if max(SWP) <= 1.05e-4 else "FAIL"])

out = ROOT / "results" / "reproduced_numbers.csv"
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh); w.writerow(["quantity", "article location", "article value", "recomputed", "status"]); w.writerows(checks)
nfail = sum(1 for c in checks if c[4] == "FAIL")
for c in checks: print(f"{c[4]:>9}  {c[0]:<55} article {c[2]:<12g} recomputed {c[3]:.6g}")
print(f"\n{len(checks)} quantities, {nfail} failed -> {out}")
