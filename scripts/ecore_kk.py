"""Independent implementation of the linear Voigt Kramers-Kronig test, Eqs. (1)-(3) of the article.

weighting="none"   : unweighted bounded linear least squares, R_inf free, R_k >= 0 (the procedure
                     used by the ElectroCore software; reproduces its exported chi2_KK).
weighting="modulus": rows weighted by 1/|Z|, all coefficients >= 0 (the |Z|^-1-weighted check).
L="off" | "free" (sign-unconstrained, as in the software) | "positive" (L >= 0).
"""
import numpy as np
from scipy.optimize import lsq_linear

def kk_test(f, Z, M=20, weighting="none", L="off"):
    f = np.asarray(f, float); Z = np.asarray(Z, complex)
    w = 2 * np.pi * f
    tau = np.logspace(np.log10(1 / (2 * np.pi * f.max())), np.log10(1 / (2 * np.pi * f.min())), M)
    cols, lo, hi = [np.ones_like(w, complex)], [-np.inf], [np.inf]       # R_inf
    if L != "off":
        cols.append(1j * w)
        lo.append(-np.inf if L == "free" else 0.0); hi.append(np.inf)
    cols += [1 / (1 + 1j * w * t) for t in tau]
    lo += [0.0] * M; hi += [np.inf] * M
    A = np.array(cols).T
    if weighting == "modulus":
        lo[0] = 0.0
        s = 1 / np.abs(Z)
    else:
        s = np.ones_like(w)
    AA = np.vstack([(A.real.T * s).T, (A.imag.T * s).T]); bb = np.concatenate([Z.real * s, Z.imag * s])
    sc = np.abs(AA).max(0); sc[sc == 0] = 1
    # columns are scaled for conditioning; the bounds are scaled accordingly
    x = lsq_linear(AA / sc, bb, bounds=(np.array(lo) * sc, np.array(hi) * sc),
                   lsmr_tol='auto', max_iter=10000).x / sc
    Zk = A @ x
    d_re = (Z.real - Zk.real) / np.abs(Z); d_im = (Z.imag - Zk.imag) / np.abs(Z)
    out = dict(chi2_kk=float(np.mean(d_re**2 + d_im**2)), delta_re=d_re, delta_im=d_im, Z_kk=Zk, R_inf=float(x[0]))
    if L != "off":
        out["L"] = float(x[1])
    return out

def window(f, Z, fmin, fmax):
    m = (f >= fmin * 0.999) & (f <= fmax * 1.001)
    return f[m], Z[m], m
