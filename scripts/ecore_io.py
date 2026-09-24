"""Readers for the ElectroCore Analyzer Excel exports used in the article."""
from pathlib import Path
import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
FILES = {
    "A": RAW / "dataset_A" / "Hasil_Fitting_1_Puncak.xlsx",
    "B": RAW / "dataset_B" / "Hasil_Fitting_2_Puncak.xlsx",
    "A_ranking": RAW / "dataset_A" / "model_comparison_1_puncak.xlsx",
    "B_ranking": RAW / "dataset_B" / "model_comparison_2_puncak.xlsx",
    10: RAW / "ferri_ferrocyanide" / "PLN_10_mM.xlsx",
    15: RAW / "ferri_ferrocyanide" / "PLN_15_mM.xlsx",
    20: RAW / "ferri_ferrocyanide" / "PLN_20_mM.xlsx",
}

def _rows(path, sheet):
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    rows = list(ws.iter_rows(values_only=True))
    return rows[0], rows[1:]

def load_fit_export(key):
    """Fitted-session export (datasets A, B): returns f, Z_exp, Z_fit sorted by frequency.
    Z = Z' + jZ'' with Z'' < 0 for a capacitive point."""
    _, rows = _rows(FILES[key], "Fitted_Data")
    a = np.array([r[:5] for r in rows if r[0] is not None], float)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1] + 1j * a[:, 2], a[:, 3] + 1j * a[:, 4]

def load_measurement(conc_mM):
    """Raw measurement export (ferri/ferrocyanide). The file stores Z_Imag as -Z''."""
    _, rows = _rows(FILES[conc_mM], "EIS_Data")
    a = np.array([r[:3] for r in rows if r[0] is not None], float)
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1] - 1j * a[:, 2]

def sheet_dict(key, sheet):
    """Two-column sheets (Advanced_Metrics, Model_Info, ...) as a dict."""
    _, rows = _rows(FILES[key], sheet)
    return {r[0]: r[1] for r in rows if r and r[0] is not None}

def table(key, sheet):
    head, rows = _rows(FILES[key], sheet)
    return list(head), [list(r) for r in rows if r and r[0] is not None]
