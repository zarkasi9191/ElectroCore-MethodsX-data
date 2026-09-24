"""Convert the raw ElectroCore exports into plain CSV files (data/processed)."""
import csv, json
import numpy as np
from ecore_io import ROOT, load_fit_export, load_measurement, table, sheet_dict

OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

def write(name, header, rows):
    with open(OUT / name, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)

for key, label in [("A", "datasetA"), ("B", "datasetB")]:
    f, Z, Zf = load_fit_export(key)
    write(f"{label}_spectrum_and_fit.csv",
          ["frequency_Hz", "Z_real_ohm", "Z_imag_ohm", "Zfit_real_ohm", "Zfit_imag_ohm"],
          [[a, b.real, b.imag, c.real, c.imag] for a, b, c in zip(f, Z, Zf)])
    head, rows = table(key, "Fitting_Parameters"); write(f"{label}_selected_fit_parameters.csv", head, rows)
    head, rows = table(key, "Auto_Fit_Ranking"); write(f"{label}_automatic_ranking.csv", head, rows)
    met = sheet_dict(key, "Advanced_Metrics"); write(f"{label}_selected_fit_metrics.csv", ["metric", "value"], met.items())

for c in (10, 15, 20):
    f, Z = load_measurement(c)
    write(f"ferri_ferrocyanide_{c}mM_spectrum.csv", ["frequency_Hz", "Z_real_ohm", "Z_imag_ohm"],
          [[a, b.real, b.imag] for a, b in zip(f, Z)])
print("CSV files written to", OUT)

# measurement settings of the ferri/ferrocyanide spectra
rows = []
for c in (10, 15, 20):
    p = sheet_dict(c, "Parameters")
    rows.append([c] + [p.get(k) for k in ("Timestamp", "Start Frequency (Hz)", "Stop Frequency (Hz)", "Amplitude (mV)",
                                            "DC Bias (mV)", "Number of Points", "TIA Resistor", "PGA Gain", "Duration (s)", "Calibration Factor")])
write("ferri_ferrocyanide_measurement_settings.csv",
      ["concentration_mM", "timestamp", "f_start_Hz", "f_stop_Hz", "amplitude_mV", "dc_bias_mV", "points", "R_TIA", "PGA", "duration_s", "calibration_factor"], rows)
