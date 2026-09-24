# Headless re-execution of the ElectroCore routines

These scripts produced `results/sensitivity_analysis.csv` and the two correlation matrices in `results/`.
They import the ElectroCore Analyzer v2.1 source (not included in this repository; available from the
corresponding author on request) and run its own `fit_circuit`, `check_identifiability` and
`rank_models_correctly` routines without the graphical interface (`tkstub.py` replaces Tkinter).

For the sensitivity analysis, copies of `main_app.py` and `plotting_export.py` were made in which the
constants ρ_sev (0.98), the relative-error threshold (150 %), the window factor (5) and the Tier-1
promotion margin (10) are read from a dictionary instead of being fixed; nothing else was changed.
The automatic search was rerun (screening: 3 restarts, tol 1e-6, 300 evaluations; refinement:
12 restarts, tol 1e-8, 1000 evaluations; M1–M15) and the ranking recomputed for each setting.
The rerun reproduces the exported sessions (AIC within 0.05 units; identical penalties for the leading
candidates).
