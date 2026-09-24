"""Regenerate the processed CSV files, the reproduced-number report and the data figures (Figs. 6-8)."""
import runpy
from pathlib import Path
HERE = Path(__file__).resolve().parent
for script in ["export_csv.py", "reproduce_numbers.py", "fig6_passive_networks.py",
               "fig7_identifiability.py", "fig8_ferri_ferrocyanide.py"]:
    print(f"\n=== {script}")
    runpy.run_path(str(HERE / script), run_name="__main__")
