"""Run the full analysis in order. Usage:  python run_all.py [--coalfield PATH_TO_KGS_REGIONS.shp]"""
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
extra = sys.argv[1:]
steps = ["01_regional_cv.py", "02_eky_transfer_test.py", "03_paired_sensor_comparison.py",
         "04_permutation_importance.py", "05_robustness_checks.py", "06_crosssensor_agreement.py",
         "07_figures.py"]
for step in steps:
    print(f"\n=== {step} ===", flush=True)
    args = [sys.executable, str(here / "scripts" / step)] + (extra if step.startswith("05") else [])
    subprocess.run(args, check=True, cwd=here)
