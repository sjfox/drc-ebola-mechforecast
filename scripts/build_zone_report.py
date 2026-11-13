"""
Build reports/<name>.html from reports/<name>.template.html, inlining the figures it
references as {{img:path}}, so the page is a single self-contained file.

  uv run python scripts/build_zone_report.py                  # zone_dynamics (after zone_series.R)
  uv run python scripts/build_zone_report.py forecast_review  # after plot_forecasts_final.R
"""

import base64
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"


def inline(match: re.Match) -> str:
    path = ROOT / match.group(1)
    data = base64.b64encode(path.read_bytes()).decode()
    return f'<img src="data:image/png;base64,{data}" alt="{path.stem}" loading="lazy">'


name = sys.argv[1] if len(sys.argv) > 1 else "zone_dynamics"
template = (REPORTS / f"{name}.template.html").read_text()
out = REPORTS / f"{name}.html"
out.write_text(re.sub(r"\{\{img:([^}]+)\}\}", inline, template))
print(f"Wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
