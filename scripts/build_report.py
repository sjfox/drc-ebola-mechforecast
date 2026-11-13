"""
Inline the verification results into the report page.

Run after scripts/verify_example.py:  uv run python scripts/build_report.py
"""

from pathlib import Path

REPORTS = Path(__file__).resolve().parents[1] / "reports"

template = (REPORTS / "model_demo.template.html").read_text()
data = (REPORTS / "verify_example.json").read_text().replace("</", "<\\/")
out = REPORTS / "model_demo.html"
out.write_text(template.replace("/*__DATA__*/", data))
print(f"Wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
