"""`python -m scripts` from a checkout, or `python -m informational_flux` when installed."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cli import main  # noqa: E402

raise SystemExit(main())
