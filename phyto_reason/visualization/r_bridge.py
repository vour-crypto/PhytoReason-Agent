"""Optional R plotting bridge.

R is intentionally not required. Deployments can add a project-specific
adapter later; Python plotting remains the deterministic fallback.
"""

from __future__ import annotations


def r_available() -> bool:
    return False


def _unavailable(*args, **kwargs):
    return None


r_heatmap = r_volcano = r_quadrant = _unavailable
