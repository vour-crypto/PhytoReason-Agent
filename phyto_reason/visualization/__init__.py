"""Publication-oriented plotting helpers used by the workflow reports.

The plotting layer is deliberately optional: analysis remains usable when
matplotlib is unavailable, while installed deployments receive PNG/SVG assets
under ``api/static/figures``.
"""

from .figure_exporter import FigureExporter

__all__ = ["FigureExporter"]
