from __future__ import annotations

from pathlib import Path


class FigureExporter:
    """Save matplotlib figures using one consistent, web-accessible layout."""

    @staticmethod
    def save(fig, path: str, dpi: int = 150) -> str:
        """Save the figure; PNG targets also get a same-basename SVG twin
        for SCI vector export (best-effort, never fails the PNG output)."""
        target = Path(path)
        if target.suffix.lower() not in {".png", ".svg", ".pdf"}:
            target = target.with_suffix(".png")
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=dpi, bbox_inches="tight", facecolor="white")
        if target.suffix.lower() == ".png":
            try:
                fig.savefig(target.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
            except Exception:
                pass
        return str(target)
