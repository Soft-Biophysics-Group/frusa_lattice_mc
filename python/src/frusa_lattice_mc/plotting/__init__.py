"""Plotting interfaces with optional 3D dependencies loaded on demand."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .plot_3d_utils import BlenderPlot

__all__ = ["BlenderPlot"]


def __getattr__(name: str):
    if name == "BlenderPlot":
        from .plot_3d_utils import BlenderPlot

        return BlenderPlot
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
