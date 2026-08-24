"""Plotting helpers.

`BlenderPlot` pulls in `bpy`, which only exists in the dedicated Python 3.11
render environment (`.venv-blender`) - bpy publishes no cp312 wheels. Importing
it lazily keeps the matplotlib-only modules in this package (chain, triangular,
square, fcc, plotting_utils) usable from the main 3.12 environment.
"""

__all__ = ["BlenderPlot"]


def __getattr__(name):
    if name == "BlenderPlot":
        from .plot_3d_utils import BlenderPlot

        return BlenderPlot
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + __all__)
