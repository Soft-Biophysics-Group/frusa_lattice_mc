"""Executing frusa_mc runs, as opposed to generating their inputs."""

# `task` is a CLI module run as its own process; importing it here would make
# `python -m frusa_lattice_mc.runners.task` warn about a double import.
from . import local, stages

__all__ = ["local", "stages"]
