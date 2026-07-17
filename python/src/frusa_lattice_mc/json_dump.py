# Copyright (c) 2024 Soft Biophysics Group LPTMS
# Part of frusa_mc, released under BSD 3-Clause License.

import json
from pathlib import Path


def make_json_file(data, address) -> None:
    """Write `data` (a dict) to `address` as JSON, creating parent dirs."""
    address = Path(address)
    address.parent.mkdir(parents=True, exist_ok=True)
    with address.open("w") as write_file:
        json.dump(data, write_file, indent=4)


def load_json(path) -> dict:
    """Read a JSON file into a dict."""
    with Path(path).open() as f:
        return json.load(f)
