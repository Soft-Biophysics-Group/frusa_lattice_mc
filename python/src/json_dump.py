# Copyright (c) 2024 Soft Biophysics Group LPTMS
# Part of frusa_mc, released under BSD 3-Clause License.

import json

def make_json_file(data,address):
    """
    Writes the data, provided as a python dictionary, to a .json file
    """
    with open(address,"w") as write_file:
        json.dump(data, write_file, indent=2)
