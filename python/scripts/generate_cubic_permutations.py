"""
Vincent Ouazan-Reboul, 2025
Script to generate the orientation permutations when viewed under different bonds for the
`include/geometry/cubic.h` file.
The function you want to use is `output_formatted_bonds(output_file)`, which puts in a text file
the bond permutations as they will be formatted in the `cubic.h` header.
By default, this script outputs the results in `python/data/geometry/cubic.txt`. To change
output file, go to the last line of the script and give a file location as argument for
`output_formatted_bonds(output_file)`.
"""

import config as cfg
from geometry.cubic import CubicParticle


def format_bond_permutations():
    all_bond_permutations = CubicParticle().gen_bond_permutations()
    permutations_str = "{\n"
    for permutations in all_bond_permutations:
        permutations_str += "{" + str(permutations)[1:-1] + "},\n"
    permutations_str += "}"
    return permutations_str

def output_formatted_bonds(output_path = cfg.python_path/"data/geometry/cubic.txt"):
    with open(output_path, "w") as f:
        f.write(format_bond_permutations())
    return

output_formatted_bonds()
