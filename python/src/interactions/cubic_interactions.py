"""
Vincent Ouazan-Reboul, 2025
Generating interaction maps for cubic particles.
"""

import config as cfg
from geometry.cubic import CubicParticle
from contact_utils import ContactMapWrapper

class CubicInteractions:
    """
    This class is essentially a wrapper for a ContactMapWrapper object.
    In ContactMapWrapper, contacts are encoded as (orientation, orientation, bond) sets.
    """
    def __init__(self, n_types):
        self.geometry = CubicParticle()
        self.cmap_wrapper = ContactMapWrapper(n_types)
        self.n_orientations = 24

    def get_contact_energy(face1, type1, face2, type2, energy):

    def assign_contact_energy(face1, type1, face2, type2, energy):
        """
        Sets contact between face1 of particles of type1 and face2 of particles of type2 to
        value energy.
        Sets all the (orientation, orientation, bond) sets verifying that contact to value
        energy in `cmap_wrapper` class member.
        """

        return
