from openff.toolkit.typing.engines.smirnoff import ImproperTorsionHandler, ParameterAttribute, IndexedParameterAttribute
from openff.toolkit.typing.engines.smirnoff.parameters import _allow_only
from simtk import openmm

# TODO: There's a lot of duplicated code in ProperTorsionHandler and ImproperTorsionHandler
class AmberImproperTorsionHandler(ImproperTorsionHandler):
    """Handle SMIRNOFF ``<AmberImproperTorsions>`` tags

    .. warning :: This API is experimental and subject to change.
    """

    # class AmberImproperTorsionType(ParameterType):
    #     """A SMIRNOFF torsion type for improper torsions.

    #     .. warning :: This API is experimental and subject to change.
    #     """
    #     _VALENCE_TYPE = 'AmberImproperTorsion'
    #     _ELEMENT_NAME = 'AmberImproper'

    #     periodicity = IndexedParameterAttribute(converter=int)
    #     phase = IndexedParameterAttribute(unit=unit.degree)
    #     k = IndexedParameterAttribute(unit=unit.kilocalorie_per_mole)
    #     idivf = IndexedParameterAttribute(default=None, converter=float)


    _TAGNAME = 'AmberImproperTorsions'  # SMIRNOFF tag name to process
    _INFOTYPE = ImproperTorsionHandler.ImproperTorsionType  # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce  # OpenMM force class to create

    potential = ParameterAttribute(
        default='k*(1+cos(periodicity*theta-phase))',
        converter=_allow_only(['k*(1+cos(periodicity*theta-phase))'])
    )
    default_idivf = ParameterAttribute(default='auto')

    from openff.toolkit.topology.topology import _TransformedDict
    class AmberImproperDict(_TransformedDict):
        """Symmetrize improper torsions."""
        
        def __keytransform__(self, key):
            """Only keep one improper per center, and assume that the central atom is the third one."""
            # Ensure key is a tuple
            key = tuple(key)
            # Retrieve connected atoms
            connectedatoms = [key[0], key[1], key[3]]
            # Sort connected atoms
            connectedatoms.sort()
            # Re-store connected atoms
            key = tuple(
                [connectedatoms[0], connectedatoms[1], key[2], connectedatoms[2]])
            return (key)

    def find_matches(self, entity):
        """Find the improper torsions in the topology/molecule matched by a parameter type.

        Parameters
        ----------
        entity : openforcefield.topology.Topology
            Topology to search.

        Returns
        ---------
        matches : ImproperDict[Tuple[int], ParameterHandler._Match]
            ``matches[atom_indices]`` is the ``ParameterType`` object
            matching the 4-tuple of atom indices in ``entity``.

        """
        return self._find_matches(entity, transformed_dict_cls=self.AmberImproperDict)

    def create_force(self, system, topology, **kwargs):
        #force = super(ImproperTorsionHandler, self).create_force(system, topology, **kwargs)
        #force = super().create_force(system, topology, **kwargs)
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [
            f for f in existing if type(f) == openmm.PeriodicTorsionForce
        ]
        if len(existing) == 0:
            force = openmm.PeriodicTorsionForce()
            system.addForce(force)
        else:
            force = existing[0]

        # Add all improper torsions to the system
        improper_matches = self.find_matches(topology)
        # The atom indices in the key of the dictionary match have been rearranged and shouldn't be trusted
        for (_, improper_match) in improper_matches.items():
        #for (atom_indices, improper_match) in improper_matches.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            # For impropers, central atom is atom 1
            # for (i, j) in [(0, 1), (1, 2), (1, 3)]:
            #     topology.assert_bonded(atom_indices[i], atom_indices[j])
            self._assert_correct_connectivity(improper_match, [(0, 2), (1, 2), (2, 3)])
            atom_indices = improper_match.environment_match.topology_atom_indices
            improper = improper_match.parameter_type

            # TODO: This is a lazy hack. idivf should be set according to the ParameterHandler's default_idivf attrib
            if improper.idivf is None:
                improper.idivf = [3 for item in improper.k]
            # Impropers are applied in three paths around the trefoil having the same handedness
            for (improper_periodicity, improper_phase, improper_k, improper_idivf) in zip(improper.periodicity,
                                               improper.phase, improper.k, improper.idivf):
                # TODO: Implement correct "auto" behavior
                if improper_idivf == 'auto':
                    improper_idivf = 3
                    #logger.warning("The OpenForceField toolkit hasn't implemented "
                    #               "support for the torsion `idivf` value of 'auto'."
                    #               "Currently assuming a value of '3' for impropers.")
                # Permute non-central atoms
                #others = [atom_indices[0], atom_indices[2], atom_indices[3]]
                ## ((0, 1, 2), (1, 2, 0), and (2, 0, 1)) are the three paths around the trefoil
                #for p in [(others[i], others[j], others[k]) for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]]:
                #    # The torsion force gets added three times, since the k is divided by three
                #    force.addTorsion(atom_indices[1], p[0], p[1], p[2],
                #                     improper_periodicity, improper_phase, improper_k/improper_idivf)
                force.addTorsion(atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3],
                                 improper_periodicity, improper_phase, improper_k / improper_idivf)

        #logger.info(
        #    '{} AMBER-style impropers added, each applied strictly in the order that the SMIRKS matched the central atom'.format(
        #        len(improper_matches)))

