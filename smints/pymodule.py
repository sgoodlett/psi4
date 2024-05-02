#
# @BEGIN LICENSE
#
# smints by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util
import smints
import numpy as np
import molsym

def run_smints(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    smints can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('smints')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.core.set_local_option('MYPLUGIN', 'PRINT', 1)

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)
    mol = ref_wfn.molecule()
    molschema = mol.to_schema('psi4')
    
    molsym_mol = molsym.Molecule.from_psi4_schema(molschema)
    symtext = molsym.Symtext.from_molecule(molsym_mol)
    molsym_basis_set = [[ref_wfn.basisset().shell(i,j).am for j in range(ref_wfn.basisset().nshell_on_center(i))] for i in range(len(molsym_mol))]
    fxn_set = molsym.salcs.SphericalHarmonics(symtext, molsym_basis_set)
    salcs = molsym.salcs.ProjectionOp(symtext, fxn_set)
    # SALCs, atom map, mult table, irrep dims, irrep mats, fxn_map, ctable
    p4salcs = psi4.core.Matrix.from_array(salcs.basis_transformation_matrix)
    p4irrepdims = psi4.core.Dimension.from_list([d.d for d in salcs.irreps])
    p4mult_table = psi4.core.Dimension.from_list(symtext.mult_table.flatten())
    p4atom_map = psi4.core.Dimension.from_list(symtext.atom_map.flatten())
    p4character_table = psi4.core.Matrix.from_array(symtext.character_table)
    big_matrix_list = []
    for irrep in salcs.irreps:
        for i in range(len(symtext)):
            big_matrix_list.append(symtext.irrep_mat[irrep.symbol][i])
    p4irrep_mat = psi4.core.Matrix.from_array(big_matrix_list)
    big_matrix_list2 = []
    for i in range(len(salcs)):
        big_matrix_list2.append(fxn_set.fxn_map[i])
    p4fxn_map = psi4.core.Matrix.from_array(big_matrix_list2)
    print(p4mult_table.to_tuple())
    print(p4atom_map.to_tuple())
    print(p4irrep_mat.nph)
    print(p4fxn_map.nph)

    smints_G = smints.smints_benbee(ref_wfn, p4salcs, p4irrepdims, p4mult_table, p4atom_map, p4character_table, p4irrep_mat, p4fxn_map)
    print(smints_G)
    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    smints_wfn = psi4.core.plugin('smints.so', ref_wfn)

    #smints_wfn = psi4.core.plugin('smints.so', ref_wfn, p4salcs, p4irrepdims, p4mult_table, p4atom_map, p4character_table, p4irrep_mat, p4fxn_map)

    return smints_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['smints'] = run_smints


def exampleFN():
    # Your Python code goes here
    pass
