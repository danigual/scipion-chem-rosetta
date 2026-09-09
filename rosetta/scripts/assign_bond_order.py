#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Assigns correct bond orders to a ligand read from a PDB file, using a reference SMILES as
template (RDKit's AssignBondOrdersFromTemplate). PDB HETATM records don't carry bond
order information, so a naive PDB read guesses bonds from geometry alone and often gets
them wrong (e.g. double vs single bonds); the SMILES acts as ground truth.

This is a standalone script, not an importable function, because it must run with the
Python interpreter of pwchem's dedicated RDKit conda env (RDKIT_DIC in
pwchem/constants.py) - RDKit is not installed in Scipion's own Python environment.
Invoked via pwchem.Plugin.runScript() from protocol_protac_model.py.

Usage: python assign_bond_order.py <input.pdb> <smiles> <output.sdf>
"""
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

if __name__ == '__main__':
    inputPdb, smiles, outputSdf = sys.argv[1], sys.argv[2], sys.argv[3]

    template = Chem.MolFromSmiles(smiles)
    mol = Chem.MolFromPDBFile(inputPdb)
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)

    writer = Chem.SDWriter(outputSdf)
    writer.write(mol)
    writer.close()
