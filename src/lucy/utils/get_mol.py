''' get_mol.py  The purpose of this program is fetch a molecule SMILES from pubchem from a name from cd ../.
    Copyright (C) 2024  github.com/brandongrahamcobb

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
from rdkit import Chem
from utils.setup_logging import logger

import pubchempy as pcp

def get_mol(arg):
    try:
        logger.info(f"Attempting to retrieve molecule for argument: {arg}")

        # Query PubChem for compounds by name
        compounds = pcp.get_compounds(arg, 'name')
        if not compounds:
            logger.warning(f"No compounds found for argument: {arg}")
            raise ValueError("Compound not found in PubChem")

        compound_data = compounds[0].to_dict(properties=['isomeric_smiles'])
        logger.debug(f"Compound data retrieved: {compound_data}")

        # Convert SMILES string to RDKit molecule
        mol = Chem.MolFromSmiles(compound_data['isomeric_smiles'])
        if mol is None:
            logger.error("Failed to generate molecule from SMILES string")
            raise ValueError('Invalid SMILES string')

        logger.info("Molecule successfully generated from PubChem data")
        return mol

    except Exception as e:
        logger.warning(f"PubChem retrieval failed for argument: {arg}, attempting direct SMILES conversion. Error: {e}")

        # Attempt to directly convert the argument into an RDKit molecule
        mol = Chem.MolFromSmiles(arg)
        if mol is None:
            logger.error(f"Failed to generate molecule from direct SMILES conversion for argument: {arg}")
            raise ValueError('Invalid SMILES string')

        logger.info("Molecule successfully generated from direct SMILES conversion")
        return mol

