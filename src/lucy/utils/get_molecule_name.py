''' get_molecule_name.py  The purpose of this program is to reverse the conversion back to a molecule name.
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
from lucy.utils.setup_logging import logger
from rdkit import Chem

import pubchempy as pcp

def get_molecule_name(molecule) -> str:
    try:
        smiles = Chem.MolToSmiles(molecule)
        if not smiles:
            return 'Unknown'
        compounds = pcp.get_compounds(smiles, 'smiles')
        if not compounds:
            raise ValueError('No compound found for the given SMILES string')
        compound_data = compounds[0].to_dict(properties=['synonyms'])
        if 'synonyms' in compound_data and compound_data['synonyms']:
            molecule_name = compound_data['synonyms'][0]
            return molecule_name
        else:
            return 'Unknown'
    except Exception as e:
        logger.error(f"An error occurred while retrieving the molecule name: {e}")
        return 'Unknown'
