''' get_molecule_name.py  The purpose of this program is to reverse the conversion back to a molecule name from cd ../.
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

def get_molecule_name(molecule) -> str:
    try:
        logger.info("Starting to retrieve molecule name from the given molecule.")

        # Convert the RDKit molecule to a SMILES string
        smiles = Chem.MolToSmiles(molecule)
        if not smiles:
            logger.warning("Failed to convert the molecule to a SMILES string.")
            return 'Unknown'

        logger.debug(f"Generated SMILES string: {smiles}")

        # Query PubChem for compound data using the SMILES string
        compounds = pcp.get_compounds(smiles, 'smiles')
        if not compounds:
            logger.warning("No compounds found for the given SMILES string.")
            raise ValueError('No compound found for the given SMILES string')

        compound_data = compounds[0].to_dict(properties=['synonyms'])
        logger.debug(f"Retrieved compound data: {compound_data}")

        # Return the first synonym as the molecule name
        if 'synonyms' in compound_data and compound_data['synonyms']:
            molecule_name = compound_data['synonyms'][0]
            logger.info(f"Molecule name retrieved: {molecule_name}")
            return molecule_name
        else:
            logger.warning("No synonyms found in compound data.")
            return 'Unknown'

    except Exception as e:
        logger.error(f"An error occurred while retrieving the molecule name: {e}")
        return 'Unknown'
