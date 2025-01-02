'''  get_proximity.py  The purpose of this program is to calculate the Tanimoto similarity between two molecules cd ../.
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
from rdkit.Chem import AllChem, DataStructs
from utils.setup_logging import logger

def get_proximity(default, input) -> float:
    try:
        logger.info('Starting proximity calculation between molecules.')

        # Generate Morgan fingerprints for the default molecule
        logger.debug('Generating Morgan fingerprint for the default molecule.')
        default_fp = AllChem.GetMorganFingerprintAsBitVect(default, 2)
        if default_fp is None:
            logger.error('Failed to generate fingerprint for the default molecule.')
            raise ValueError('Invalid default molecule.')

        # Generate Morgan fingerprints for the input molecule
        logger.debug('Generating Morgan fingerprint for the input molecule.')
        input_fp = AllChem.GetMorganFingerprintAsBitVect(input, 2)
        if input_fp is None:
            logger.error('Failed to generate fingerprint for the input molecule.')
            raise ValueError('Invalid input molecule.')

        # Calculate fingerprint similarity
        logger.debug('Calculating similarity between the fingerprints.')
        similarity = DataStructs.FingerprintSimilarity(default_fp, input_fp)
        logger.info(f'Calculated similarity: {similarity}')

        return similarity

    except Exception as e:
        logger.error(f'An error occurred during proximity calculation: {e}')
        raise
