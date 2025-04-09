'''  get_proximity.py  The purpose of this program is to calculate the Tanimoto similarity between two molecules.
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
from rdkit.Chem import AllChem, DataStructs

def get_proximity(default, input) -> float:
    try:
        default_fp = AllChem.GetMorganFingerprintAsBitVect(default, 2)
        if default_fp is None:
            raise ValueError('Invalid default molecule.')
        input_fp = AllChem.GetMorganFingerprintAsBitVect(input, 2)
        if input_fp is None:
            raise ValueError('Invalid input molecule.')
        similarity = DataStructs.FingerprintSimilarity(default_fp, input_fp)
        return similarity
    except Exception as e:
        logger.error(f'An error occurred during proximity calculation: {e}')
        raise
