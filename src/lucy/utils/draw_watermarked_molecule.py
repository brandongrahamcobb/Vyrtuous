''' draw_watermarked_molecule.py  The purpose of this program is to generate an rdkit drawing with a watermark.
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
from io import BytesIO
from lucy.utils.add_watermark import add_watermark
from lucy.utils.get_molecule_name import get_molecule_name
from lucy.utils.setup_logging import logger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

import numpy as np

def draw_watermarked_molecule(molecule, rotation=0, rdkit_bool=True) -> BytesIO:
    rdDepictor.SetPreferCoordGen(rdkit_bool)
    try:
        logger.info('Starting to draw watermarked molecule.')
        resolved_name = get_molecule_name(molecule)
        logger.debug(f'Resolved molecule name: {resolved_name}')
        d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
        Options = d2d.drawOptions()
        Options.prepareMolsBeforeDrawing = False
        Options.includeMetadata = False
        Options.bondLineWidth = 4.0
        d2d.SetDrawOptions(Options)
        logger.debug('Drawing options configured.')
        rdMolDraw2D.SetDarkMode(Options)
        mol = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=True)
        mol = rotate_molecule(mol, rotation)
        mol.UpdatePropertyCache(False)
        logger.debug('Molecule prepared for drawing.')
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        logger.info('Molecule drawing completed.')
        drawing = d2d.GetDrawingText()
        output = add_watermark(BytesIO(drawing), resolved_name, False)
        logger.info('Watermark added successfully.')
        return output

    except Exception as e:
        logger.error(f'An error occurred while drawing the watermarked molecule: {e}')
        raise

from rdkit.Geometry import Point3D

def rotate_molecule(mol, angle):
    """Rotate a molecule by a given angle (in degrees)."""
    conf = mol.GetConformer()
    rad_angle = float(np.radians(angle))  # Convert degrees to radians
    cos_theta = float(np.cos(rad_angle))
    sin_theta = float(np.sin(rad_angle))

    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x_new = float(cos_theta * pos.x - sin_theta * pos.y)
        y_new = float(sin_theta * pos.x + cos_theta * pos.y)
        z = float(pos.z)
        conf.SetAtomPosition(i, Point3D(x_new, y_new, z))

    return mol
