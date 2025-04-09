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
from rdkit.Geometry import Point3D

import numpy as np

def draw_watermarked_molecule(molecule, rotation=0, rdkit_bool=True) -> BytesIO:
    rdDepictor.SetPreferCoordGen(rdkit_bool)
    try:
        resolved_name = get_molecule_name(molecule)
        d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
        Options = d2d.drawOptions()
        Options.prepareMolsBeforeDrawing = False
        Options.includeMetadata = False
        Options.bondLineWidth = 4.0
        d2d.SetDrawOptions(Options)
        rdMolDraw2D.SetDarkMode(Options)
        mol = rdMolDraw2D.PrepareMolForDrawing(molecule, kekulize=True)
        mol = rotate_molecule(mol, rotation)
        mol.UpdatePropertyCache(False)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        drawing = d2d.GetDrawingText()
        output = add_watermark(BytesIO(drawing), resolved_name, False)
        return output
    except Exception as e:
        logger.error(f'An error occurred while drawing the watermarked molecule: {e}')
        raise

def rotate_molecule(mol, angle):
    conf = mol.GetConformer()
    rad_angle = float(np.radians(angle))
    cos_theta = float(np.cos(rad_angle))
    sin_theta = float(np.sin(rad_angle))
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x_new = float(cos_theta * pos.x - sin_theta * pos.y)
        y_new = float(sin_theta * pos.x + cos_theta * pos.y)
        z = float(pos.z)
        conf.SetAtomPosition(i, Point3D(x_new, y_new, z))
    return mol
