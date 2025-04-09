''' draw_fingerprint.py  The purpose of this program is to draw rdkit generated graphs of the MorganFingerprint for combination and presentation.
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
from lucy.utils.adjust_hue_and_saturation import adjust_hue_and_saturation
from lucy.utils.get_molecule_name import get_molecule_name
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, Draw, rdMolAlign, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.Geometry import Point3D

import numpy as np

def draw_fingerprint(pair, rdkit_bool=True, rotation=0) -> BytesIO:
    rdDepictor.SetPreferCoordGen(rdkit_bool)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    d2d.prepareMolsBeforeDrawing = False
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    mol1, mol2 = pair[0], pair[1]
    mol1, mol2 = standardize_molecule(mol1), standardize_molecule(mol2)
    mol1 = rotate_molecule(mol1, rotation)
    mol2 = rotate_molecule(mol2, rotation)
    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
        mol1, mol2,
        lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192),
        draw2d=d2d, drawingOptions=Options, forceCoords=True
    )
    name = get_molecule_name(mol2)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = BytesIO(drawing)
    output.seek(0)
    img = Image.open(output).convert('RGBA')
    inverted_img = Image.eval(img, lambda x: 255 - x)
    adjusted_output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)
    adjusted_output.seek(0)
    watermarked_image_buffer = add_watermark(adjusted_output, name, False)
    return watermarked_image_buffer

def standardize_molecule(mol):
    mol = Chem.Mol(mol)
    Chem.rdDepictor.Compute2DCoords(mol)
    return mol

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

