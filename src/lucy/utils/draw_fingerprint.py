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
#from io import BytesIO
#from lucy.utils.setup_logging import logger
#from rdkit.Chem import rdFingerprintGenerator
#from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
#
#def draw_fingerprint(pair) -> BytesIO:
#    # Create a Morgan Fingerprint Generator with specified parameters
#    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, countSimulation=True)
#
#    # Define a function to get the fingerprint using mfpgen
#    def get_fp(mol, *args, **kwargs):
#        return mfpgen.GetFingerprint(mol)
#
#    # Set up the drawing
#    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
#    d2d.prepareMolsBeforeDrawing = False
#    Options = d2d.drawOptions()
#    Options.prepareMolsBeforeDrawing = False
#    Options.includeMetadata = False
#    Options.bondLineWidth = 4.0
#    d2d.SetDrawOptions(Options)
#
#    # Prepare molecules for drawing
#    mol1 = rdMolDraw2D.PrepareMolForDrawing(pair[0], kekulize=True)
#    mol1.UpdatePropertyCache(False)
#    mol2 = rdMolDraw2D.PrepareMolForDrawing(pair[1], kekulize=True)
#    mol2.UpdatePropertyCache(False)
#
#    # Get similarity map for fingerprints
#    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
#        mol1, 
#        mol2, 
#        get_fp, 
#        draw2d=d2d, 
#        drawingOptions=Options
#    )
#
#    # Finish drawing and return the image as a BytesIO object
#    d2d.FinishDrawing()
#    drawing = d2d.GetDrawingText()
#    output = BytesIO(drawing)
#    return output

from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, Draw, rdMolAlign, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
import numpy as np
from io import BytesIO
from PIL import Image
from lucy.utils.add_watermark import add_watermark
from lucy.utils.adjust_hue_and_saturation import adjust_hue_and_saturation
from lucy.utils.get_molecule_name import get_molecule_name

def draw_fingerprint(pair, rotation=0, rdkit_coords=True) -> BytesIO:
    rdDepictor.SetPreferCoordGen(rdkit_coords)
    """Draws similarity maps for a molecular pair, aligning based on a core structure and allowing individual rotation."""
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    d2d.prepareMolsBeforeDrawing = False
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)

    mol1, mol2 = pair[0], pair[1]

    # Align molecules if a core is provided
    mol1, mol2 = standardize_molecule(mol1), standardize_molecule(mol2)
    # Apply independent rotations
    rotate_molecule(mol1, rotation)
    rotate_molecule(mol2, rotation)

    # Compute similarity map
    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
        mol1, mol2, 
        lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), 
        draw2d=d2d, drawingOptions=Options
    )

    name = get_molecule_name(mol2)

    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = BytesIO(drawing)
    output.seek(0)
    img = Image.open(output).convert("RGBA")
    inverted_img = Image.eval(img, lambda x: 255 - x)  # Invert colors
    adjusted_output = adjust_hue_and_saturation(inverted_img, hue_shift=-180, saturation_shift=160)  
    adjusted_output.seek(0)
    watermarked_image_buffer = add_watermark(adjusted_output, name, True)
    
    return watermarked_image_buffer

def standardize_molecule(mol):
    """Generate 2D coordinates for a molecule to ensure consistent orientation."""
    mol = Chem.Mol(mol)  # Create a copy
    Chem.rdDepictor.Compute2DCoords(mol)  # Ensure 2D coordinates
    return mol

def rotate_molecule(mol, angle):
    """Rotate a molecule by a given angle (in degrees)."""
    conf = mol.GetConformer()
    rad_angle = np.radians(angle)  # Convert degrees to radians
    cos_theta, sin_theta = np.cos(rad_angle), np.sin(rad_angle)

    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x_new = cos_theta * pos.x - sin_theta * pos.y
        y_new = sin_theta * pos.x + cos_theta * pos.y
        conf.SetAtomPosition(i, (x_new, y_new, pos.z))

    return mol
