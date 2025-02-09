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
rdDepictor.SetPreferCoordGen(True)
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

from io import BytesIO
from lucy.utils.setup_logging import logger
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps

def draw_fingerprint(pair) -> BytesIO:
#    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, countSimulation=True)
#    def get_fp(mol, *args, **kwargs):
#        return mfpgen.GetFingerprint(mol)
    d2d = rdMolDraw2D.MolDraw2DCairo(1024, 1024)
    d2d.prepareMolsBeforeDrawing = False
    Options = d2d.drawOptions()
    Options.prepareMolsBeforeDrawing = False
    Options.includeMetadata = False
    Options.bondLineWidth = 4.0
    d2d.SetDrawOptions(Options)
    mol1 = rdMolDraw2D.PrepareMolForDrawing(pair[0], kekulize=True)
    mol1.UpdatePropertyCache(False)
    mol2 = rdMolDraw2D.PrepareMolForDrawing(pair[1], kekulize=True)
    mol2.UpdatePropertyCache(False)
    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(mol1, mol2, lambda m, i: SimilarityMaps.GetMorganFingerprint(m, i, radius=2, fpType='bv', nBits=8192), draw2d=d2d, drawingOptions=Options)
#    fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(pair[1], pair[0], get_fp, draw2d=d2d) #colorMap=brighter_color, draw2d=d2d)
    d2d.FinishDrawing()
    drawing = d2d.GetDrawingText()
    output = BytesIO(drawing)
    return output
