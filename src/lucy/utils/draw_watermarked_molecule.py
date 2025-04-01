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

def draw_watermarked_molecule(molecule, rdkit_coords=True) -> BytesIO:
    rdDepictor.SetPreferCoordGen(rdkit_coords)

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
