''' get_mol.py  The purpose of this program is fetch a molecule SMILES from pubchem from a name.
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
from pyPept.sequence import Sequence, correct_pdb_atoms
from pyPept.molecule import Molecule
from pyPept.converter import Converter
from rdkit import Chem

import pubchempy as pcp
import re

def three_to_fasta(three_letter_seq):
    """Convert a three-letter peptide sequence (e.g., ArgGlyAla) to FASTA format (e.g., RGA)."""
    three_to_one = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
    }

    # Ensure the sequence is split correctly (handles CamelCase format)
    tokens = re.findall(r'[A-Z][a-z]{2}', three_letter_seq)
    
    # Convert to FASTA
    fasta_seq = ''.join(three_to_one.get(aa, '?') for aa in tokens)

    if '?' in fasta_seq:
        return None  # Return None if any invalid amino acids are found

    return fasta_seq

def fasta_to_helm(fasta_seq):
    """Convert FASTA sequence (AGP) to HELM format for pyPept."""
    helm = f"PEPTIDE1{{{'.'.join(fasta_seq)}}}$$$$V2.0"  # HELM format expected by pyPept
    return helm

def get_mol(arg):
    """Retrieve an RDKit molecule object from PubChem, SMILES, or a peptide sequence."""
    try:
        logger.info(f"🔍 Attempting to retrieve molecule for: {arg}")

        # 1️⃣ Try PubChem Lookup
        compounds = pcp.get_compounds(arg, 'name')
        if compounds:
            smiles = compounds[0].to_dict(properties=['isomeric_smiles']).get('isomeric_smiles')
            if smiles:
                logger.info(f"🧪 PubChem SMILES: {smiles}")
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    logger.info(f"✅ Molecule '{arg}' retrieved from PubChem")
                    return mol
        
        logger.warning(f"⚠️ No valid compound found for '{arg}' in PubChem.")

        # 2️⃣ Try Direct SMILES Conversion
        mol = Chem.MolFromSmiles(arg)
        if mol:
            logger.info(f"✅ Molecule '{arg}' generated from direct SMILES input.")
            return mol

        # 3️⃣ Convert Three-Letter Peptide to FASTA
        fasta_seq = three_to_fasta(arg)
        if fasta_seq:
            logger.info(f"🔄 Converted '{arg}' to FASTA sequence: {fasta_seq}")
            helm_seq = fasta_to_helm(fasta_seq)  # Convert to HELM
            logger.info(f"🔄 Converted FASTA '{fasta_seq}' to HELM '{helm_seq}'")

        # 4️⃣ Extract SMILES from pyPept Instead of Using RDKit
        try:
            conv = Converter(helm=helm_seq)  # ✅ Use HELM format instead of FASTA
            biln_seq = conv.get_biln()
            seq = Sequence(biln_seq)
            seq = correct_pdb_atoms(seq)

            mol_obj = Molecule(seq)  # Create Molecule object

            # ✅ Get SMILES from pyPept
            smiles = mol_obj.get_molecule(fmt="Smiles")
            if smiles:
                logger.info(f"🧪 pyPept SMILES: {smiles}")  # Debugging SMILES output
                rdkit_mol = Chem.MolFromSmiles(smiles)
                if rdkit_mol:
                    logger.info(f"✅ Peptide '{arg}' successfully converted to RDKit molecule using pyPept SMILES.")
                    return rdkit_mol
                else:
                    logger.error(f"❌ RDKit failed to generate molecule from pyPept SMILES '{smiles}'.")
            else:
                logger.error(f"❌ pyPept did not return a SMILES string for peptide '{arg}'.")

        except Exception as peptide_error:
            logger.error(f"❌ Failed to convert peptide sequence '{arg}': {peptide_error}")

        # 5️⃣ If Everything Fails, Return None
        logger.error(f"❌ '{arg}' is not a valid chemical name, SMILES, or peptide sequence.")
        return None

    except Exception as e:
        logger.error(f"🚨 Failed to generate molecule for '{arg}'. Error: {e}")
        return None
