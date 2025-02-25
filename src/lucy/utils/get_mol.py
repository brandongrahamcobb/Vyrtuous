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
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
import pubchempy as pcp
import re
import shlex
import os
import discord
from discord.ext import commands

# ------------------------
# Helper Functions: Peptide Cleaning and Conversion
# ------------------------

def parse_helm_for_residues(helm_str):
    """
    Given a HELM string (assumed linear), extract the residue list.
    For example, from "PEPTIDE1{H.P.K}$$$$V2.0" return ["H", "P", "K"].
    """
    m = re.search(r'PEPTIDE1\{([^}]+)\}', helm_str)
    if not m:
        raise ValueError("Could not parse HELM polymer definition.")
    residue_str = m.group(1)  # e.g. "H.P.K"
    residues = residue_str.split('.')
    return residues

def clean_peptide_sequence(seq):
    """
    Clean a peptide sequence (three-letter codes) by:
      - Removing whitespace and parentheses.
      - Ensuring the length is a multiple of 3.
      - Verifying that every 3-letter token is one of the allowed residues.
    
    Allowed residues: Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys,
    Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val.
    
    Returns the cleaned sequence (with no spaces/parentheses) or None if invalid.
    """
    seq = seq.replace("(", "").replace(")", "").replace(" ", "").strip()
    if len(seq) % 3 != 0:
        logger.error(f"Peptide sequence length is not a multiple of 3: '{seq}'")
        return None
    allowed = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
               "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",
               "Thr", "Trp", "Tyr", "Val"}
    tokens = [seq[i:i+3] for i in range(0, len(seq), 3)]
    for token in tokens:
        if token not in allowed:
            logger.error(f"Unknown residue: '{token}' in sequence '{seq}'")
            return None
    return "".join(tokens)

def three_to_fasta(three_letter_seq):
    """
    Convert a three-letter peptide sequence (e.g. "HisProLys") to a one-letter FASTA string.
    Example: "HisProLys" -> "HPK"
    """
    cleaned = clean_peptide_sequence(three_letter_seq)
    if not cleaned:
        return None
    three_to_one = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
    }
    tokens = [cleaned[i:i+3] for i in range(0, len(cleaned), 3)]
    fasta = ''.join(three_to_one.get(token, '?') for token in tokens)
    if '?' in fasta:
        logger.error(f"Unknown residue in sequence: '{three_letter_seq}'")
        return None
    return fasta

def reverse_peptide_sequence(seq):
    """
    Reverse the order of residues in a peptide sequence (three-letter codes).
    Example: "HisProLys" -> "LysProHis"
    """
    cleaned = seq.replace("(", "").replace(")", "").strip()
    tokens = [cleaned[i:i+3] for i in range(0, len(cleaned), 3)]
    reversed_tokens = list(reversed(tokens))
    return "".join(reversed_tokens)

def construct_helm_from_peptide(seq):
    """
    Given a peptide sequence in three-letter codes (assumed linear), convert it to a HELM string.
    """
    fasta = three_to_fasta(seq)
    if not fasta:
        return None
    helm = construct_linear_helm(fasta)
    logger.info(f"Constructed HELM: {helm}")
    return helm

def construct_linear_helm(fasta_str):
    """
    Construct a HELM string for a linear peptide from a one-letter FASTA string.
    For example, if fasta_str is "HPK", returns "PEPTIDE1{H.P.K}$$$$V2.0"
    """
    polymer = "PEPTIDE1{" + ".".join(list(fasta_str)) + "}"
    sections = [polymer, "", "", "", "V2.0"]
    helm = "$".join(sections)
    return helm

def construct_helm(sequence, reverse=False):
    """
    Construct a HELM string from a peptide sequence (assumed linear).
    If reverse is True, the peptide sequence is reversed manually.
    
    Example:
      Input: "HisProLys" with reverse=False yields "PEPTIDE1{H.P.K}$$$$V2.0"
      Input: "HisProLys" with reverse=True yields "PEPTIDE1{K.P.H}$$$$V2.0"
    """
    if reverse:
        sequence = reverse_peptide_sequence(sequence)
        logger.info(f"Reversed peptide sequence: {sequence}")
    fasta_str = three_to_fasta(sequence)
    if not fasta_str:
        logger.error(f"Failed to convert peptide sequence '{sequence}' to FASTA.")
        return None
    helm = construct_linear_helm(fasta_str)
    logger.info(f"Constructed HELM: {helm}")
    return helm

# ------------------------
# Residue Fragments Dictionary (20 Natural Amino Acids)
# These are simplified fragment SMILES with dummy atoms [*:1] (N-terminal) and [*:2] (C-terminal).
residue_fragments = {
    "A": "[*:1]NC(=O)[C@@H](C)[*:2]",
    "R": "[*:1]NC(=O)[C@@H](CCCNC(=N)N)[*:2]",
    "N": "[*:1]NC(=O)[C@@H](CC(=O)N)[*:2]",
    "D": "[*:1]NC(=O)[C@@H](CC(=O)O)[*:2]",
    "C": "[*:1]NC(=O)[C@@H](CSH)[*:2]",
    "Q": "[*:1]NC(=O)[C@@H](CCC(=O)N)[*:2]",
    "E": "[*:1]NC(=O)[C@@H](CCC(=O)O)[*:2]",
    "G": "[*:1]NCC(=O)[*:2]",
    "H": "[*:1]NC(=O)[C@@H](Cc1cnc[nH]1)[*:2]",
    "I": "[*:1]NC(=O)[C@@H](C(C)CC)[*:2]",
    "L": "[*:1]NC(=O)[C@@H](CC(C)C)[*:2]",
    "K": "[*:1]NC(=O)[C@@H](CCCCN)[*:2]",
    "M": "[*:1]NC(=O)[C@@H](CCSC)[*:2]",
    "F": "[*:1]NC(=O)[C@@H](Cc1ccccc1)[*:2]",
    "P": "[*:1]NC(=O)[C@@H](C1CCCN1)[*:2]",
    "S": "[*:1]NC(=O)[C@@H](CO)[*:2]",
    "T": "[*:1]NC(=O)[C@@H]([C@H](O)C)[*:2]",
    "W": "[*:1]NC(=O)[C@@H](Cc1c[nH]c2ccccc12)[*:2]",
    "Y": "[*:1]NC(=O)[C@@H](Cc1ccc(O)cc1)[*:2]",
    "V": "[*:1]NC(=O)[C@@H](C(C)C)[*:2]"
}

def build_peptide_from_residues(res_list, cyclic=False):
    """
    Given a list of one-letter residue codes (e.g. ["H", "P", "K"]),
    build a peptide molecule manually using the residue_fragments dictionary.
    For linear peptides, join the residues sequentially.
    Finally, remove all dummy atoms.
    
    Note: This is a simplified joining process.
    """
    if not res_list:
        raise ValueError("Residue list is empty.")
    mol_list = []
    for r in res_list:
        frag = residue_fragments.get(r)
        if frag is None:
            raise ValueError(f"Residue '{r}' is not defined in residue_fragments.")
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            raise ValueError(f"Failed to parse SMILES for residue '{r}'.")
        mol_list.append(mol)
    current = mol_list[0]
    for i in range(1, len(mol_list)):
        current = join_two_molecules(current, mol_list[i])
    final_mol = remove_dummy_atoms(current)
    return final_mol

def join_two_molecules(mol1, mol2):
    """
    Join two molecules by connecting the dummy atom with map number 2 in mol1
    to the dummy atom with map number 1 in mol2.
    """
    dummy1 = next((atom for atom in mol1.GetAtoms() if atom.GetAtomMapNum() == 2), None)
    if dummy1 is None:
        raise ValueError("No dummy atom with map number 2 found in first molecule.")
    dummy2 = next((atom for atom in mol2.GetAtoms() if atom.GetAtomMapNum() == 1), None)
    if dummy2 is None:
        raise ValueError("No dummy atom with map number 1 found in second molecule.")
    combined = Chem.CombineMols(mol1, mol2)
    n1 = mol1.GetNumAtoms()
    idx1 = dummy1.GetIdx()
    idx2 = n1 + dummy2.GetIdx()
    emol = Chem.EditableMol(combined)
    emol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
    new_mol = emol.GetMol()
    return new_mol

def remove_dummy_atoms(mol):
    """
    Remove all dummy atoms (atomic number 0) from the molecule.
    """
    emol = Chem.RWMol(mol)
    dummy_idxs = [atom.GetIdx() for atom in emol.GetAtoms() if atom.GetAtomicNum() == 0]
    for idx in sorted(dummy_idxs, reverse=True):
        emol.RemoveAtom(idx)
    return emol.GetMol()

def helm_to_smiles_manual(helm_str):
    """
    Convert a HELM string (for a linear peptide) to a SMILES string by manually building the peptide
    molecule using the residue_fragments dictionary.
    """
    # Expect HELM of the form: "PEPTIDE1{H.P.K}$$$$V2.0"
    m = re.search(r'PEPTIDE1\{([^}]+)\}', helm_str)
    if not m:
        raise ValueError("Could not parse HELM polymer definition.")
    residue_str = m.group(1)  # e.g. "H.P.K"
    residues = residue_str.split('.')
    mol = build_peptide_from_residues(list(residues), cyclic=False)
    smiles = Chem.MolToSmiles(mol, canonical=False)
    return smiles

# ------------------------
# Peptide Conversion Function
# ------------------------

def get_linear_peptide_mol(sequence, reverse=False):
    """
    Convert a linear peptide sequence (three-letter codes) into an RDKit molecule using pyPept.
    Optionally reverse the sequence if reverse=True.
    """
    helm = construct_helm(sequence, reverse=reverse)
    if helm is None:
        return None
    logger.info(f"Constructed HELM for linear peptide: {helm}")
    try:
        conv = Converter(helm=helm)
        biln_seq = conv.get_biln()
        seq_obj = Sequence(biln_seq)
        seq_obj = correct_pdb_atoms(seq_obj)
        mol_obj = Molecule(seq_obj)
        rdkit_mol = mol_obj.get_molecule(fmt="ROMol")
        return rdkit_mol
    except Exception as e:
        logger.error(f"Error converting linear peptide: {e}")
        return None

def get_mol(arg, reverse=False):
    """
    Retrieve an RDKit molecule object from the input string.
    
    1. Try PubChem lookup.
    2. Try direct SMILES conversion.
    3. If those fail and if the input appears to be a peptide (contains known peptide abbreviations),
       use peptide conversion via pyPept.
    The reverse flag (default False) reverses the peptide sequence (if applicable).
    """
    logger.info(f"🔍 Attempting to retrieve molecule for: {arg}")
    # 1. PubChem lookup.
    try:
        compounds = pcp.get_compounds(arg, 'name')
        if compounds:
            smiles = compounds[0].to_dict(properties=['isomeric_smiles']).get('isomeric_smiles')
            if smiles:
                logger.info(f"🧪 PubChem SMILES: {smiles}")
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    logger.info(f"✅ Molecule '{arg}' retrieved from PubChem")
                    return mol
    except Exception as e:
        logger.warning(f"PubChem lookup failed for '{arg}': {e}")
    # 2. Direct SMILES conversion.
    try:
        mol = Chem.MolFromSmiles(arg)
        if mol:
            logger.info(f"✅ Molecule '{arg}' generated from direct SMILES input.")
            return mol
    except Exception as e:
        logger.warning(f"Direct SMILES conversion failed for '{arg}': {e}")
    # 3. Assume input is a peptide.
    if any(x in arg for x in ["Ala", "Arg", "His", "Lys", "Pro", "Gly", "Cys"]):
        logger.info("Assuming input is a peptide sequence; using peptide conversion.")
        try:
            return get_linear_peptide_mol(arg, reverse=reverse)
        except Exception as e:
            logger.error(f"Peptide conversion error: {e}")
            return None
    logger.error(f"❌ '{arg}' is not a valid chemical name, SMILES, or peptide sequence.")
    return None

# ------------------------
# Standardization Helper
# ------------------------
from rdkit.Chem.MolStandardize import rdMolStandardize

def standardize_smiles(smiles_str):
    """
    Standardize a SMILES string using RDKit's rdMolStandardize.
    Returns the standardized SMILES or the original if standardization fails.
    """
    try:
        std_smiles = rdMolStandardize.StandardizeSmiles(smiles_str)
        return std_smiles
    except Exception as e:
        logger.error(f"Error standardizing SMILES: {e}")
        return smiles_str
## ------------------------
#
#def three_to_fasta(three_letter_seq):
#    """Convert a three-letter peptide sequence (e.g., ArgGlyAla) to FASTA format (e.g., RGA)."""
#    three_to_one = {
#        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
#        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
#        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
#        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
#    }
#
#    # Ensure the sequence is split correctly (handles CamelCase format)
#    tokens = re.findall(r'[A-Z][a-z]{2}', three_letter_seq)
#    
#    # Convert to FASTA
#    fasta_seq = ''.join(three_to_one.get(aa, '?') for aa in tokens)
#
#    if '?' in fasta_seq:
#        return None  # Return None if any invalid amino acids are found
#
#    return fasta_seq
#
#def fasta_to_helm(fasta_seq):
#    """Convert FASTA sequence (AGP) to HELM format for pyPept."""
#    helm = f"PEPTIDE1{{{'.'.join(fasta_seq)}}}$$$$V2.0"  # HELM format expected by pyPept
#    return helm
#
#def get_mol(arg):
#    """Retrieve an RDKit molecule object from PubChem, SMILES, or a peptide sequence."""
#    try:
#        logger.info(f"🔍 Attempting to retrieve molecule for: {arg}")
#
#        # 1️⃣ Try PubChem Lookup
#        compounds = pcp.get_compounds(arg, 'name')
#        if compounds:
#            smiles = compounds[0].to_dict(properties=['isomeric_smiles']).get('isomeric_smiles')
#            if smiles:
#                logger.info(f"🧪 PubChem SMILES: {smiles}")
#                mol = Chem.MolFromSmiles(smiles)
#                if mol:
#                    logger.info(f"✅ Molecule '{arg}' retrieved from PubChem")
#                    return mol
#        
#        logger.warning(f"⚠️ No valid compound found for '{arg}' in PubChem.")
#
#        # 2️⃣ Try Direct SMILES Conversion
#        mol = Chem.MolFromSmiles(arg)
#        if mol:
#            logger.info(f"✅ Molecule '{arg}' generated from direct SMILES input.")
#            return mol
#
#        # 3️⃣ Convert Three-Letter Peptide to FASTA
#        fasta_seq = three_to_fasta(arg)
#        if fasta_seq:
#            logger.info(f"🔄 Converted '{arg}' to FASTA sequence: {fasta_seq}")
#            helm_seq = fasta_to_helm(fasta_seq)  # Convert to HELM
#            logger.info(f"🔄 Converted FASTA '{fasta_seq}' to HELM '{helm_seq}'")
#
#        # 4️⃣ Extract SMILES from pyPept Instead of Using RDKit
#        try:
#            conv = Converter(helm=helm_seq)  # ✅ Use HELM format instead of FASTA
#            biln_seq = conv.get_biln()
#            seq = Sequence(biln_seq)
#            seq = correct_pdb_atoms(seq)
#
#            mol_obj = Molecule(seq)  # Create Molecule object
#
#            # ✅ Get SMILES from pyPept
#            smiles = mol_obj.get_molecule(fmt="Smiles")
#            if smiles:
#                logger.info(f"🧪 pyPept SMILES: {smiles}")  # Debugging SMILES output
#                rdkit_mol = Chem.MolFromSmiles(smiles)
#                if rdkit_mol:
#                    logger.info(f"✅ Peptide '{arg}' successfully converted to RDKit molecule using pyPept SMILES.")
#                    return rdkit_mol
#                else:
#                    logger.error(f"❌ RDKit failed to generate molecule from pyPept SMILES '{smiles}'.")
#            else:
#                logger.error(f"❌ pyPept did not return a SMILES string for peptide '{arg}'.")
#
#        except Exception as peptide_error:
#            logger.error(f"❌ Failed to convert peptide sequence '{arg}': {peptide_error}")
#
#        # 5️⃣ If Everything Fails, Return None
#        logger.error(f"❌ '{arg}' is not a valid chemical name, SMILES, or peptide sequence.")
#        return None
#
#    except Exception as e:
#        logger.error(f"🚨 Failed to generate molecule for '{arg}'. Error: {e}")
#        return None
def manual_helm_to_smiles(helm_str):
    """
    Convert a HELM string (for a linear peptide) to a SMILES string manually.
    
    For demonstration purposes only, we define a very simplified residue fragment dictionary.
    We assume that the peptide bond is represented by the string "NC(=O)" between residues.
    
    WARNING: This simplistic joining does not yield chemically correct SMILES for real peptides!
    """
    # Parse the HELM string to get the residue sequence.
    try:
        residues = parse_helm_for_residues(helm_str)
    except Exception as e:
        logger.error(f"Error parsing HELM: {e}")
        return None

    # Define a simplified dictionary mapping one-letter codes to fragment strings.
    # These are not real residue SMILES; they are placeholders to illustrate the concept.
    frag_dict = {
        "A": "NCC",   # Placeholder for Ala
        "R": "NC(CCCNC(=N)N)",  # Placeholder for Arg
        "N": "NC(CC(=O)N)",     # Placeholder for Asn
        "D": "NC(CC(=O)O)",     # Asp
        "C": "NC(CS)",          # Cys
        "Q": "NC(CCC(=O)N)",     # Gln
        "E": "NC(CCC(=O)O)",     # Glu
        "G": "NC",              # Gly
        "H": "NC(CC1=CN=CN1)",   # His
        "I": "NC(C(C)CC)",       # Ile
        "L": "NC(CC(C)C)",       # Leu
        "K": "NC(CCCCN)",        # Lys
        "M": "NC(CCSC)",         # Met
        "F": "NC(CC1=CC=CC=C1)",  # Phe
        "P": "NC(C1CCCN1)",       # Pro
        "S": "NC(CO)",           # Ser
        "T": "NC(C(C)O)",        # Thr
        "W": "NC(CC1=CNC2=CC=CC=C12)",  # Trp
        "Y": "NC(CC1=CC=C(O)C=C1)",      # Tyr
        "V": "NC(C(C)C)"         # Val
    }
    # Now, build the SMILES string by joining fragments with a simple peptide-bond connector.
    # For demonstration, assume the peptide bond connector is "NC(=O)".
    connector = "NC(=O)"
    # For the first residue, we use its fragment as-is.
    smiles_parts = []
    for r in residues:
        frag = frag_dict.get(r)
        if not frag:
            logger.error(f"Residue {r} not defined in frag_dict.")
            return None
        smiles_parts.append(frag)
    # Now, join the fragments with the connector.
    # This simplistic method just concatenates the strings.
    manual_smiles = connector.join(smiles_parts)
    # Optionally, add terminal groups (this is left as an exercise).
    logger.info(f"Manually built SMILES: {manual_smiles}")
    return manual_smiles
