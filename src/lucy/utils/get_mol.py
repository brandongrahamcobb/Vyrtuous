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
from rdkit.Chem.MolStandardize import rdMolStandardize

def parse_helm_for_residues(helm_str):
    m = re.search(r'PEPTIDE1\{([^}]+)\}', helm_str)
    if not m:
        raise ValueError("Could not parse HELM polymer definition.")
    residue_str = m.group(1)
    residues = residue_str.split('.')
    return residues

def clean_peptide_sequence(seq):
    seq = seq.replace("(", "").replace(")", "").replace(" ", "").strip()
    if len(seq) % 3 != 0:
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
        return None
    return fasta

def reverse_peptide_sequence(seq):
    cleaned = seq.replace("(", "").replace(")", "").strip()
    tokens = [cleaned[i:i+3] for i in range(0, len(cleaned), 3)]
    reversed_tokens = list(reversed(tokens))
    return "".join(reversed_tokens)

def construct_helm_from_peptide(seq):
    fasta = three_to_fasta(seq)
    if not fasta:
        return None
    helm = construct_linear_helm(fasta)
    logger.info(f"Constructed HELM: {helm}")
    return helm

def construct_linear_helm(fasta_str):
    polymer = "PEPTIDE1{" + ".".join(list(fasta_str)) + "}"
    sections = [polymer, "", "", "", "V2.0"]
    helm = "$".join(sections)
    return helm

def construct_helm(sequence, reverse=False):
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
    emol = Chem.RWMol(mol)
    dummy_idxs = [atom.GetIdx() for atom in emol.GetAtoms() if atom.GetAtomicNum() == 0]
    for idx in sorted(dummy_idxs, reverse=True):
        emol.RemoveAtom(idx)
    return emol.GetMol()

def helm_to_smiles_manual(helm_str):
    m = re.search(r'PEPTIDE1\{([^}]+)\}', helm_str)
    if not m:
        raise ValueError("Could not parse HELM polymer definition.")
    residue_str = m.group(1)
    residues = residue_str.split('.')
    mol = build_peptide_from_residues(list(residues), cyclic=False)
    smiles = Chem.MolToSmiles(mol, canonical=False)
    return smiles

def get_linear_peptide_mol(sequence, reverse=False):
    helm = construct_helm(sequence, reverse=reverse)
    if helm is None:
        return None
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
    try:
        mol = Chem.MolFromSmiles(arg)
        if mol:
            return mol
    except Exception as e:
        logger.warning(f"Direct SMILES conversion failed for '{arg}': {e}")
    try:
        compounds = pcp.get_compounds(arg, 'name')
        if compounds:
            smiles = compounds[0].to_dict(properties=['isomeric_smiles']).get('isomeric_smiles')
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return mol
    except Exception as e:
        logger.warning(f"PubChem lookup failed for '{arg}': {e}")
    if any(x in arg for x in ["Ala", "Arg", "His", "Lys", "Pro", "Gly", "Cys"]):
        try:
            return get_linear_peptide_mol(arg, reverse=reverse)
        except Exception as e:
            logger.error(f"Peptide conversion error: {e}")
            return None
    return None


def standardize_smiles(smiles_str):
    try:
        std_smiles = rdMolStandardize.StandardizeSmiles(smiles_str)
        return std_smiles
    except Exception as e:
        logger.error(f"Error standardizing SMILES: {e}")
        return smiles_str

def manual_helm_to_smiles(helm_str):
    try:
        residues = parse_helm_for_residues(helm_str)
    except Exception as e:
        logger.error(f"Error parsing HELM: {e}")
        return None
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
    connector = "NC(=O)"
    smiles_parts = []
    for r in residues:
        frag = frag_dict.get(r)
        if not frag:
            logger.error(f"Residue {r} not defined in frag_dict.")
            return None
        smiles_parts.append(frag)
    manual_smiles = connector.join(smiles_parts)
    return manual_smiles
