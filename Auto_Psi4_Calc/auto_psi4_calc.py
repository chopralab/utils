# Portions taken from https://iwatobipen.wordpress.com/2018/08/27/calculate-homo-and-lumo-with-psi4-reviced-rdkit-psi4/

# Import statements
import pandas as pd
import psi4
from psi4.driver.qcdb import molecule
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser()

parser.add_argument('--infile', type=str, default='infile.txt', help='Name of the input file')
parser.add_argument('--outfile', type=str, default='outfile.csv', help='Name of the output file')
parser.add_argument('--corefile', type=str, default='output.dat', help='Name of the Psi4 core output file')

parser.add_argument('--energy', action='store_true', default=False, help='Argument flag for calculating energy info')
parser.add_argument('--homolumo', action='store_true', default=False, help='Argument flag for calculating homo and lumo info')
parser.add_argument('--geomopt', action='store_true', default=False, help='Arugment flag for calculation of geometry optimization')
parser.add_argument('--vibfreq', action='store_true', default=False, help='Argument flag for calculation of vibrational frequency')

parser.add_argument('--method', type=str, default='HF', help='Method used for energy calculation (i.e. Hartree Fock)')
parser.add_argument('--basis', type=str, default='6-31G*', help='Basis set used for energy calculation (i.e. 6-31G*)')

args = parser.parse_args()

# Method for converting a RDKit mol object to a xyz coordinate file
def mol_to_xyz(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
    AllChem.UFFOptimizeMolecule(mol)
    atoms = mol.GetAtoms()
    string = string = "\n"
    for i, atom in enumerate(atoms):
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
    string += "units angstrom\n"
    return string, mol

# Method for calculating molecule energy, HOMO, and LUMO
def calc_homo_lumo(smiles):
    mol = Chem.MolFromSmiles(smiles)
    xyz, mol=mol_to_xyz(mol)
    psi4.set_memory('4 GB')
    psi4.set_num_threads(4)
    geom = psi4.geometry(xyz)
    method_basis = args.method + '/' + args.basis
    e, wfn = psi4.energy(method_basis, return_wfn=True)
    homo = wfn.epsilon_a_subset('AO', 'ALL').np[wfn.nalpha()-1]
    lumo = wfn.epsilon_a_subset('AO', 'ALL').np[wfn.nalpha()]
    print("SMILES: " + smiles + ", Energy: " + str(e) + " H, HOMO: " + str(homo) + " H, LUMO: " + str(lumo) + " H")
    return homo, lumo, e

# Method for calculation of geometry optimization 
def opt_geom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    xyz, mol=mol_to_xyz(mol)
    psi4.set_memory('4 GB')
    psi4.set_num_threads(4)
    method_basis = args.method + '/' + args.basis
    geom = psi4.geometry(xyz)
    opt_e = psi4.optimize(method_basis, molecule=geom)
    print("SMILES: " + smiles + ", Optimization Energy: " + str(opt_e) + " H")
    return opt_e
    
# Method for vibrational frequency calculcation
def calc_vib_freq(smiles):
    mol = Chem.MolFromSmiles(smiles)
    xyz, mol=mol_to_xyz(mol)
    psi4.set_memory('4 GB')
    psi4.set_num_threads(4)
    method_basis = args.method + '/' + args.basis
    geom = psi4.geometry(xyz)
    e, wfn = psi4.frequency(method_basis, molecule=geom, return_wfn=True)
    print("SMILES: " + smiles + ", Vibrational Frequency Energy: " + str(e))
    return e

# Main method for calculating energy and HL information for all molecules in the input file
def main():
    psi4.core.set_output_file(args.corefile, True)

    with open(args.infile, 'r') as infile:
        cols = ['smiles']
        if args.energy:
            cols.append('energy')
        if args.homolumo:
            cols.append('homo')
            cols.append('lumo')
        if args.geomopt:
            cols.append('geom_opt_energy')
        if args.vibfreq:
            cols.append('vib_freq_energy')

        df  = pd.DataFrame(columns=cols)

        for smiles in infile:
            row = [smiles]
            hl_flag = args.homolumo
            if args.energy:
                homo, lumo, e = calc_homo_lumo(smiles)
                row.append(e)
                if hl_flag:
                    row.append(homo)
                    row.append(lumo)
                    hl_flag = False
            if hl_flag:
                homo, lumo, e = calc_homo_lumo(smiles)
                row.append(homo)
                row.append(lumo)
            if args.geomopt:
                opt_e = opt_geom(smiles)
                row.append(opt_e)
            if args.vibfreq:
                e = calc_vib_freq(smiles)
                row.append(e)
        
        df.to_csv(args.outfile, index=False, sep=',')
    infile.close()

# Run the main method
if __name__ == "__main__":
    main()