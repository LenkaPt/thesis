import argparse
from typing import List, Tuple, TextIO
from structures import Atom, Residue, Chain, Model, Protein
from sdf import skip_non_structural_data


def get_atom_pdb(line: str) -> Atom:
    """Returns one atom from pdb file (as instance of Atom class)"""
    atom_name = line[13:16]
    x, y, z = float(line[31:38]), float(line[39:46]), float(line[47:54])
    atom = Atom(x, y, z, atom_name)
    return atom


def get_residue_pdb(file: TextIO, line: str) -> Tuple[Residue, str, str]:
    """Returns one residuum from pdb file, chain ID and last line."""
    chain_id = line[21]
    residue_seq_number = line[23:26]
    residue_name = line[17:20]
    # collect all atoms of one residue
    atoms = []
    while residue_seq_number == line[23:26]:
        atom = get_atom_pdb(line)
        atoms.append(atom)
        line = file.readline()
        if line.startswith('TER'):
            # TER -> termination of chain (ignored in this program)
            line = file.readline()
            break

    residue = Residue(residue_name, residue_seq_number, atoms)
    return residue, chain_id, line


def get_chain(file: TextIO, line: str) -> Tuple[Chain, str]:
    """Returns one chain of biomolecule and last executed line."""
    residues = []
    last_chain_id = line[21]
    while last_chain_id == line[21]:
        residue, chain_id, line = get_residue_pdb(file, line)
        residues.append(residue)
        last_chain_id = chain_id

    chain = Chain(last_chain_id, residues)
    return chain, line


def get_models(file: TextIO) -> List[Model]:
    """Returns all models in pdb file."""
    models = []
    line = skip_non_structural_data(file, ['MODEL', 'ATOM'])
    # It is possible, that no keyword MODEL is in file - in that case
    # line starts with ATOM keyword. This program saves data as
    # model with name 1
    while line.startswith('MODEL') or line.startswith('ATOM'):
        if line.startswith('MODEL'):
            model_name = int(line[11:14])
            line = file.readline()
        else:
            model_name = 1

        chains = []
        while line.startswith('ATOM'):
            chain, line = get_chain(file, line)
            chains.append(chain)

        model = Model(model_name, chains)
        models.append(model)

        # line is now ENDMDL
        line = file.readline()  # next MODEL or end of structural part
    return models


def read_pdb_file(path_to_file: str) -> Protein:
    with open(path_to_file) as file:
        models = get_models(file)
        protein = Protein(models)

        print(protein.models[-1].chains[-1].residues[-1].atoms[-1].name)
        return protein


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Input file in .sdf or .pdb format')
    args = parser.parse_args()

    input_file = args.file
    read_pdb_file(input_file)