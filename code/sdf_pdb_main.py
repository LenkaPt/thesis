import argparse
import sys
import os
from sdf import read_sdf_file
from pdb import read_pdb_file
from pdb_to_sdf import V2000, V3000


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Input file in .sdf or .pdb format')
    parser.add_argument('-c', '--conversion', dest='conversion_file',
                        help='Allows to convert protein from input pdb file to output sdf file. '
                             'Please specify name of output .sdf file.',
                        default=False)
    args = parser.parse_args()

    input_file = args.file
    try:
        if input_file.endswith('.sdf'):
            read_sdf_file(input_file)
        if input_file.endswith('.pdb'):
            protein = read_pdb_file(input_file)
    except ValueError as e:
        print(e)
        sys.exit(1)

    if args.conversion_file:
        with open(args.conversion_file, mode='w', encoding='utf8') as output_file:
            if protein[0].get_atom_count() <= 999:
                data = V2000(protein, output_file)
            else:
                # TODO - formát V3000 špatně - nelze zoobrazit pomocí Pymolu
                data = V3000(protein, output_file)
            data.write_to_file()
        os.startfile(args.conversion_file)
