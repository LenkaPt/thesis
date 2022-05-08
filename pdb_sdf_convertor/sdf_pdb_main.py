import argparse
import sys
from sdf import read_sdf_file
from pdb import read_pdb_file
from conversion import convert_pdb_to_sdf
from pathlib import Path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sdf', '--sdf_file', help='Input file in .sdf format')
    parser.add_argument('-pdb', '--pdb_file', help='Input file in .pdb format')
    parser.add_argument('-c', '--conversion', dest='conversion_file',
                        help='Allows to convert protein from input pdb file to '
                             'output sdf file. Please specify name of output .sdf '
                             'file and input .pdb file.',
                        default=False)
    parser.add_argument('-d', '--directory', help='Specify directory to '
                                                  'convert all .pdb files in directory.')
    args = parser.parse_args()

    if args.sdf_file:
        if not args.sdf_file.endswith('.sdf'):
            print('Please specify file in .sdf format.')
            sys.exit(1)
        try:
            read_sdf_file(args.sdf_file)
        except ValueError as e:
            print(e)
            sys.exit(1)

    if args.pdb_file and args.conversion_file:
        try:
            protein = read_pdb_file(args.pdb_file)
            convert_pdb_to_sdf(protein, args.conversion_file)
        except ValueError as e:
            print(e)
            sys.exit(1)

    elif args.pdb_file:
        if not args.pdb_file.endswith('.pdb'):
            print('Please specify file in .pdb format.')
            sys.exit(1)
        try:
            protein = read_pdb_file(args.pdb_file)
        except ValueError as e:
            print(e)
            sys.exit(1)

    if args.directory and args.conversion_file:
        my_path = Path(args.directory)
        # print(f'Absolute path: {my_path.absolute()}')
        # print(f'Absolute path os: {os.path.abspath(my_path)}')
        for file in sorted(my_path.rglob('*.pdb')):
            print(file.name)
            if file.name == '3onz.pdb':
                continue
            protein = read_pdb_file(file)

            convert_pdb_to_sdf(protein, args.conversion_file)
