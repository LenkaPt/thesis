import argparse
import sys
from sdf import read_sdf_file
from pdb import read_pdb_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Input file in .sdf or .pdb format')
    args = parser.parse_args()

    input_file = args.file
    try:
        if input_file.endswith('.sdf'):
            read_sdf_file(input_file)
        elif input_file.endswith('.pdb'):
            read_pdb_file(input_file)
    except ValueError as e:
        print(e)
        sys.exit(1)
