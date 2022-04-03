from flask import Flask, render_template, request, send_from_directory, jsonify
from sdf import read_sdf_file
from pdb import read_pdb_file
from conversion import convert_pdb_to_sdf
import tempfile
import os
import csv
import json
import chargefw2_python
import requests
from time import sleep
import subprocess


app = Flask(__name__)


@app.route('/')
def how_to_use():
    return render_template('index.html')


def json_ok(data):
    return jsonify({"status_code": 200,
                    "message": "OK",
                    "data": data})


def json_error(message, status_code=404):
    return jsonify({"status_code": status_code,
                    "message": message})


@app.route('/available_methods')
def available_methods():
    return json_ok(chargefw2_python.get_available_methods())


@app.route('/available_parameters')
def available_parameters():
    method = request.args.get('method')
    if not method:
        return json_error(f'You have not specified method. '
                          f'Add to URL following, please: ?method=your_chosen_method')

    if method not in chargefw2_python.get_available_methods():
        return json_error(f'Method {method} is not available.')

    return json_ok(chargefw2_python.get_available_parameters(method))


def valid_suffix(files):
    for file in files:
        if not file.filename.endswith(('sdf', 'pdb', 'mol2', 'cif')):
            return False
    return True


def save_file_identifiers(identifiers):
    # output_file_path = "/tmp/file_identifiers.json"
    # with open(output_file_path, mode="w+") as file:
    #     if os.stat(output_file_path).st_size == 0:
    #         file_data = []
    #     else:
    #         file_data = json.load(file)
    #
    # with open(output_file_path, mode="w") as file:
    #     for identifier, path_to_file in identifiers.items():
    #         file_data.append({"identifier": identifier,
    #                           "path_to_file": path_to_file})
    #     json.dump(file_data, file, indent=4)
    hidden_file = '/tmp/.uploading'
    while os.path.exists(hidden_file):
        sleep(0.01)
    open_hidden_file = open(hidden_file, 'w')
    with open("/tmp/file_identifiers.csv", mode="a") as file:
        writer = csv.writer(file)
        for identifier, path_to_file in identifiers.items():
            writer.writerow([identifier, path_to_file])

    open_hidden_file.close()
    os.remove(hidden_file)


@app.route('/send_files', methods=['GET', 'POST'])
def send_files():
    if request.method == 'POST':
        files = request.files.getlist('file[]')

        if not files:
            return json_error(f'No file sended. Add to URL following, please: ?file[]=path_to_file')

        if not valid_suffix(files):
            return json_error(f'Unsupported file. Send only .sdf and .pdb files.')

        tmpdir = tempfile.mkdtemp()
        identifiers = {}
        for file in files:
            path_to_file = os.path.join(tmpdir, file.filename)
            file.save(path_to_file)

            identifier = os.path.basename(tempfile.NamedTemporaryFile(prefix=file.filename.rsplit('.')[0]).name)
            identifiers[identifier] = path_to_file

        save_file_identifiers(identifiers)

        return json_ok(list(identifiers.keys()))

    return '''
                 <form method=post enctype=multipart/form-data>
                 <input type=file name=file[] multiple=true>
                 <input type=submit value=Upload>
                 </form>'''


@app.route('/pdb_id', methods=['GET', 'POST'])
def pdb_id():
    pdb_identifiers = request.args.getlist('pid[]')

    if not pdb_identifiers:
        return json_error('No pdb id specified. Add to URL following, please: ?pid[]=pdb_id')

    tmpdir = tempfile.mkdtemp()
    identifiers = {}
    for pdb_id in pdb_identifiers:
        r = requests.get('https://files.rcsb.org/download/' + pdb_id + '.pdb', stream=True)
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            return json_error(f'{e}', r.status_code)

        # save requested pdb structure
        path_to_file = os.path.join(tmpdir, pdb_id + '.pdb')
        with open(path_to_file, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

        identifier = os.path.basename(tempfile.NamedTemporaryFile(prefix=pdb_id).name)
        identifiers[identifier] = path_to_file

    save_file_identifiers(identifiers)

    return json_ok(list(identifiers.keys()))


def get_path_based_on_identifier(identifier):
    # with open("/tmp/file_identifiers.json", mode="r") as file:
    #     data = json.load(file)
    #     for record in data:
    #         if record["identifier"] == identifier:
    #             return record["path_to_file"]
    #     return None
    with open("/tmp/file_identifiers.csv", mode="r") as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] == identifier:
                return row[1]
    return None


def get_molecules(identifier, read_hetatm = True, ignore_water = False):
    path_to_file = get_path_based_on_identifier(identifier)
    if not path_to_file:
        raise ValueError(f'Identifier {identifier} does not exist.')
    try:
        return chargefw2_python.Molecules(path_to_file, read_hetatm, ignore_water)
    except RuntimeError as e:
        raise ValueError(e)


@app.route('/add_hydrogens')
def add_hydrogens():
    identifier = request.args.get('identifier')
    if not identifier:
        return json_error(f'You have not specified identifier obtained after uploading your file. '
                          f'Add to URL following, please: ?identifier=obtained_identifier')

    pH = request.args.get('pH')
    if not pH:
        pH = 7.0

    input_file = get_path_based_on_identifier(identifier)
    output_file = os.path.join(tempfile.mkdtemp(), 'result.pqr')

    # hydrogen bond optimalization
    noopt = request.args.get('noopt')

    if not noopt:
        subprocess.run(['pdb2pqr30', f'--noopt', f'--pH', f'{pH}', f'{input_file}', f'{output_file}'])
    else:
        subprocess.run(['pdb2pqr30', f'--pH', f'{pH}', f'{input_file}', f'{output_file}'])

    output_identifier = os.path.basename(tempfile.NamedTemporaryFile(prefix='hydro').name)
    save_file_identifiers({output_identifier: output_file})
    return json_ok(output_identifier)



def get_suitable_methods(identifier, read_hetatm, ignore_water):
    try:
        molecules = get_molecules(identifier, read_hetatm, ignore_water)
    except ValueError as e:
        raise ValueError(e)
    return chargefw2_python.get_suitable_methods(molecules)


@app.route('/suitable_methods',  methods=['GET', 'POST'])
def suitable_methods():
    identifiers = request.args.getlist('identifier[]')
    read_hetatm = request.args.get('read_hetatm')
    ignore_water = request.args.get('ignore_water')
    if not identifiers:
        return json_error(f'You have not specified identifier obtained after uploading your file. '
                          f'Add to URL following, please: ?identifier[]=obtained_identifier')

    read_hetatm = read_hetatm is None
    ignore_water = ignore_water is None
    data = {}
    for identifier in identifiers:
        try:
            suitable_methods = get_suitable_methods(identifier, read_hetatm, ignore_water)
        except ValueError as e:
            return json_error(str(e))
        data[identifier] = suitable_methods

    return json_ok(data)


def get_calculated_charges(identifier, method, parameter, read_hetatm, ignore_water):
    try:
        molecules = get_molecules(identifier, read_hetatm, ignore_water)
    except ValueError as e:
        raise ValueError(e)
    return chargefw2_python.calculate_charges(molecules, method, parameter)


def get_txt(charges, identifier):
    with tempfile.TemporaryDirectory() as tmpdir:
        filename = os.path.basename(get_path_based_on_identifier(identifier))
        filename = filename + '.txt'
        path_to_output_file = os.path.join(tmpdir, filename)
        with open(path_to_output_file, mode="w") as output_file:
            for molecule in charges:
                output_file.write(molecule + '\n')
                output_file.write(" ".join([str(charge) for charge in charges[molecule]]))
                output_file.write('\n')
        return send_from_directory(tmpdir, filename, as_attachment=True)


@app.route('/calculate_charges',  methods=['GET', 'POST'])
def calculate_charges():
    identifier = request.args.get('identifier')
    method = request.args.get('method')
    parameters = request.args.get('parameters')

    read_hetatm = request.args.get('read_hetatm')
    ignore_water = request.args.get('ignore_water')

    read_hetatm = read_hetatm is None
    ignore_water = ignore_water is None

    if not identifier:
        return json_error(f'You have not specified identifier obtained after uploading your file. '
                          f'Add to URL following, please: identifier=obtained_identifier')

    if not method:
        try:
            suitable_methods = get_suitable_methods(identifier, read_hetatm, ignore_water)
        except ValueError as e:
            return json_error(str(e))
        method = list(suitable_methods)[1]
        parameters = None        # if method does not require parameters
        if suitable_methods[method]:
            parameters = suitable_methods[method][0]

    try:
        charges = get_calculated_charges(identifier, method, parameters, read_hetatm, ignore_water)
    except ValueError as e:
        return json_error(str(e))
    return json_ok(charges)


@app.route('/testing_long_answer',  methods=['GET', 'POST'])
def testing_long_answer():
    identifier = request.args.get('identifier')
    if not identifier:
        return json_error("Identifier missing.")
    sleep(300)
    return identifier



if __name__ == '__main__':
    app.run(host='0.0.0.0')
