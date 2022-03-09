import subprocess
import json


def available_methods():
    available_methods = subprocess.run(['curl', 'http://page'],
                                       capture_output=True, text=True)
    return json.loads(available_methods.stdout)


def success(json):
    return json['code'] == 200:


def available_parameters(method):
    available_parameters = subprocess.run(['curl', 'http://page?method=' + method],
                                          capture_output=True, text=True)

    return json.loads(available_parameters.stdout)


def trying_available_parameters_route(method):
    print(f'Trying available_parameters route\n')
    data = available_parameters(method)
    if success(data):
        print(f'Method: {method}\n{data["data"]}\n')
    else:
        print(f'Method: {method}\n{data["message"]}\n')
    print(f'------------\n')


def send_file(file):
    identifiers = subprocess.run(['curl', '-X', 'POST', '-F',
                                  'file[]=@' + file,
                                  'http://page'],
                                 capture_output=True, text=True)
    return json.loads(identifiers.stdout)


def trying_send_file_route(file):
    print(f'Trying send_file route ({file})\n')
    data = send_file(file)
    if success(data):
        print(f'Identifiers: {data["data"]}')
    else:
        print(f'Error: {data["message"]}')
    print(f'----------------\n')


def suitable_methods(identifier):
    suitable_methods = subprocess.run(['curl', 'http://page?identifier[]=' + identifier],
                                          capture_output=True, text=True)
    return json.loads(suitable_methods.stdout)


def trying_suitable_method_route(identifiers):
    for identifier in identifiers:
        print(f'Trying suitable_methods route ({identifier})\n')
        data = suitable_methods(identifier)
        if success(data):
            print(f'Suitable methods: {data["data"]}')
        else:
            print(f'Error: {data["message"]}')
        print(f'----------------\n')


def calculate_charges(identifier):
    calculated_charges = subprocess.run(['curl', 'http://page?identifier=' + identifier],
                                      capture_output=True, text=True)
    return json.loads(calculated_charges.stdout)


def trying_calculate_charges_route(identifiers):
    for identifier in identifiers:
        print(f'Trying calculate_charges route ({identifier})\n')
        data = calculate_charges(identifier)
        if success(data):
            print(f'Calculated charges: {data["data"]}')
        else:
            print(f'Error: {data["message"]}')
        print(f'----------------\n')
