import requests
from concurrent.futures import ThreadPoolExecutor
import concurrent
import pytest


def available_methods():
    return requests.get('http://78.128.250.156/available_methods')


def test_available_methods():
    response = available_methods().json()
    assert 'eem' in response['available_methods']


def available_parameters(method):
    return requests.get('http://78.128.250.156/available_parameters',
                        params={'method': method})


@pytest.mark.parametrize('method, expected', [
    ('eem', 'OK'),
    ('aam', 'Method aam is not available.')
])
def test_available_parameters(method, expected):
    response = available_parameters(method).json()
    assert response['message'] == expected


def send_file(file):
    return requests.post('http://78.128.250.156/send_files',
                         files={'file[]': open(file)})


valid_file_identifier = send_file('./../../1ner.pdb').json()['structure_ids']['1ner']
invalid_file_identifier = send_file('./../../code/6a5j.pdb').json()['structure_ids']['6a5j']


@pytest.mark.parametrize('structure_id, expected', [
    ('./../../1ner.pdb', 'OK'),
    ('./../../test.txt', 'Unsupported format'),
    ('./../../4wfb.pdb', 'exceeds the capacity limit')
])
def test_send_file(structure_id, expected):
    response = send_file(structure_id).json()
    assert expected in response['message']


def suitable_methods(structure_id):
    return requests.get('http://78.128.250.156/suitable_methods',
                        params={'structure_id': [structure_id]})


@pytest.mark.parametrize('structure_id, expected', [
    (valid_file_identifier, 'OK'),
    (invalid_file_identifier, 'No molecules were loaded'),
    (None, 'not specified structure ID'),
    ('jhskhk', 'Structure ID jhskhk does not exist')
])
def test_suitable_methods(structure_id, expected):
    response = suitable_methods(structure_id).json()
    assert expected in response['message']


def calculate_charges(structure_id, method, parameters):
    return requests.get('http://78.128.250.156/calculate_charges',
                        params={'structure_id': structure_id,
                                'method': method,
                                'parameters': parameters})

@pytest.mark.parametrize('structure_id, method, parameters, expected', [
    (valid_file_identifier, 'denr', 'DENR_00_from_QEq', 'OK'),
    (valid_file_identifier, None, None, 'OK'),
    (valid_file_identifier, 'eqeq', 'gdfnxrg', 'Parameters gdfnxrg are not available for method eqeq'),
    (valid_file_identifier, 'eqeq', 'DENR_00_from_QEq', 'Parameters DENR_00_from_QEq are not available for method eqeq'),
    (valid_file_identifier, 'dgjybgcakgr', 'DENR_00_from_QEq', 'Method dgjybgcakgr is not available'),
    (valid_file_identifier, 'eem', None, 'Method eem is not suitable'),
    (valid_file_identifier, 'denr', 'EEM_00_NEEMP_ccd2016_npa', 'Parameters EEM_00_NEEMP_ccd2016_npa are not available for method denr'),
    (invalid_file_identifier, 'eem', None, 'No molecules were loaded'),
    (None, None, None, 'not specified structure ID'),
    ('jhskhk', None, None, 'Structure ID jhskhk does not exist')
])
def test_calculate_charges(structure_id, method, parameters, expected):
    response = calculate_charges(structure_id, method, parameters).json()
    assert expected in response['message']



def pdb_id(identifier):
    return requests.post('http://78.128.250.156/pdb_id', params={'pid[]': identifier})


@pytest.mark.parametrize('structure_id, expected', [
    ('1ner', 'OK'),
    ('4wfb', 'bigger than 10 Mb'),
    ('hgkcgyargy', 'Not Found for url'),
    (None, 'No pdb id specified')
])
def test_pdb_id(structure_id, expected):
    response = pdb_id(structure_id).json()
    assert expected in response['message']


def cid(identifier):
    return requests.post('http://78.128.250.156/pubchem_cid', params={'cid[]': identifier})

sdf_identifier = cid('1').json()['structure_ids']['1']

@pytest.mark.parametrize('structure_id, expected', [
    ('1', 'OK'),
    ('hgkcgyargy', 'BadRequest'),
    (None, 'No pubchem cid specified')
])
def test_cid(structure_id, expected):
    response = cid(structure_id).json()
    assert expected in response['message']


def add_hydrogens(identifier):
    return requests.post('http://78.128.250.156/add_hydrogens', params={'structure_id': identifier})


@pytest.mark.parametrize('structure_id, expected', [
    (valid_file_identifier, 'OK'),
    (sdf_identifier, 'not in .pdb or .cif format'),
    (invalid_file_identifier, 'Error occurred when using pdb2pqr30')
])
def test_add_hydrogens(structure_id, expected):
    response = add_hydrogens(structure_id).json()
    assert expected in response['message']






#######################################################################################################
identifiers = ['12as7e10b014', '12asdfjsaofa', '12asseh_hs9c', '12asct92u8an', '12as9cotg6v2']
def concurrent_tasks(identifiers):
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(suitable_methods, identifier) for identifier in identifiers}

        for future in concurrent.futures.as_completed(futures):
            print(f'{future.result().json()["data"]}')
            print('--------------')

# concurrent_tasks(identifiers)
