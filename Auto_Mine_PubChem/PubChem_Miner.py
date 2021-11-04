# Import statements
from time import sleep
import requests
import pandas as pd
import json
import xml.etree.ElementTree as ET
from io import StringIO
from rdkit import Chem
from time import sleep

# Exception for invalid command line arguments
class ArgumentError(Exception):
    pass

# Exception for PubChem REST API error status codes
class StatusCodeError(Exception):

    def __init__(self, status_code):
        message = 'Unknown Status Code Received: ' + str(status_code)

        if status_code == 400:
            message = 'Status Code 400 - PUGREST.BadRequest - Request is improperly formed (syntax error in the URL, POST body, etc.)'
        elif status_code == 404:
            message = 'Status Code 404 - PUGREST.NotFound - The input record was not found (e.g. invalid CID)'
        elif status_code == 405:
            message = 'Status Code 405 - PUGREST.NotAllowed - Request not allowed (such as invalid MIME type in the HTTP Accept header)'
        elif status_code == 504:
            message = 'Status Code 504 - PUGREST.Timeout - The request timed out, from server overload or too broad a request'
        elif status_code == 503:
            message = 'Status Code 503 - PUGREST.ServerBusy - Too many requests or server is busy, retry later'
        elif status_code == 501:
            message = 'Status Code 501 - PUGREST.Unimplemented - The requested operation has not (yet) been implemented by the server'
        elif status_code == 500:
            message = 'Status Code 500 - PUGREST.ServerError - Some problem on the server side (such as a database server down, etc.)'
 
        super().__init__(message)

class PubChem_Miner():

    def __init__(self):
        pass

    @staticmethod
    def get_compound_info(input_type,input,data_type,data,output_type='TXT'):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        domain = "compound/"

        if input_type is 'smiles':
            input = 'smiles/' + input + '/'
        elif input_type is 'inchikey':
            input = 'inchikey/' + input + '/'
        elif input_type is 'name':
            input = 'name/' + input + '/'
        else:
            print("Input type not supported")
            return None

        if data_type is 'property':
            if isinstance(data, list):
                dcopy = 'property/'
                for elem in data:
                    if elem is 'smiles':
                        dcopy += 'canonicalSMILES' + ','
                    else:
                        dcopy += elem + ','
                dcopy = dcopy[:-1]
                data = dcopy + '/'
            else:
                if data is 'smiles':
                    data = 'property/' + 'canonicalSMILES' + '/'
                else:
                    data = 'property/' + data + '/'
        elif data_type is 'synonyms':
            data = 'synonyms' + '/'
        elif data_type is 'assay':
            if output_type not in ['XML','CSV','JSON']:
                print("Invalid output type for this query")
                return None
            else:
                data = 'assaysummary/' 
        else:
            print("Data type not supported")
            return None

        if output_type not in ['TXT','XML','CSV','JSON']:
            print('Output type not supported')
            return None

        url = url + domain + input + data + output_type
        response = requests.get(url)
        if response:
            print("Request on " + url + " succeded!")
            if output_type is 'TXT':
                return response.text
            if output_type is 'CSV':
                return pd.read_csv(StringIO(response.text))
            if output_type is 'XML':
                return ET.ElementTree(ET.fromstring(StringIO(response.text)))
            if output_type is 'JSON':
                return json.loads(StringIO(response.text))
        else:
            try:
                raise StatusCodeError(response.status_code)
            except StatusCodeError:
                print("Request on " + url + " failed")
                return None

    @staticmethod
    def get_json_info(smiles,delay=0.2,printing=False):
        try:
            mol = Chem.MolFromSmiles(smiles)
            inchi_key = Chem.inchi.MolToInchiKey(mol)
        except Exception as e:
            if printing:
                print("RDKit failure on " + smiles)
        
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"+str(inchi_key)+"/cids/TXT"

        response = requests.get(url)
        sleep(delay)
        if response:
            cid = response.text.split("\n")[0]
            if printing:
                print("Request on " + smiles + " Succeded! - CID Accessed: " + cid)
        else:
            if printing:
                print("Request on " + smiles + " Failed! - No CID Accessed")
            print(response.text)
            return smiles, None

        json_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/"+ str(cid) +"/JSON/?response_type=display"
        response = requests.get(json_url)
        sleep(delay)
        if response:
            if printing:
                print("Requeset on " + cid + " Succeded - JSON Accssed")
            data = json.loads(response.text)
            return smiles, data
        else:
            if printing:
                print("Request on " + cid + " Failed! - No JSON Accessed")
            if response == 503:
                print(response.headers)
                print(response.text)
            return smiles, None