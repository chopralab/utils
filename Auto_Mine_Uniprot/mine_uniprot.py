# Import statements
import requests
import argparse
import os
import pandas as pd


'''
Uniprot_Miner class used to mine FASTA information from the UniProt database

Static Methods:
    get_fasta_from_uniprot(id): Gets FASTA information for a specific protein based on provided UniProt ID
    get_all_fasta_from_list(fname): Gets the FASTA information from all of the proteins UniProt IDs provided in the input file
'''
class Uniprot_Miner():
    
    def __init__(self):
        pass

    '''
    Gets FASTA information for a specific protein based on provided UniProt ID

    Arguments:
        id: the UniProt ID of the protein

    Returns:
        The FASTA file text corresponding to the protein
    '''
    @staticmethod
    def get_fasta_from_uniprot(id):
        url_base = 'https://www.uniprot.org/uniprot/'
        url = url_base + str(id) + '.fasta'
        response = requests.get(url)
        if response:
            return response.text
        else:
            if response.status_code >= 500:
                print('Status Code ' + str(response.status_code)  + ': Server Error')
            elif response.status_code >= 400:
                print('Status Code' + str(response.status_code) + ': Client Error')
            return None
    

    '''
    Gets the FASTA information from all of the proteins UniProt IDs provided in the input file

    Arguments:
        fname: The path to the input file containing a list of the UniProt IDs
    '''
    @staticmethod
    def get_all_fasta_from_list(fname):
        if not os.path.isfile(fname):
            print("File Not Found: " + fname)
            return None

        df = pd.DataFrame(columns=['uniprot_id','fasta'])
    
        if fname.split('csv')[1] == 'txt':
            with open(fname, 'r') as file:
                for line in file:
                    fasta = Uniprot_Miner.get_fasta_from_uniprot(line)
                    if fasta.startswith('>'):
                        fasta = fasta.split('\n',1)[1]
                        fasta = fasta.replace('\n','')

                    df.loc[len(df)] = [line,fasta]

        elif fname.split('.')[1] == 'csv':
            csv_df = pd.read_csv(fname)
            uniprot_list = csv_df['UniprotID'].tolist()
            for elem in uniprot_list:
                fasta = Uniprot_Miner.get_fasta_from_uniprot(elem)
                if fasta.startswith('>'):
                    fasta = fasta.split('\n',1)[1]
                    fasta = fasta.replace('\n','')
                
                df.loc[len(df)] = [elem,fasta]

        return df


'''
Main method for running the module as a script

Command Line Arguments:
    --infile: 
    --outfile:

Output:
    CSV file with the following columns [uniprot_id, fasta]
'''
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, required=True)
    parser.add_argument('--outfile', type=str, default='outfile.csv')
    args = parser.parse_args()
        
    df = Uniprot_Miner.get_all_fasta_from_list(args.infile)
    df.to_csv(args.outfile, sep=',', index=False)

# Run the main method if the script is called
if __name__ == '__main__':
    main()             