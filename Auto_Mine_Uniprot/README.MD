# Automated Mining of FASTA Sequence Information from the UniProt Database

## Overview

This module and script is used to mine FASTA file sequence information given a UniProt ID. Individual sequences can be mined from the database using the Uniprot_Miner class or the script can be called with command line arguments to mine them from a file containing a list of UniProt IDs.

## Required Packages

[`pandas`](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)

## Input File Format

The input file should be structured with one UniProt ID per line and no additional information. See the example format below:

```txt
UniProtID1
UniProtID2
UniProtID3
...
```

## Command Line Arguments

### `--infile`

Name/path of the input file.

This argument is **Required**

### `--outfile`

Name/path of the output file.

Default is `outfile.csv`

## Output File Format

The output of the script will be a CSV file with two columns, the UniProt ID and the FASTA sequence. The file will have a header, see the example format below:

```csv
uniprot_id,fasta
UniProtID1,FASTA1
UniProtID2,FASTA2
UniProtID3,FASTA3
...
```

## Module Usage

To use the `mine_uniprot.py` script as a module with access to the `Uniprot_Miner` class you must import and use it into your script in the following manner:

```py
import sys
sys.path.insert(0, path_to_mine_uniprot_directory)
from mine_uniprot import Uniprot_Miner

print(Uniprot_Miner.get_fasta_from_uniprot(my_uniprot_id))
...
```

Both methods are static methods so you must call them from a class object or import them directly.

## Testing

Tested on Python 3.7.1
