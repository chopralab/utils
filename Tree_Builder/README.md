# Tree Builder Utility

## Overview

Creates a tree builder class that can submit queries to a deployed instance
of ASKCOS (CASP tool) to build a retrosynthetic route for given smiles strings.
    
## Required Packages / Scripts

## Usage

Example usage:
```python
import requests
from pprint import pprint

HOST = 'http://XX.XX.XX.XX' # Replace with your address with ASKCOS

# Create the tree builder
tb = TreeBuilder(HOST)
smi = 'CC1=C(COC2=CC(OC)=C(CN3CCCC[C@H]3C(O)=O)C(OC)=C2)C=CC=C1C4=CC=CC=C4'

# Parameters for the query
params = {
    'max_depth': 9,
    'max_branching': 25,
    'expansion_time': 60,
    'max_ppg': 100,
    'template_count': 1000,
    'max_cum_prob': 0.9999,
    'chemical_property_logic': 'none',
    'max_chemprop_c': 0,
    'max_chemprop_n': 0,
    'max_chemprop_o': 0,
    'max_chemprop_h': 0,
    'chemical_popularity_logic': 'none',
    'min_chempop_reactants': 5,
    'min_chempop_products': 5,
    'filter_threshold': 0.0001,
    'return_first': 'true' # default is false
}

# Build the tree for the given smiles and parameters
js = tb.build_tree(smi, params)
# Print the tree as a json
pprint(js)
```

## Testing
Tested on Python 3.6