# Automated Mining of PubChem Database using RESTful API (PUG REST)

## Overview

This scirpt is used to automate mining of the [PubChem Database](https://pubchem.ncbi.nlm.nih.gov/) using their implement [PUG REST RESTful API](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest). A RESTful API is a way to directly access formatted database information using a simple URL call. Visit [PUG REST tutorial](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial) for more information on how to use the PUG REST API.

## Required Packages

[`pandas`](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)

All other packages used are built into Python 3.X

## Input File Format

The input file should contain one PubChem recognized identifier (see **Command Line Arguments** for details) per line. This will identify the compound of interest when accessing the database. **Note that all of the identifiers must be of the same type for a single input file**. An example input file for the compound name identifier is shown below.

```txt
glucose
sucrose
fructose
...
```

## Command Line Arguments

This script uses command line arguments to specify information on the `input`, `data`, and `output` of the PUG REST formatted URL call. For more information on how PUG REST URL's are formatted, see the [following page](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial) and go to the section labeled `Input: Design of the URL`.

### General Arguments

These arguments are used to provide general information when running the script.

### `--infile`

Name/path of the input file.

Default value is `infile.txt`

### `--outfile`

Name/path of the output file.

Default value is `outfile.txt`

### `--no-append`

Flag for appending output. If Flagged (`True`) then output for each URL call will be placed into seperate files. If unflagged (`Flase`) then all output will be appended into a single file.

Default value is `False`

### `--url`

Flag for input file containing a list of pre formatted PUG REST URL's.

Default values is `False`

### `--domain-override`

Flag for overriding the domain of access in the PubChem database. Do not use this flag unless your are sure the domain of access needs to be changed. You can read more on PubChem domains [here](https://pubchemblog.ncbi.nlm.nih.gov/2014/06/19/what-is-the-difference-between-a-substance-and-a-compound-in-pubchem/).

Choices are `substance`, `compound`, or`assay`

### Input Arguments

These arguments are used to provide information on the input used to identify compound(s), substances(s), or assay(s) of interest in the PubChem database.

### `--in-cid`

Flag to `True` for input being PubChem compound ID's.

Default value is `False`

### `--in-sid`

Flag to `True` for input being PubChem substance ID's.

Default value is `False`

### `--in-name`

Flag to `True` for input being compound names.

Default value is `False`

### `in-smiles`

Flag to `True` for input being SMILES strings.

Default value is `False`

### `--in-formula`

Flag to `True` for input being chemical formulas.

Default value is `False`

### `--in-inchikey`

Flag to `True` for input being InChi keys.

Default value is `False`

### `--in-xref`

Flag to `True` for input being information from a cross referenced database.

Default value is `False`

### `--in-xtype`

Must be used with `--in-xref`. This argument specifies the type of cross referenced information provided as input.

Choices are `RegistryID`, `RN`, `PubMedID`, `MMDBID`, `ProteinGI`, `NucleotideGI`, `TaxonomyID`, `MIMID`, `GeneID`, `ProbeID`, and `PatentID`

### `--in-structure`

Must be used with `--in-smiles` or `--in-cid`. This argument specifies any type of modified structure search (such as substructure).

Choices are `substructure`, `superstructure`, `similarity`, and `identity`

### `--direct-input`

Argument for provided the input segment of the URL directly.

Example usage would be `--direct-input=compound/name/glucose/`

### Data Arguments

These arguments are used to provide information on the data accessed from the PubChem database.

### `--data-record`

Flag to `True` for accessing PubChem full data record.

Default value is `False`

### `--data-summary`

Flag to `True` for accessing PubChem compound summary.

Default value is `Flase`

### `--data-property`

Argument for accessing a physical or chemical property (or properties) of a compound

Example usage would be `--data-property=MolecularWeight,MolecularFormulat,InChi`

### `--data-synonyms`

Flag to `True` for accessing synonyms of a compound

Default value is `False`

### `--data-xref`

Flag to `True` for accessing cross referenced data for a compound

Default values is `False`

### `--data-xtype`

Must be used with `--data-xref`. This argument specifies the type of cross referenced information provided as input.

Choices are `RegistryID`, `RN`, `PubMedID`, `MMDBID`, `ProteinGI`, `NucleotideGI`, `TaxonomyID`, `MIMID`, `GeneID`, `ProbeID`, and `PatentID`

### BioAssay Arguments

These arguments are used to access the BioAssay portion of the PubChem database.

### `--assay`

Flag to `True` for use of the `assay` domain and other BioAssay flags.

Default value is `False`

### `--aid`

Flag to `True` for use of BioAssay ID as input.

Default value is `False`

### `--atarget`

Argument for providing specific PubChem BioAssay target(s) to mine from the database.

Example usage is `--atarget=ProteinGI,ProteinName,GeneID`

### `--adescription`

Flag to `True` to access description of selected BioAssays.

Default value is `False`

### `--asummary`

Flag to `True` to access summary of selected BioAssays.

Default values is `False`

### Output Arguments

These arguments are used to provide information on the output format of the data accessed from the PubChem database.

### `--outtype`

Specify the type of output used (reccomended it matches up with the suffix of the outfile).

Choices are `XML`, `JSON`, `JSONP`, `ASNB`, `ASNT`, `SDF`, `CSV`, `PNG`, and `TXT`

### `--callback`

Used in conjuection with JSONP outtype. You can read more about JSONP articles [here](https://en.wikipedia.org/wiki/JSONP).

## Output File Format

The output file format varies based on the specific output file type and data selected. Please see the **Examples** section below for example output file formats.

## Examples
