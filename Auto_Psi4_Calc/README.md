# Automated Quantum Mechanical (QM) Calculations Using Psi4

## Overview

This scipt is used to automate a range of Psi4 QM calculations for a set of molecules and provide a standardized output file with calculation information.

Portions of the code use [this webpage](https://iwatobipen.wordpress.com/2018/08/27/calculate-homo-and-lumo-with-psi4-reviced-rdkit-psi4/
) as a reference.

## Required Packages
[`pandas`](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)

[`psi4`](https://anaconda.org/psi4/psi4)

[`rdkit`](https://www.rdkit.org/docs/Install.html)

All other packages used are built into Python 3.X

## Input File Format

The input file should be structure with one SMILES string per line and no additional information. See the example format below:

```txt
smiles1
smiles2
smiles3
...
```

## Command Line Arguments

### `--infile`

Name/path of the input file.

Default value is `infile.txt`

### `--outfile`

Name/path of the CSV output file. This file contains formatted information on energy/homo/lumo information.

Default value is `outfile.csv`

### `--corefile`

Name/path of the Psi4 core output file. This file contains detailed information on the calculations run, include geometry, vibrational frequency, etc.

Default value is `output.dat`

### `--energy`

Flag for running single-point electronic energy calculations.

### `--homolumo`

Flag for running single-point electronic energy calculations and then accessing homo and lumo orbital values.

### `--geomopt`

Flag for running geometry optimization.

### `--vibfreq`

Flag for running vibrational frequency analysis

### `--method`

Name of the method that will be used in the calculations. For an example of available methods for energy see [here](https://psicode.org/psi4manual/1.3.2/energy.html)

Default value is `HF`

### `--basis`

Name of the basis set that will be used in the QM calculations. Basis sets supported by Psi4 can be found [here](https://psicode.org/psi4manual/master/basissets_tables.html#apdx-basistables).

Default value is `6-31G*`

## Output File Format

The CSV output file has a header and then smiles/energy/homo/lumo information on each line, corresponding to the calculation on the molecule for that line. See the example format below:

```
smiles,energy,homo,lumo
smiles1,energy1,homo1,lumo1
smiles2,energy2,homo2,lumo2
smiles3,energy3,homo3,lumo3
...
```

Remeber that there is also the `.dat` file that corresponds to the Psi4 core output.

## Testing

Tested on Python 3.6.8
