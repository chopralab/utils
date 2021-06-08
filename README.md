# ChopraLab Computational Utilities

## Overview

Welcome to the Chopra Lab utilities repository. This repo contains information on commonly used data-parsing, cheminformatic, bioinformatic, machine learning, and other utility scripts. The scripts are generalized and all runtime parameters are specified as runtime arguments unless stated otherwise. Check the `List of Utilities` section below to see what utilities are available and check the `Submission Instructions` for instruction on how to submit your own utility to the repo.

## List of Utilities

### `Auto_Mine_PubChem` - Automated Mining of PubChem Database using RESTful API (PUG REST)

### `Auto_Psi4_Energy_Cal` - Automated Quantum Mechanical (QM) Calculations Using Psi4

### `Auto_Mine_Uniprot` - Automated Mining of FASTA Information from the UniProt Database

## Submission Instructions

1. Clone the repository to your own personal machine
2. If you plan for the utility to be a work in progress, it is suggested that you create a new branch
    - `git checkout -b [branch_name]`
3. Create a `new folder` in the `base directory`
4. Add your scripts to the `new folder`
5. Create a `README` in the `new folder` and describe your scripts and how to use them
6. Create a `.gitignore` in the `new folder` if you want to ignore any files
    - Please use this to ignore any testing files/scripts
    - **Do Not create a global `.gitignore` file as this can affect other files in other script folders**
7. Update the `README` in the `base directory` by adding the appropriate information under the `List of Utilities` section
8. Regardless of the branch you are working on (`master` or `[branch_name]`) you will need to run the following:
    - `git add .`
    - `git commit`
    - `git push origin [branch_name]`
9. If you are working on a seperate branch you will need to run the above commands in addition to the following to merge to master:
    - `git checkout master`
    - `git pull origin master`
    - `git merge [branch_name]`
    - Resolve any conflicts present
    - `git push origin master`
