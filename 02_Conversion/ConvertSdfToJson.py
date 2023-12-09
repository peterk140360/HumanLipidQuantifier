###############################################################################
# Author:      Sebastian Peterka
# Company:     -
# File:        ConvertSdfToJson.py
# created:     07.12.2023 11:22
# edited:      -
#
# Description: This Python program, using the RDKit library, converts (.sdf)
#              structural information from a Chemical Structure DataBase
#              file containing lipid data into a structured JSON format,
#              retaining key properties such as LM_ID, NAME, SYSTEMATIC_NAME,
#              CATEGORY, FORMULA, INCHI_KEY, and SMILES.
#
# Comments:    Used Datasets:
#              * LIPID MAPS -> LMSD Structure-data file (SDF)
#                              Last updated: 2023-12-06
#                              Size on disk: 267 MB
#              Run "pip install rdkit"
#
# Usage:       1. Download the mentioned datasets from the resource pages
#              2. Extract them
#              3. Move them to the data directory
#              4. Run the script
#
# Todo:        * add the path to the .sdf file per user input
#              * don't show the warnings in the terminal caused by rdkit
#              * add information about date/ time in the header of the file
#
# Resources:   https://lipidmaps.org/databases/lmsd/download
#              https://www.rdkit.org/docs/Install.html
#
###############################################################################

from rdkit import Chem
import json


def sdf_to_json(sdf_file, json_file):
    # Read the .sdf file
    suppl = Chem.SDMolSupplier(sdf_file)

    data = []

    # Iterate through the molecules in the .sdf file
    for mol in suppl:
        if mol:
            mol_data = {}
            mol_data['LM_ID'] = (
                mol.GetProp('LM_ID')
                if mol.HasProp("LM_ID") else None)
            mol_data['NAME'] = (
                mol.GetProp('NAME')
                if mol.HasProp("NAME") else None)
            mol_data['SYSTEMATIC_NAME'] = (
                mol.GetProp('SYSTEMATIC_NAME')
                if mol.HasProp("SYSTEMATIC_NAME") else None)
            mol_data['CATEGORY'] = (
                mol.GetProp('CATEGORY')
                if mol.HasProp("CATEGORY") else None)
            mol_data['FORMULA'] = (
                mol.GetProp('FORMULA')
                if mol.HasProp("FORMULA") else None)
            mol_data['INCHI_KEY'] = (
                mol.GetProp('INCHI_KEY')
                if mol.HasProp("INCHI_KEY") else None)
            mol_data['SMILES'] = (
                mol.GetProp('SMILES')
                if mol.HasProp("SMILES") else None)
            data.append(mol_data)

    # Save the data as a JSON file
    with open(json_file, 'w') as outfile:
        json.dump(data, outfile, indent=2)


if __name__ == "__main__":
    sdf_file = "data/structures.sdf"
    json_file = "lipids.json"

    sdf_to_json(sdf_file, json_file)

    print(f"Conversion completed. JSON file saved as {json_file}")
