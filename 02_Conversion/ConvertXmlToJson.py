###############################################################################
# Author:      Sebastian Peterka
# Company:     -
# File:        ConvertXmlToJson.py
# created:     07.12.2023 11:36
# edited:      -
#
# Description: The script extracts the attributes "accession," "name"
#              "iupac_name," "chemical_formula," "inchikey," and "smiles"
#              from <metabolite> elements in an XML file using ElementTree's
#              iterparse, ensuring memory efficiency, and converts the
#              extracted data into a JSON file.
#
# Comments:    Used Datasets:
#              * HMDB ->  All Metabolites (in XML format)
#                              Last updated: 2021-10-13
#                              Size: 6,04 GB
#
# Usage:       1. Download the mentioned datasets from the resource pages
#              2. Extract them
#              3. Move them to the data directory
#              4. Run the script
#
# Todo:        * add the path to the .xml file per user input
#              * add information about date/ time in the header of the file
#
# Resources:   https://hmdb.ca/downloads
#
###############################################################################

import xml.etree.ElementTree as ET
import json


def parse_metabolite(metabolite_element):
    key_order = ['accession', 'iupac_name',
                 'chemical_formula', 'inchikey', 'smiles']

    metabolite_data = {}

    for key in key_order:
        element = metabolite_element.find(f".//hmdb:{key}",
                                          namespaces={'hmdb':
                                                      'http://www.hmdb.ca'})
        if element is not None:
            metabolite_data[key] = element.text

    return metabolite_data


def save_as_json(data, json_file):
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    xml_file = "data/hmdb_metabolites.xml"
    json_file = "metabolites.json"

    metabolite_data = []

    for event, elem in ET.iterparse(xml_file, events=('start', 'end')):
        if elem.tag.endswith('metabolite') and event == 'end':
            metabolite_data.append(parse_metabolite(elem))
            elem.clear()  # Clear the element to free memory

    save_as_json(metabolite_data, json_file)

    print(f"Conversion completed. JSON file saved as {json_file}")
