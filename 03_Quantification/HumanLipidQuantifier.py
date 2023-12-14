###############################################################################
# Author:      Sebastian Peterka
# Company:     -
# File:        HumanLipidQuantifier.py
# created:     06.12.2023 10:32
# last edited: 14.12.2023
#
# Description: HumanLipidQuantifier is a specialized software designed for
#              accurately quantifying lipids within human cells.
#              The program employs LIPID MAPS as a reference database and
#              systematically analyzes input data sets from the
#              Human Metabolome Database (HMDB).
#
# Comments:    Used Datasets:
#              * LIPID MAPS -> LMSD Structure-data file (SDF)
#                              Last updated: 2023-12-06
#                              Size on disk: 267 MB
#              * HMDB ->       All Metabolites (XML)
#                              Last updated: 2021-10-13
#                              Size: 6,04 GB
#              * Formula 1: only change the values for X, r and v
#                X = size_common; r = radius_metabolite; v = radius_lipid
#              * Formula 2: only change values for d, r_1 and r_2
#                d = pos_common[0]; r_1 = radius_metabolite; r_2 = radius_lipid
#              * To check the results the geogebra-export.ggb file can be
#                imported by using the link provided in the ressources section.
#
# Usage:       1. Make sure the metabolites.json and the lipids.json files
#                 are located in the data-directory
#              2. Run the script by opening the cli and typing
#                 "python HumanLipidQuantifier.py"
#              3. Deside if you want to change save_fig; with_logo or save_txt
#
# Todo:        * solve the formula for the distance d and the intersection
#                point s in this script instead of solving in wolframalpha
#                Formulas as seen in formula.txt
#              * Add user input for save_fig; with_logo or save_txt
#
# Resources:   https://lipidmaps.org/databases/lmsd/download
#              https://hmdb.ca/downloads
#              https://www.arndt-bruenner.de/mathe/scripts/kreissehnen.htm
#              https://kratochwill.lima-city.de/PDFLinAlg2/K06.pdf
#              https://www.geogebra.org/geometry
#
###############################################################################

import os
import json
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math


def extract_values(data, key):
    """
    Extract unique values of a specified key from data entries and
    exclude entries where the value is None

    Parameters:
    - data (list): List of data entries from the json
    - key (str): The key to extract from the json

    Returns:
    set: A set of unique values for the specified key.
    """
    return set(entry.get(key) for entry in data if entry.get(key) is not None)


def save_as_txt(data_set, file_name):
    """
    Save the elements of a set to a text file in the 'txt' folder and
    print the number of elements.

    Parameters:
    - data_set (set): The set to be saved.
    - file_name (str): The name of the file to save the set.

    Returns:
    None
    """
    folder_name = 'txt'
    os.makedirs(folder_name, exist_ok=True)  # Create the 'txt' folder
    file_path = os.path.join(folder_name, f"{file_name}.txt")

    with open(file_path, "w") as output_file:
        output_file.write("\n".join(data_set))

    print(f"{len(data_set)} elements saved to {file_path}")


def count_matching_inchikeys(lipids_data, metabolites_data, save_txt):
    """
    Compare and find common InChiKeys, SMILES, and chemical formulas
    between lipid and metabolite datasets loaded from JSON files.

    Parameters:
    - lipids_file (str): Path to the .json containing lipids data.
    - metabolites_file (str): Path to the .json containing metabolites data.
    - save_txt (bool): Save all InChiKeys as .txt file into a new folder.

    Returns:
    tuple: A tuple containing sets of common InChiKeys, SMILES and
    chemical formulas found in both lipid and metabolite datasets.
    """
    # Extract InChiKey values from data
    lipid_inchikeys = extract_values(lipids_data, "INCHI_KEY")
    metabolite_inchikeys = extract_values(metabolites_data, "inchikey")

    # Extract Smiles values from lipids data
    lipid_smiles = extract_values(lipids_data, "SMILES")
    metabolite_smiles = extract_values(metabolites_data, "smiles")

    # Extract chemical formula values from lipids data
    lipid_formula = extract_values(lipids_data, "FORMULA")
    metabolite_formula = extract_values(metabolites_data, "chemical_formula")

    # Save lipid_inchikeys to a text file
    if (save_txt):
        save_as_txt(lipid_inchikeys, 'lipid_inchikeys')
        save_as_txt(metabolite_inchikeys, 'metabolite_inchikeys')

    # Find common keys
    common_inchikeys = lipid_inchikeys.intersection(metabolite_inchikeys)
    common_simles = lipid_smiles.intersection(metabolite_smiles)
    common_formula = lipid_formula.intersection(metabolite_formula)

    return common_inchikeys, common_simles, common_formula


def plot_inchikey_cycles(lipids_data, metabolites_data, save_fig, with_logo):
    """
    Generate a plot with two cycles representing the number of InChiKeys
    in lipids.json and metabolites.json.The Intersection area represents
    the common lipids which are found in both datasets.

    Parameters:
    - lipids_data (list): List of lipids data entries.
    - metabolites_data (list): List of metabolites data entries.
    - save_fig (bool): Save as plot as png.
    - with_logo (bool): Plot intersection diagram with logo.

    Returns:
    None
    """
    common_keys, _, _ = count_matching_inchikeys(lipids_data, metabolites_data,
                                                 save_txt=False)

    # Calculate the sizes of each set
    total_metabolite_inchikeys = len(extract_values(metabolites_data,
                                                    'inchikey'))
    total_lipid_inchikeys = len(extract_values(lipids_data,
                                               'INCHI_KEY'))
    total_common_inchikeys = len(common_keys)

    # Calculate the proportional sizes by setting the size for metabolites
    size_metabolite = 10
    size_lipid = (size_metabolite * total_lipid_inchikeys /
                  total_metabolite_inchikeys)

    size_common = (size_metabolite * total_common_inchikeys /
                   total_metabolite_inchikeys)

    # Calculate the radius of the circles based on the area
    radius_metabolite = math.sqrt(size_metabolite / math.pi)
    radius_lipid = math.sqrt(size_lipid / math.pi)
    radius_common = math.sqrt(size_common / math.pi)

    # print("\nRadius")
    # print(radius_metabolite)
    # print(radius_lipid)
    # print(radius_common)
    # print("\nSize")
    # print(size_metabolite)
    # print(size_lipid)
    # print(size_common)

    # Calculate the positions by using the formula for the distance d
    # and wolframalpha as seen in the formula.txt file
    pos_metabolite = (0, 0)
    pos_lipid = (0, 2.23718218257649)   # (2.23718, 0)
    pos_common = (0, 1.67411)           # (1.67411, 0)

    # Plotting
    plt.figure(figsize=(7, 8))          # (10, 7)

    # Draw the first circle representing metabolite InChiKeys
    plt.gca().add_patch(plt.Circle(pos_metabolite, radius_metabolite,
                                   color='#D17E1A', alpha=0.6,
                                   label='HMDB Dataset'))
    plt.text(*pos_metabolite,
             f'Metabolites\n'
             f'{total_metabolite_inchikeys:,.0f}'.replace(",", " "),
             ha='center', va='center',
             color='white', size=radius_metabolite*10)

    # Draw the second circle representing lipid InChiKeys
    plt.gca().add_patch(plt.Circle(pos_lipid, radius_lipid,
                                   color='#1F5CA8', alpha=0.5,
                                   label='LIPID MAPS Dataset'))
    plt.text(*pos_lipid,
             f'Lipids\n{total_lipid_inchikeys:,.0f}'.replace(",", " "),
             ha='center', va='center', color='white', size=radius_lipid*15)

    # Draw the intersection representing common InChiKeys
    # plt.gca().add_patch(plt.Circle(pos_common, radius_common,
    #                                color='red', alpha=0.5))
    plt.text(*(pos_common[0], pos_common[1] - 0.08),
             f'Common\n{total_common_inchikeys}',
             ha='center', va='center', color='white',  size=radius_common*25)

    if (with_logo):
        # Load the logo image
        logo_hmdb_img = mpimg.imread('00_Ressources/logo/hmdb_logo.png')
        logo_lm_img = mpimg.imread('00_Ressources/logo/lipid_maps_logo.png')

        # Draw the logo (left, right, bottom, top)
        plt.imshow(logo_hmdb_img, extent=[
            pos_metabolite[0] - 0.6, pos_metabolite[0] + 0.8,
            pos_metabolite[1] - 1.5, pos_metabolite[1] - 0.5],
            zorder=10, alpha=0.8)
        plt.imshow(logo_lm_img, extent=[
            pos_lipid[0] - 0.35, pos_lipid[0] + 0.35,
            pos_lipid[1] + 0.75, pos_lipid[1] + 0.25])  # change + to -

    # Set plot limits and show the plot
    plt.xlim(-2, 2)         # (-2, 3.3)
    plt.ylim(-2.5, 3.3)     # (-2, 2)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('on')  # This ensures the axis is turned on
    plt.tick_params(axis='both', which='both', length=0)  # Hide ticks
    plt.xticks([])  # Hide x-axis labels
    plt.yticks([])  # Hide y-axis labels
    plt.legend(loc='lower center')
    # plt.grid(alpha=0.5)
    if (save_fig):
        figures_folder = 'figures'
        os.makedirs(figures_folder, exist_ok=True)  # Create the folder
        if (with_logo):
            plt.savefig(os.path.join(figures_folder,
                                     'intersection_plot_with_logo.png'),
                        format='png', dpi=500, bbox_inches='tight')
        else:
            plt.savefig(os.path.join(figures_folder,
                                     'intersection_plot_without_logo.png'),
                        format='png', dpi=500, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    # Inputs for saving the plot, using a logo, and saving the InChiKeys as txt
    save_fig = False
    with_logo = False
    save_txt = False

    lipids_file_path = "data/lipids.json"
    metabolites_file_path = "data/metabolites.json"

    # Load data from JSON files
    with open(lipids_file_path, 'r') as lipids_json:
        lipids_data = json.load(lipids_json)

    with open(metabolites_file_path, 'r') as metabolites_json:
        metabolites_data = json.load(metabolites_json)

    # Count the sizes of each set
    total_metabolite_inchikeys = len(extract_values(metabolites_data,
                                                    'inchikey'))
    total_lipid_inchikeys = len(extract_values(lipids_data,
                                               'INCHI_KEY'))

    # Count matching keys
    common_inchikeys, _, _ = count_matching_inchikeys(lipids_data,
                                                      metabolites_data,
                                                      save_txt)

    # Calulate percentage of common lipids relative to all metabolites
    percent = round((100 * len(common_inchikeys) /
                     total_metabolite_inchikeys), 2)

    # Print number of all InChiKeys of both datasets
    print("\n#################### NUMBER OF InChiKeys ####################")
    print(f"\t... in metabolite dataset: \t"
          f"{total_metabolite_inchikeys:,.0f}".replace(",", " "))
    print(f"\t... in lipid dataset: \t\t"
          f"{total_lipid_inchikeys:,.0f}".replace(",", " "))

    # RESULTS
    print("\n########################## RESULTS ##########################")
    # Print the number of matching keys
    print(f"Number of matching InChiKeys: {len(common_inchikeys)}")
    # Print the percentage of common lipids reelative to all human metabolites.
    print(f"Percentage of lipids relative to all metabolites: {percent} %")

    # Generate and display the plot
    plot_inchikey_cycles(lipids_data, metabolites_data,
                         save_fig, with_logo)
