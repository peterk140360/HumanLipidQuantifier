###############################################################################
# Author:      Sebastian Peterka
# Company:     -
# File:        MetaboliteVsInchikey.py
# created:     06.12.2023 11:30
# edited:      -
#
# Description: The Python analyzes two text files containing line numbers,
#              calculating and visualizing the distances between
#              corresponding lines, along with individual data points,
#              using matplotlib subplots for a comprehensive overview.
#
# Comments:    Used Datasets:
#              * HMDB ->  All Metabolites (in XML format)
#                              Released on: 2021-11-17
#                              File Size: 910 MB
#              * First Plot:  The first plot displays the distances between
#                             corresponding lines in "Inchikey" and "Metabolite
#                             using purple markers to represent the distances.
#              * Second Plot: The second plot visualizes individual data points
#                             in the "Inchikey" and "Metabolite" text files,
#                             with orange and blue markers.
#
# Usage:       1. Download the mentioned datasets from the resource pages
#              2. Extract them
#              3. Open notepad
#              4. Search for "<metabolite>"
#              5. Click on "Find all in current Document"
#              6. Copy the text into a .txt File
#              7. Delete the first two lines, so that the file starts
#                 with "Line ..."
#              8. Save it and copy the filename in line 87 of this file.
#              9. Do the same for the "<inchikey>" keyword
#             10. Edit the filename in line 86 of this file.
#             11. Run the script
#
# Todo:        *
#              *
#              *
#
# Resources:   https://hmdb.ca/downloads
#
###############################################################################

import matplotlib.pyplot as plt
import os


def calculate_distances(file_path1, file_path2):
    with open(file_path1, 'r') as file1, open(file_path2, 'r') as file2:
        lines1 = file1.readlines()
        lines2 = file2.readlines()

    line_numbers1 = [int(line.split(':')[0].split()[-1]) for line in lines1]
    line_numbers2 = [int(line.split(':')[0].split()[-1]) for line in lines2]

    # Calculate distances
    distances = [line_numbers1[i] - line_numbers2[i]
                 for i in range(min(len(line_numbers1), len(line_numbers2)))]

    return distances


def plot_distances_and_points(distances,
                              inchikey_file_path,
                              metabolite_file_path,
                              save_fig):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    # Plot distances
    ax1.plot(distances, marker='o', color='purple',
             label='Distance between Metabolite and Inchikey')
    ax1.set_title('Distances Between Corresponding Lines '
                  'of Metabolite and Inchikey')
    ax1.set_ylabel('Distance in Lines')
    ax1.legend()

    # Plot individual points
    with open(inchikey_file_path, 'r') as inchikey_file, \
         open(metabolite_file_path, 'r') as metabolite_file:
        inchikey_lines = inchikey_file.readlines()
        metabolite_lines = metabolite_file.readlines()

    inchikey_line_numbers = [int(line.split(':')[0].split()[-1])
                             for line in inchikey_lines]
    metabolite_line_numbers = [int(line.split(':')[0].split()[-1])
                               for line in metabolite_lines]

    ax2.plot(inchikey_line_numbers, 'o', color='orange',
             label='Inchikey Datapoints')
    ax2.plot(metabolite_line_numbers, 'o', color='blue',
             label='Metabolite Datapoints')
    ax2.set_title('Individual Datapoints for Metabolite and Inchikey')
    ax2.set_xlabel('Line Numbers')
    ax2.set_ylabel('')

    # Adjust y-axis ticks to display labels without exponential notation
    ax2.ticklabel_format(axis='y', style='plain')

    ax2.legend()
    if (save_fig):
        figures_folder = 'figures'
        os.makedirs(figures_folder, exist_ok=True)  # Create the folder
        plt.savefig(os.path.join(figures_folder, 'distances.png'),
                    format='png', dpi=500, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    inchikey_file_path = "data/inchikey_217899_hits.txt"
    metabolite_file_path = "data/metabolite_217920_hits.txt"

    distances = calculate_distances(inchikey_file_path, metabolite_file_path)

    # Calculate mean, max, and min distances
    mean_distance = sum(distances) / len(distances)
    max_distance = max(distances)
    min_distance = min(distances)

    # Print results
    print(f"Mean distance between Metabolite and Inchikey:"
          f"{mean_distance} Lines.")
    print(f"Max distance between Metabolite and Inchikey:"
          f"{max_distance} Lines.")
    print(f"Min distance between Metabolite and Inchikey:"
          f"{min_distance} Lines.")

    # Print results for the first 5 distances
    print("Results for the first 5 distances:")
    for i, distance in enumerate(distances[:5], start=1):
        print(f"Distance {i}: {distance} Lines.")

    # Plot distances and individual points
    plot_distances_and_points(distances,
                              inchikey_file_path, metabolite_file_path,
                              save_fig=False)
