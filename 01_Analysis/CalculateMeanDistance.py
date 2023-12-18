###############################################################################
# Author:      Sebastian Peterka
# Company:     -
# File:        CalculateMeanDistance.py
# created:     06.12.2023 11:30
# edited:      -
#
# Description: The Python script analyzes a text file containing structured
#              data, calculates and prints various statistics such as
#              mean, maximum, and minimum distances between consecutive lines,
#              generates a line chart representing these distances,
#              and counts the number of distances over 5000.
#
# Comments:    Used Datasets:
#              * HMDB ->  All Metabolites (in XML format)
#                              Released on: 2021-11-17
#                              File Size: 910 MB
#
# Usage:       1. Download the mentioned datasets from the resource pages
#              2. Extract them
#              3. Open notepad
#              4. Search for "<metabolite>"
#              5. Click on "Find all in current Document"
#              6. Copy the text into a .txt File
#              7. Delete the first two lines, so that the file starts
#                 with "Line ..."
#              8. Save it and copy the filename in line 71 of this file.
#              9. Run the script
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


def calculate_distances(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    line_numbers = [int(line.split(':')[0].split()[-1]) for line in lines]

    # Calculate distances
    distances = [line_numbers[i+1] - line_numbers[i]
                 for i in range(len(line_numbers)-1)]

    # Calculate mean, max, and min distances
    mean_distance = sum(distances) / len(distances)
    max_distance = max(distances)
    min_distance = min(distances)

    # Calculate distances between the first 5 elements separately
    first_5_distances = [line_numbers[i+1] - line_numbers[i] for i in range(4)]

    # Find lines for min and max distances
    min_distance_line = lines[distances.index(min_distance)].strip()
    max_distance_line = lines[distances.index(max_distance)].strip()

    return (mean_distance, max_distance, min_distance, first_5_distances,
            min_distance_line, max_distance_line, distances)


def plot_distances(distances, save_fig):
    plt.plot(distances, marker='o')
    plt.title('Distances between Metabolites')
    plt.xlabel('Number of matches in "hmdb_matabolites.xml"'
               'on "<metabolites>"')
    plt.ylabel('Distance in Lines')
    if (save_fig):
        figures_folder = 'figures'
        os.makedirs(figures_folder, exist_ok=True)  # Create the folder
        plt.savefig(os.path.join(figures_folder, 'distances_metabolites.png'),
                    format='png', dpi=500, bbox_inches='tight')
    plt.show()


def plot_boxplot(ax, data, title, xlabel, ylabel):
    # Create a vertical boxplot with whisker caps and outliers
    ax.boxplot(data, sym='k+', showcaps=True, showfliers=True)

    # Set x-axis to logarithmic scale
    ax.set_yscale('log')

    # Add grid lines
    ax.grid(True, linestyle='--', alpha=0.7)

    # Add labels and title
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Customize box and whisker colors
    box_colors = ['#1f77b4']
    whisker_colors = ['#7570b3']
    caps_colors = ['#7570b3']
    median_colors = ['#b2df8a']

    ax.boxplot(data, sym='k+', showcaps=True, showfliers=True,
               boxprops=dict(color=box_colors[0]),
               whiskerprops=dict(color=whisker_colors[0]),
               capprops=dict(color=caps_colors[0]),
               medianprops=dict(color=median_colors[0]))

    # Remove ticks on the x-axis
    ax.set_xticks([])


if __name__ == "__main__":
    save_fig = True
    file_path_met = "data/metabolite_217920_hits.txt"
    file_path_lip = "data/lipid_47751_hits.txt"

    # Calculate distances for metabolite dataset
    (mean_distance, max_distance, min_distance,
     first_5_distances, min_distance_line,
     max_distance_line, all_distances) = calculate_distances(file_path_met)

    # Calculate distances for lipid dataset
    (mean_distance_lip, max_distance_lip, min_distance_lip,
     first_5_distances_lip, min_distance_line_lip,
     max_distance_line_lip, all_distances_lip) = calculate_distances(file_path_lip)

    # Print statistics for metabolite dataset
    print("Metabolite Dataset:")
    print(f"Mean distance is {mean_distance} Lines.")
    print(f"Max distance is {max_distance} Lines.")
    print(f"Min distance is {min_distance} Lines.")
    for i, distance in enumerate(first_5_distances, start=1):
        print(f"Distance between lines {i} and {i + 1} is {distance} Lines.")
    print(f"Line with the minimum distance: {min_distance_line}")
    print(f"Line with the maximum distance: {max_distance_line}")
    print(f"Number of distances over 5000: {len([d for d in all_distances if d > 5000])}")

    # Print statistics for lipid dataset
    print("\nLipid Dataset:")
    print(f"Mean distance is {mean_distance_lip} Lines.")
    print(f"Max distance is {max_distance_lip} Lines.")
    print(f"Min distance is {min_distance_lip} Lines.")
    for i, distance in enumerate(first_5_distances_lip, start=1):
        print(f"Distance between lines {i} and {i + 1} is {distance} Lines.")
    print(f"Line with the minimum distance: {min_distance_line_lip}")
    print(f"Line with the maximum distance: {max_distance_line_lip}")
    print(f"Number of distances over 5000: {len([d for d in all_distances_lip if d > 5000])}")

    # Create a subplot with 1 row and 2 columns (arranged horizontally)
    fig, axes = plt.subplots(1, 2, figsize=(7, 9))

    # Set a common title for the entire figure
    # fig.suptitle('Length of Data Entries per Data Record', fontsize=11)
    # fig.patch.set_edgecolor('black')  # Add a black frame around the entire plot
    # fig.patch.set_linewidth(2)  # Set the linewidth of the frame

    # Plot the first boxplot (metabolite)
    plot_boxplot(axes[0], all_distances, '', 'HMDB', 'Length in Lines')

    # Plot the second boxplot (lipid)
    plot_boxplot(axes[1], all_distances_lip, '', 'LIPID MAPS', '')

    # Adjust layout to prevent clipping of labels
    plt.tight_layout()

    if (save_fig):
        figures_folder = 'figures'
        os.makedirs(figures_folder, exist_ok=True)  # Create the folder
        plt.savefig(os.path.join(figures_folder, 'boxplot_vertical.png'), format='png', dpi=500, bbox_inches='tight')

    # Show the plot
    plt.show()
