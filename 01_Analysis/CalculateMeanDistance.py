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


if __name__ == "__main__":
    file_path = "data/metabolite_217920_hits.txt"
    (mean_distance, max_distance, min_distance,
     first_5_distances, min_distance_line,
     max_distance_line, all_distances) = calculate_distances(file_path)

    print(f"Mean distance is {mean_distance} Lines.")
    print(f"Max distance is {max_distance} Lines.")
    print(f"Min distance is {min_distance} Lines.")

    for i, distance in enumerate(first_5_distances, start=1):
        print(f"Distance between lines {i} and {i + 1} is {distance} Lines.")

    print(f"Line with the minimum distance: {min_distance_line}")
    print(f"Line with the maximum distance: {max_distance_line}")

    # Find distances over 5000
    over_5000_distances = [distance for distance in all_distances
                           if distance > 5000]
    # Count the number of distances over 5000
    count_over_5000 = len(over_5000_distances)
    print(f"Number of distances over 5000: {count_over_5000}")

    # Plot distances
    plot_distances(all_distances, save_fig=False)
