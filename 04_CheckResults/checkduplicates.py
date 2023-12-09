def find_duplicates(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Strip newline characters and create a set to check for duplicates
        stripped_lines = [line.strip() for line in lines]
        seen = set()

        for line_number, line in enumerate(stripped_lines, start=1):
            if line in seen:
                print(f'Duplicate found on line {line_number}: {line}')
            else:
                seen.add(line)

    except FileNotFoundError:
        print(f"Error: File not found - {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")


find_duplicates('data/metabolite_inchikeys.txt')
print("\n")
find_duplicates('data/lipid_inchikeys.txt')
