"""Summarising each results page into
the correct format"""

from pathlib import Path

import pandas as pd


def find_filepath():
    """Finds the correct filepaths for all folders with Pinhole results"""

    folder_path = Path(__file__).parent / r"outputs/PinHoleDataset_SIMPER_output_run2"
    filepath_list = []
    for child_folderpath in Path.iterdir(folder_path):
        # find the filepaths in the folder
        if child_folderpath.is_dir():
            for filepath in Path.iterdir(child_folderpath):
                filepath_list.append(filepath)

    return filepath_list


def get_results_info(sample_path):
    """Gets the correct result information
    and adds it to a list or dataframe"""

    file_name = sample_path.name
    sample_name = file_name.split("_re")[0]

    old_id = sample_path.parent.stem

    results_df = pd.read_csv(sample_path)
    # get top score
    top_score = results_df.loc[0, "Match"]

    # find the top species/ multiple top species
    for index, row in results_df.iterrows():
        if index == 0:
            species_string = row["SPECIES"]
            species_score = row["Match"]
        elif row["Match"] == species_score:
            species_string += f"/{row['SPECIES']}"
        elif row["Match"] < species_score:
            break
    return (sample_name, old_id, species_string, top_score)


def main():
    """Runs the commands to get the results
    in the correct format"""

    fp_list = find_filepath()
    results_df = pd.DataFrame(
        columns=["Sample Name", "Previous ID", "SIF Top Species", "SIF Top Score"]
    )
    for filepath in fp_list:
        results_df.loc[len(results_df)] = get_results_info(filepath)
    print(results_df.shape)
    print(results_df.head())
    print(results_df.tail())
    results_df.to_csv("test_results2.csv")



if __name__ == "__main__":
    main()
