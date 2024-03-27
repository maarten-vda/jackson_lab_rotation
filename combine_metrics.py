import argparse
import os
import pandas as pd
import json

def extract_scores(file_path, input_directory):
    # Read the JSON file
    file_path = os.path.join(input_directory, file_path)
    try:
        with open(file_path, "r") as file:
            data = json.load(file)

        # Access values if the keys exist
        if 'ptm' in data and 'iptm' in data:
            return({'ptm': data['ptm'], 'iptm': data['iptm']})
        else:
            return("Keys 'ptm' and 'ipTM' not found in the JSON data.")
    except Exception as e:
        return(f"Error reading the JSON file: {e}")


def get_score_files(directory_path):
    # Get a list of all files in the directory
    all_files = os.listdir(directory_path)

    # Filter files that contain "scores" and end with ".json"
    score_files = [file for file in all_files if "scores" in file and file.endswith(".json")]
    return score_files


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Process data from CSV and input directory")
    parser.add_argument("input_csv", help="Input CSV file")
    parser.add_argument("input_directory", help="Input directory")
    parser.add_argument("output_csv", help="Output CSV file")
    args = parser.parse_args()

    json_files = get_score_files(args.input_directory)
    complexes_dict = {}

    # Extract scores from JSON files and store them in a complexes_dict
    for i in json_files:

        name = i.split('_scores_rank')[0]


        ptm_score = (extract_scores(i, args.input_directory)['ptm'])
        iptm_score = (extract_scores(i, args.input_directory)['iptm'])
        ranking_confidence = 0.2*ptm_score + 0.8*iptm_score
        model_number = int(i.split('_rank')[1].split('_')[1])

        #Added the model number to the key since the same complex can have multiple models
        complexes_dict[f"{name}_{model_number}"] = {'ptm': ptm_score, 'iptm': iptm_score, 'model': model_number, 'ranking_confidence': ranking_confidence}

    # Read the input CSV file
    df = pd.read_csv(args.input_csv)

    # Add ptm and iptm columns to the dataframe
    for key, value in complexes_dict.items():
        complex_name = key.rsplit('_', 1)[0]
        print(complex_name)
        df.loc[(df['complex_name'] == complex_name) & (df['model_num'] == value['model']), 'pTM'] = value['ptm']
        df.loc[(df['complex_name'] == complex_name) & (df['model_num'] == value['model']), 'ipTM'] = value['iptm']
        df.loc[(df['complex_name'] == complex_name) & (df['model_num'] == value['model']), 'ranking_confidence'] = value['ranking_confidence']

    df.rename(columns={'model_num': 'rank'}, inplace=True)
    df.to_csv(args.output_csv, index=False)
    print(f"Modified data saved to {args.output_csv}")

