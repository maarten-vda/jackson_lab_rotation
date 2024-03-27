import argparse
import pandas as pd

def calculate_average(input_path):
    # Read the input CSV into a Pandas DataFrame
    input_data = pd.read_csv(input_path)

    # Process values in column 'complex_name'
    input_data['complex_name'] = input_data['complex_name'].apply(lambda value: '_'.join(sorted(value.split('_'))))

    # Group by the processed values in column 'complex_name' and calculate the average for each group
    grouped_data = input_data.groupby('complex_name').mean()

    grouped_data.drop('rank', axis=1, inplace=True)

    return grouped_data

def split_dataframe_by_confidence(df):
    # Assuming 'ranking_confidence' is the column name
    # Create two dataframes based on the condition
    interactors = df[df['ranking_confidence'] >= 0.7]
    non_interactors = df[df['ranking_confidence'] < 0.7]

    return interactors, non_interactors

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Split input CSV into two output CSVs.')
    parser.add_argument('input_csv', help='Path to the input CSV file')
    parser.add_argument('output1_csv', help='Path to the first output CSV file')
    parser.add_argument('output2_csv', help='Path to the second output CSV file')
    args = parser.parse_args()

    # Split the CSV and write to output files

    grouped_data = calculate_average(args.input_csv)
    interactors, non_interactors = split_dataframe_by_confidence(grouped_data)
    interactors.to_csv(args.output1_csv)
    non_interactors.to_csv(args.output2_csv)

    print(f"CSV successfully split by threshold and saved to {args.output1_csv} and {args.output2_csv}")

if __name__ == "__main__":
    main()
