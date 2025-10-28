import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create ordering for color assignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ordering', type=str, required=True, help="output ordering file")
    parser.add_argument('--metadata', type=str, required=True, help="input metadata file for dates extraction")
    args = parser.parse_args()

    column = "date"
    metadata = pd.read_csv(args.metadata, delimiter='\t')

    if column not in metadata.columns:
        print(f"The column '{column}' does not exist in the file.")
        sys.exit(1)

    deflist = metadata[column].dropna().tolist()
    # Store unique values (ordered)
    deflist = sorted(set(deflist))
    if "XXXX-XX-XX" in deflist:
        deflist.remove("XXXX-XX-XX")

    result_df = pd.DataFrame({
        'column': ['date'] * len(deflist),
        'value': deflist
    })

    result_df.to_csv(args.ordering, sep='\t', index=False, header=False)
