import pandas as pd

# Define the file paths
file1_path = '/Users/chintansoni/Desktop/NGS/TWIST analysis/Rb2NH2/Results/Analysis, ncAA>0, noncAA>0 for out/Malu Test/Rb2NH2>0_noncAAox_avg_en.csv'
#file2_path = '/Users/chintansoni/Desktop/NGS/TWIST analysis/Rb2NH2/Results/Analysis, ncAA>0/noncAAox_avg_en.csv'
file2_path = '/Users/chintansoni/Desktop/Run1,2_merged_aaRS_df.csv'

outputfile = '/Users/chintansoni/Desktop/NGS/TWIST analysis/Rb2NH2/Results/Analysis, ncAA>0, noncAA>0 for out/Malu Test/Rb2NH2>0_noncAAox_avg_en_PB.csv'
# Read both CSV files into DataFrames
df1 = pd.read_csv(file1_path)
df2 = pd.read_csv(file2_path)

# Drop duplicates from df2 based on the 'barcode' column
df2 = df2.drop_duplicates(subset='barcode')

# Find common barcodes between the two DataFrames
common_barcodes = set(df1['barcode']).intersection(set(df2['barcode']))

# Filter rows in df2 that have common barcodes
filtered_df2 = df2[df2['barcode'].isin(common_barcodes)]

# Merge data from df1 and filtered_df2 based on the 'barcode' column
merged_df = df1.merge(filtered_df2, on='barcode', how='inner')

# Save the merged DataFrame to a new CSV file
merged_df.to_csv(outputfile, index=False)

print("Merged data saved successfully.")
