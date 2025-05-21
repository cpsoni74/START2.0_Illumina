import pandas as pd

# Load the CSV file into a DataFrame
input_csv = '/Users/chintansoni/Desktop/NGS/TWIST analysis/peptoid/ncAA>0, noncAA>=0/stdev/peptoid_avg_en_noncAA_PB.csv'  
output_csv = '/Users/chintansoni/Desktop/NGS/TWIST analysis/peptoid/ncAA>0, noncAA>=0/stdev/peptoid_avg_en_noncAA_PB_grouped_stdev.csv'
output2_csv = '/Users/chintansoni/Desktop/NGS/TWIST analysis/peptoid/ncAA>0, noncAA>=0/stdev/peptoid_avg_en_noncAA_PB_avgd_stdev.csv'

df = pd.read_csv(input_csv)

# Group by 'aaRS' and aggregate data into lists for calculating std
grouped_df = df.groupby('aaRS').agg(
    barcodes=('barcode', lambda x: list(x)),
    rawbccountsin=('Rawbccounts_in', lambda y: list(y)),
    aaRS_count=('aaRS', 'size'),
    sum_average_ncAA=('average_ncAA', 'sum'),
    sum_average_noncAA=('average_noncAA', 'sum'),
    list_average_ncAA=('average_ncAA', lambda x: list(x)),
    list_average_noncAA=('average_noncAA', lambda x: list(x))
).reset_index()

# Calculate average values
grouped_df['average_ncAA'] = grouped_df['sum_average_ncAA'] / grouped_df['aaRS_count']
grouped_df['average_noncAA'] = grouped_df['sum_average_noncAA'] / grouped_df['aaRS_count']

# Calculate standard deviation for ncAA and noncAA
grouped_df['std_ncAA'] = grouped_df['list_average_ncAA'].apply(lambda x: pd.Series(x).std())
grouped_df['std_noncAA'] = grouped_df['list_average_noncAA'].apply(lambda x: pd.Series(x).std())

# Save the full grouped data (including std columns) to CSV
grouped_df.to_csv(output_csv, index=False)

# Prepare the reduced columns version (including std) for output2
df2 = grouped_df[['aaRS', 'barcodes', 'rawbccountsin', 'aaRS_count', 'average_ncAA', 'std_ncAA', 'average_noncAA', 'std_noncAA']]

df2.to_csv(output2_csv, index=False)

print("Grouped data with standard deviation saved to", output_csv)
print("Averaged data with standard deviation saved to", output2_csv)
