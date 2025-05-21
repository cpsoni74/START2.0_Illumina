# Wizard: Chintan Soni
# at least n-fold diff between +ncM and -ncM

import pandas as pd

# Input and output file paths
inputfile = '/Users/chintansoni/Desktop/NGS/TWIST analysis/peptoid/ncAA>0, noncAA>=0/stdev/peptoid_avg_en_noncAA_PB_avgd_stdev.csv'
atn_file = '/Users/chintansoni/Desktop/NGS/TWIST analysis/peptoid/ncAA>0, noncAA>=0/stdev/peptoid_avg_en_noncAA_PB_avgd_stdev_at5.csv'
atnminn_file = '/Users/chintansoni/Desktop/NGS/TWIST analysis/peptoid/ncAA>0, noncAA>=0/stdev/peptoid_avg_en_noncAA_PB_avgd_stdev_at5min20.csv'

# Load the data
df = pd.read_csv(inputfile)

# Apply initial filter: average_noncAA <= average_ncAA / 10
filtered_df = df[df['average_noncAA'] <= (df['average_ncAA'] / 10)]

# NEW: Filter to only keep rows where std_ncAA <= 2 * average_ncAA
filtered_df = filtered_df[filtered_df['std_ncAA'] <= 1.5 * filtered_df['average_ncAA']]

filtered_df = filtered_df[filtered_df['std_ncAA'] <= 200]
# Save intermediate file (after applying the std_ncAA filter)
filtered_df.to_csv(atn_file, index=False)

# Apply additional filter: average_ncAA >= 10
df2 = filtered_df[filtered_df['average_ncAA'] >= 20]

# Save final file
df2.to_csv(atnminn_file, index=False)

print(f"Filtered data saved to: {atn_file}")
print(f"Filtered data (ncAA >= 10) saved to: {atnminn_file}")
