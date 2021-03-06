import pandas as pd
from scipy.io import mmread

data_mtx = mmread("combinedSRR_output/alevin/quants_mat.mtx").toarray()

cols = open("combinedSRR_output/alevin/quants_mat_cols.txt", "r")
rows = open("combinedSRR_output/alevin/quants_mat_rows.txt", "r")

cols_list = []
rows_list = []

for line in cols:
	cols_list.append(line.strip("\n"))

for line in rows:
	rows_list.append(line.strip("\n"))

df = pd.DataFrame(data_mtx, columns = cols_list, index = rows_list)

df = df.T

x = open("combined_count_mtx.csv", "w")

x.write(df.to_csv())

