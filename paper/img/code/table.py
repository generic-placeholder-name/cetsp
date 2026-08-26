import pandas as pd

# --- 1. Load CSV ---
df = pd.read_csv("comparison.csv")

# --- 2. Extract instance name from filename ---
df['Instance'] = df['filename'].str.replace(r'\.txt$', '', regex=True)

# --- 3. Select and round columns ---
df = df[['Instance', 'N', 'OR%', 'benchmark_value', 'output_value', 'percent_diff']]
df['OR%'] = df['OR%'].map("{:.2f}".format)           # 2 digits
df['percent_diff'] = df['percent_diff'].map("{:.2f}".format)  # 2 digits
df['benchmark_value'] = df['benchmark_value'].map("{:.3f}".format)  # 3 digits
df['output_value'] = df['output_value'].map("{:.3f}".format)        # 3 digits

# --- 4. Rename columns for LaTeX ---
df.columns = ['Instance', 'N', 'OR\\%', 'Benchmark Value', 'Output Value', 'Gap\\%']

# --- 5. Generate LaTeX table rows only ---
# index=False avoids row numbers, header=False avoids a nested tabular header
table_rows = df.to_latex(index=False, header=False, escape=True)

# --- 6. Build longtable preamble and footer ---
latex_table = r"""
\begin{longtable}{lrrrrr}
\caption{Comparison of benchmark vs output values for all instances} \label{tab:comparison} \\
\hline
Instance & N & OR\% & Benchmark Value & Output Value & Gap\% \\
\hline
\endfirsthead

\multicolumn{6}{c}{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
\hline
Instance & N & OR\% & Benchmark Value & Output Value & Gap\% \\
\hline
\endhead

\hline \multicolumn{6}{r}{{Continued on next page}} \\
\endfoot

\hline
\endlastfoot
""" + table_rows + r"\end{longtable}"

# --- 7. Write to .txt file ---
with open("comparison_table.txt", "w") as f:
    f.write(latex_table)

print("LaTeX longtable written to comparison_table.txt successfully.")
