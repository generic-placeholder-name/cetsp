import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
import numpy as np

# Load CSV
df = pd.read_csv("comparison.csv")

# Filter out bubbles and rotatingDiamonds
df_filtered = df[~df['filename'].str.contains('bubbles|rotatingDiamonds')]

# Scatter plot
plt.figure(figsize=(7,5))
sns.scatterplot(data=df_filtered, x='OR%', y='percent_diff')

# Fit linear trendline
X = df_filtered[['OR%']].values
y = df_filtered['percent_diff'].values
reg = LinearRegression().fit(X, y)
plt.plot(X, reg.predict(X), color='red', label=f'Trendline: slope={reg.coef_[0]:.3f}')

plt.xlabel('Overlap Ratio (%)')
plt.ylabel('Percent Gap (%)')
plt.title('CETSP: Overlap Ratio vs Percent Gap (excluding Bubbles & RotatingDiamonds)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
