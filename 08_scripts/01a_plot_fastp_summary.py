import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the summary CSV file
df = pd.read_csv("02_fastqc_reports/fastp_summary.csv", sep="\t")

# Convert all relevant columns to numeric
target_cols = ['Total_Reads', 'Passed_Reads', 'Q20_Rate', 'Q30_Rate', 'GC_Content']
for col in target_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

sns.set(style="whitegrid")

# Plot Q20 vs Q30
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df, x="Q20_Rate", y="Q30_Rate", hue="Sample", s=60)
plt.title("Q20 vs Q30")
plt.savefig("07_plots/q20_vs_q30.png")
plt.close()

# Plot Total Reads vs Passed Reads
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df, x="Total_Reads", y="Passed_Reads", hue="Sample", s=60)
plt.title("Total Reads vs Passed Reads")
plt.savefig("07_plots/total_vs_passed.png")
plt.close()

# Plot GC content per sample
plt.figure(figsize=(10, 6))
sns.barplot(x="Sample", y="GC_Content", data=df)
plt.xticks(rotation=90)
plt.title("GC Content per Sample")
plt.tight_layout()
plt.savefig("07_plots/gc_content.png")
plt.close()
