import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

# 1. Parse command line arguments from Nextflow
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help="Input merged assoc file")
parser.add_argument('--out', required=True, help="Output image file")
args = parser.parse_args()

# 2. Load and clean data
print(f"Loading data from {args.input}...")
df = pd.read_csv(args.input, sep='\s+')
df = df.dropna(subset=['P'])
df['minus_log10_p'] = -np.log10(df['P'])

# 3. Process X-axis for all 22 chromosomes
# We need to calculate a continuous X-axis so chromosomes don't overlap
df['CHR'] = df['CHR'].astype(int)
df = df.sort_values(['CHR', 'BP'])

df['ind'] = range(len(df)) # Create a continuous index for plotting
df_grouped = df.groupby('CHR')

# 4. Plotting
fig = plt.figure(figsize=(16, 6))
ax = fig.add_subplot(111)
colors = ['royalblue', 'lightsteelblue']
x_labels = []
x_labels_pos = []

for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='minus_log10_p', 
               color=colors[num % len(colors)], ax=ax, s=8, alpha=0.8)
    # Record the middle position of each chromosome for the X-axis label
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))
    x_labels.append(str(name))

# 5. Styling
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
ax.set_ylim([0, max(df['minus_log10_p']) + 5])

plt.axhline(y=7.3, color='red', linestyle='--', label='Significance Threshold (5e-8)')
plt.xlabel('Chromosome', fontsize=12)
plt.ylabel('-log10(P-value)', fontsize=12)
plt.title('GWAS Manhattan Plot (Chr1-22): EAS vs EUR Population', fontsize=14)
plt.legend()
plt.grid(axis='y', alpha=0.3)

# 6. Save
plt.tight_layout()
plt.savefig(args.out, dpi=300)
print(f"✅ Full Genome Manhattan plot saved as {args.out}")