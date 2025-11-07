# Stage One Task: Python Deep Dive

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 
from scipy.stats import ttest_ind
import numpy as np

# Translating DNA to protein
'''
Generate a dictionary containing codon to amino acid mappings.
Translate a given DNA sequence into its corresponding amino acid sequence.
'''

dic = {'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
'UAC':'Y', 'UAU':'Y', 'UAA':'STOP', 'UAG':'STOP',
'UGC':'C', 'UGU':'C', 'UGA':'STOP', 'UGG':'W'}

def translate_dna_protein(codon_table, dna_seq):
    protein = ''
    # Convert the input sequence from  DNA to RNA
    rna_seq = dna_seq.replace('T', 'U')
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i+3]
        # break the loop if the codon length is less than 3
        if len(codon) < 3:
            break

        aa = codon_table.get(codon)

        # Stop translation if it is a STOP  codon 
        if aa == 'STOP':
            break

        # Skip unknown codons
        if aa:
            protein += aa

    return protein


# 2.Calculating the hamming distance
'''
Calculate the Hamming distance between your username in HacBio and another username from social media  of equal length.
'''

def measure_hamming_distance(seq1, seq2):
    # Raise error if lengths are not equal
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    # Initialize distance counter
    distance = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            distance += 1

    return distance

# --- Script execution ---
if __name__ == "__main__":
    # DNA to protein translation
    dna_seq = 'ATGGCCTTTGACGATGTTGACCTGACCGGCTACTTCTTGAAGGCGTTCGTACTGGGCCATGAACCGTGGCGTATTTGAGGCTTAA'
    protein_seq = translate_dna_protein(dic, dna_seq)
    print("Protein sequence:", protein_seq)

    # Hamming distance
    first = 'anyavala'
    second = 'ozlemkal'
    distance = measure_hamming_distance(first, second)
    print(f"Hamming distance between '{first}' and '{second}': {distance}")


print('---------------------------------------------------------------------')
print('Gene Expression Analysis - Heatmap of Top Differentially Expressed Genes')
#Part A – Gene Expression Analysis
#a. Heatmap
"""
    Plot a clustered heatmap of top differentially expressed genes.
    
    Parameters:
    - file_path: str, path to CSV file (genes x samples)
    - groups: dict, e.g., {"HBR": ["HBR_1","HBR_2"], "UHR": ["UHR_1","UHR_2"]}
    - top_n: int, number of top genes to select
    - cmap: str, color palette
    - normalize_rows: bool, whether to normalize rows (genes) to 0-1
"""

df=pd.read_csv('normalized_counts_for_HBR_and_UHR.csv', index_col=0)
# Sample groups
groups = {
    "HBR": ["HBR_1", "HBR_2", "HBR_3"],
    "UHR": ["UHR_1", "UHR_2", "UHR_3"]
}

# Compute mean expression per group
df["HBR_mean"] = df[groups["HBR"]].mean(axis=1)
df["UHR_mean"] = df[groups["UHR"]].mean(axis=1)

# Compute fold change (absolute difference or log2)
df["diff"] = np.abs(df["HBR_mean"] - df["UHR_mean"])

# Select top N genes (e.g., top 20)
top_genes = df.nlargest(20, "diff")
heatmap_data = top_genes[groups["HBR"] + groups["UHR"]]
# Set a color palette, e.g., Blues
sns.set(font_scale=1)
cg = sns.clustermap(
    heatmap_data,
    cmap="Blues",
    standard_scale=1,   # normalize rows to 0-1
    row_cluster=True,   # cluster genes
    col_cluster=True,   # cluster samples
    linewidths=0.5,
    figsize=(10, 10)
)

# Add labels
cg.ax_heatmap.set_xlabel("Samples")
cg.ax_heatmap.set_ylabel("Genes")
plt.show()  

print('---------------------------------------------------------------------')
print('Volcano Plot of Differentially Expressed Genes (DEGs)')

'''
    Create a volcano plot to visualize DEGs between HBR and UHR samples.
    Highlight significantly upregulated and downregulated genes.
'''
df1=pd.read_csv('differential_expression_results.csv')
def significance(row):
    if row["PAdj"] < 0.05 and row["log2FoldChange"] > 1:
        return "Upregulated"
    elif row["PAdj"] < 0.05 and row["log2FoldChange"] < -1:
        return "Downregulated"
    else:
        return "Not significant"

df1["significance"] = df1.apply(significance, axis=1)
# Prepare colors
colors = {
    "Upregulated": "green",
    "Downregulated": "orange",
    "Not significant": "grey"
}

# Scatter plot
plt.figure(figsize=(10,6))
for sig in ["Upregulated", "Downregulated", "Not significant"]:
    subset = df1[df1["significance"] == sig]
    plt.scatter(
        subset["log2FoldChange"], 
        -np.log10(subset["PAdj"]),
        c=colors[sig],
        label=sig,
        alpha=0.7,
        edgecolors="k",
        s=50
    )

# Add dashed vertical lines at log2FC = ±1
plt.axvline(x=1, color='black', linestyle='--')
plt.axvline(x=-1, color='black', linestyle='--')

# Labels and title
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(adjusted p-value)")
plt.title("Volcano Plot of DEGs (HBR vs UHR)")
plt.legend()
plt.show()

print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - Scatter Plot of Radius vs Texture')
df2=pd.read_csv('breast_cancer_wisconsin.csv')
df.head()
# Map diagnosis to colors
colors = {'M': 'red', 'B': 'blue'}
point_colors = df2['diagnosis'].map(colors)

plt.figure(figsize=(8,6))
plt.scatter(
    df2['radius_mean'], 
    df2['texture_mean'], 
    c=point_colors, 
    alpha=0.7,   # slightly transparent for overlapping points
    edgecolor='k' # black edge for better visibility
)

# Add legend manually
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Malignant', markerfacecolor='red', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Benign', markerfacecolor='blue', markersize=8)
]
plt.legend(handles=legend_elements, title="Diagnosis")

# Labels and title
plt.xlabel('Radius Mean')
plt.ylabel('Texture Mean')
plt.title('Breast Cancer: Radius vs Texture')
plt.show()

print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - Correlation Heatmap of Key Features')


# Select the six key features
features = ['radius_mean', 'texture_mean', 'perimeter_mean', 
            'area_mean', 'smoothness_mean', 'compactness_mean']

#Calculate correlation matrix

df_features = df2[features]
corr_matrix = df_features.corr()
plt.figure(figsize=(8,6))
sns.heatmap(
    corr_matrix, 
    annot=True,         # show correlation values
    cmap='Blues',    # color gradient: blue (-1) → red (+1)
    fmt=".2f",          # format numbers to 2 decimal places
    linewidths=0.5      # grid lines between cells
)

plt.title("Correlation Heatmap of Key Breast Cancer Features")
plt.show()
print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - Scatter Plot of Compactness vs Smoothness')

# Map diagnosis to colors
colors = {'M': 'red', 'B': 'blue'}
point_colors = df2['diagnosis'].map(colors)

plt.figure(figsize=(8,6))
plt.scatter(
    df2['smoothness_mean'],       # X-axis
    df2['compactness_mean'],      # Y-axis
    c=point_colors,
    alpha=0.7,
    edgecolor='k'                # black edge for better visibility
)

# Add gridlines
plt.grid(True, linestyle='--', alpha=0.5)

# Axis labels and title
plt.xlabel('Smoothness Mean')
plt.ylabel('Compactness Mean')
plt.title('Breast Cancer: Compactness vs Smoothness')

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Malignant', markerfacecolor='red', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Benign', markerfacecolor='blue', markersize=8)
]
plt.legend(handles=legend_elements, title="Diagnosis")
plt.show()

print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - KDE of Area Mean by Diagnosis')

plt.figure(figsize=(8,6))

# KDE for malignant
sns.kdeplot(
    data=df2[df2['diagnosis']=='M'],
    x='area_mean',
    fill=True,             # fill under the curve
    color='red',
    alpha=0.5,
    label='Malignant'
)

# KDE for benign
sns.kdeplot(
    data=df2[df2['diagnosis']=='B'],
    x='area_mean',
    fill=True,
    color='blue',
    alpha=0.5,
    label='Benign'
)

# Labels and title
plt.xlabel('Area Mean')
plt.ylabel('Density')
plt.title('Distribution of Area Mean by Diagnosis')
plt.legend(title='Diagnosis')
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()




