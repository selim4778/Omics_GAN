import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os

# Load your PSF and DCD files
os.chdir('C:/Users/HP/OneDrive - Tulane University/Desktop\Other_project/Ashad_alam/4th_project/colon/feature_selection')


data_original = pd.read_csv('mRNA_T.csv', index_col=0)
data_original = data_original.T
data_synthetic  = pd.read_csv('run3omics1_4.csv', index_col=0)
label = pd.read_csv('label.csv', index_col=0)


x1 = pd.concat([data_original, label], axis=1)
x2 = pd.concat([data_synthetic, label], axis=1)

## original data 
df = pd.DataFrame(x1)
# Perform Wilcoxon test for each feature
# Perform Wilcoxon rank-sum test (Mann-Whitney U test)  for original data 
selected_features = []
p_values = {}
for feature in df.columns[:-1]:  # Exclude 'CaseControl' column
    case_values = df[df['EVENT'] == 1][feature]
    control_values = df[df['EVENT'] == 0][feature]
    
    # Mann-Whitney U test (Wilcoxon rank-sum test)
    stat, p_value = mannwhitneyu(case_values, control_values, alternative='two-sided')
    p_values[feature] = p_value
    if p_value < 0.005:  # Feature selection threshold
        selected_features.append(feature)

# Convert to DataFrame
p_values_df = pd.DataFrame.from_dict(p_values, orient='index', columns=['p-value'])
print(p_values_df)

selected_features_set = set(selected_features)

# Original data significant features with p-values
selected_features_df = pd.DataFrame({
    'Feature': list(selected_features),
    'p-value': [p_values_df.loc[f, 'p-value'] for f in selected_features]
})
selected_features_df.to_csv("colon_selected_features_original.csv", index=False)

## synthetic data 
df = pd.DataFrame(x2)
# Perform Wilcoxon test for each feature
# Perform Wilcoxon rank-sum test (Mann-Whitney U test)  for synthetic data  
selected_features_syn = []
p_values = {}
for feature in df.columns[:-1]:  # Exclude 'CaseControl' column
    case_values = df[df['EVENT'] == 1][feature]
    control_values = df[df['EVENT'] == 0][feature]
    
    # Mann-Whitney U test (Wilcoxon rank-sum test)
    stat, p_value = mannwhitneyu(case_values, control_values, alternative='two-sided')
    p_values[feature] = p_value
    if p_value < 0.005:  # Feature selection threshold
        selected_features_syn.append(feature)

# Convert to DataFrame
p_values_df = pd.DataFrame.from_dict(p_values, orient='index', columns=['p-value'])
print(p_values_df)

# Draw Venn Diagram for common features
sig_F_ori  = set(selected_features)
sig_F_syn = set(selected_features_syn)


plt.figure(figsize=(5, 5))
venn = venn2([sig_F_ori, sig_F_syn], ('Original-data', 'Synthetic-data'))
plt.title("Venn Diagram of Selected Features")
plt.savefig('colon_feature.png', dpi=300)
plt.show()


# Find and print common features
common_features = sig_F_ori.intersection(sig_F_syn)
print("Common features between both datasets:", common_features)

# Save common features to a CSV file
common_features_df = pd.DataFrame(list(common_features), columns=['Common_Features'])
common_features_df.to_csv("colon_features.csv", index=False)


####################### save significant gene with p-value

# Synthetic data significant features with p-values
selected_features_syn_df = pd.DataFrame({
    'Feature': list(selected_features_syn),
    'p-value': [p_values_df.loc[f, 'p-value'] for f in selected_features_syn]
})
selected_features_syn_df.to_csv("colon_selected_features_synthetic.csv", index=False)
