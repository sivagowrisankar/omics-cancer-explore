import pandas as pd
from scipy.stats import ttest_ind, false_discovery_control


def perform_dge(expr_df):
	"""
	Performs differential gene expression analysis between matched tumor and normal samples.
	Returns a DataFrame with DGE results.
	"""
	
	print("Running DGE analysis")
	# Separate primary tumor (marked by 01 string) and matched normal (marked by 11 string) samples
	# Ignoring Recurrent and Metastatic Tumors for this exercise
	tumor_samples = [s for s in expr_df.columns if s.split("-")[3].startswith('01')]
	normal_samples = [s for s in expr_df.columns if s.split("-")[3].startswith('11')]	

	# Find and sort matched pairs
	patient_ids = set(s.split('-')[2] for s in tumor_samples) & set(s.split('-')[2] for s in normal_samples)
	matched_tumor = sorted(s for s in tumor_samples if s.split("-")[2] in patient_ids)
	matched_normal = sorted(s for s in normal_samples if s.split("-")[2] in patient_ids)

	expr_tumor = expr_df[matched_tumor]
	expr_normal = expr_df[matched_normal]

	# Perform Welch's t-test for DGE
	results = []
	for gene in expr_df.index:
		t_stat, p_val = ttest_ind(expr_tumor.loc[gene], expr_normal.loc[gene], equal_var=False, nan_policy='omit')
		log2fc = expr_tumor.loc[gene].mean() - expr_normal.loc[gene].mean()
		results.append({'gene': gene, 'log2fc': log2fc, 'p_value': p_val})

	dge_results_df = pd.DataFrame(results).dropna().set_index('gene')	

	# Adjust for multiple testing
	dge_results_df['fdr'] = false_discovery_control(dge_results_df['p_value'])	

	return dge_results_df, matched_tumor, matched_normal

def prepare_survival_data(expr_df, clinical_df):
	"""
	Align expression and clinical data for survival analysis.
	Returns data frame for plotting.
	"""
	print("Preparing data for survival analysis...")
