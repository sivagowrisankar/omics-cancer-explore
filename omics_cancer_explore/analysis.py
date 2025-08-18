import pandas as pd
from scipy.stats import ttest_ind, false_discovery_control
import numpy as np
import logging

# Get the logger instance
logger = logging.getLogger(__name__)

def perform_dge(expr_df):
	"""
	Performs differential gene expression analysis between matched tumor and normal samples.
	Returns a DataFrame with DGE results.
	"""
	
	logger.info("Running DGE analysis")
	# Separate primary tumor (marked by 01 string) and matched normal (marked by 11 string) samples
	# Ignoring Recurrent and Metastatic Tumors for this exercise
	tumor_samples = [s for s in expr_df.columns if s.split("-")[3].startswith('01')]
	normal_samples = [s for s in expr_df.columns if s.split("-")[3].startswith('11')]	

	# Find and sort matched pairs
	patient_ids = set(s.split('-')[2] for s in tumor_samples) & set(s.split('-')[2] for s in normal_samples)
	matched_tumor = sorted(s for s in tumor_samples if s.split("-")[2] in patient_ids)
	matched_normal = sorted(s for s in normal_samples if s.split("-")[2] in patient_ids)

	if not matched_tumor or not matched_normal:
		raise ValueError("Error: No matched tumor-normal pairs found in the dataset. Cannot perform differential expression")

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
	logger.info("Preparing data for survival analysis...")
	required_cols = ['demographic.vital_status','demographic.days_to_death','diagnoses.days_to_last_follow_up']
	
	if not all(col in clinical_df.columns for col in required_cols):
		raise KeyError(f"Error: Clinical data is missing one of the required columns: {required_cols}")

	tumor_samples = [s for s in expr_df.columns if s.split('-')[3].startswith('01')]
	expr_tumor = expr_df[tumor_samples].T
	expr_tumor.index = expr_tumor.index.str.slice(0,12)

	aligned_df = expr_tumor.join(clinical_df, how='inner')

	# Create survival points
	aligned_df['deceased'] = (aligned_df['demographic.vital_status'] == 'Dead').astype(int)
	aligned_df[['demographic.days_to_death', 'diagnoses.days_to_last_follow_up']] = (
		aligned_df[['demographic.days_to_death', 'diagnoses.days_to_last_follow_up']]
		.replace('\'--', np.nan)
	)

	aligned_df['demographic.days_to_death'] = pd.to_numeric(
		aligned_df['demographic.days_to_death'], errors="coerce"
	)
	aligned_df['diagnoses.days_to_last_follow_up'] = pd.to_numeric(
		aligned_df['diagnoses.days_to_last_follow_up'], errors="coerce"
	)

	aligned_df['time_to_event'] = aligned_df.apply(
        	lambda row: row['demographic.days_to_death'] if row['deceased'] == 1 else row['diagnoses.days_to_last_follow_up'],
        	axis=1
	)
	aligned_df = aligned_df.dropna(subset=['time_to_event'])

	aligned_df['time_to_event'] = aligned_df['time_to_event'].astype(float)

	return aligned_df
