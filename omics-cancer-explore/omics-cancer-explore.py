import pandas as pd
import os

import analysis,visualization

def main():
	"""Main function to run the entire analysis pipeline."""
	#Ensure output directory exists
	os.makedirs('output', exist_ok = True)


	# Load Data 
	print("Loading Data...")
	expr_df = pd.read_csv('../data/TCGA.COADREAD.sampleMap_HiSeqV2.tsv',sep='\t', index_col=0)
	clinical_df = pd.read_csv('../data/clinical.tsv', sep='\t')
	clinical_df = clinical_df.drop_duplicates(subset=['cases.submitter_id'], keep='first')
	clinical_df = clinical_df.set_index('cases.submitter_id')

	# Step 1: DGE and Heatmap
	print("\n--- Running Step 1: DGE and Heatmap--")
	dge_results, matched_tumor, matched_normal = analysis.perform_dge(expr_df)	
	visualization.create_interactive_heatmap(expr_df, dge_results, matched_tumor, matched_normal, 100)

	# Step 2: Survival Analysis
	print("\n--- Running Step 2: Survival Analysis ---")
	survival_data = analysis.prepare_survival_data(expr_df, clinical_df)
	visualization.create_KM_plot(survival_data, gene_symbol='PMS2')

if __name__ == '__main__':
	main()
