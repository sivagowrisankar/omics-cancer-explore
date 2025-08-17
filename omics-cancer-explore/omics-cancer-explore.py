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
	clinical_df = pd.read_csv('../data/TCGA_MetaData.tsv', sep='\t', index_col=0)

	# Step 1: DGE and Heatmap
	print("\n--- Running Step 1: DGE and Heatmap--")
	dge_results, matched_tumor, matched_normal = analysis.perform_dge(expr_df)	
	visualization.create_interactive_heatmap(expr_df, dge_results, matched_tumor, matched_normal, 100)

if __name__ == '__main__':
	main()
