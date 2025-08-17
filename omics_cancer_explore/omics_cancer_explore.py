import pandas as pd
import os
import argparse

import analysis,visualization

def run_analysis(outputdir, survgene):
	"""Main function to run the entire analysis pipeline."""

	#Ensure output directory exists
	os.makedirs('output', exist_ok = True)

	# Load Data
	print("Loading Data...")
	expr_df = pd.read_csv('data/TCGA.COADREAD.sampleMap_HiSeqV2.tsv',sep='\t', index_col=0)
	clinical_df = pd.read_csv('data/clinical.tsv', sep='\t')
	clinical_df = clinical_df.drop_duplicates(subset=['cases.submitter_id'], keep='first')
	clinical_df = clinical_df.set_index('cases.submitter_id')

	# Step 1: DGE and Heatmap
	print("\n--- Running Step 1: DGE and Heatmap--")
	dge_results, matched_tumor, matched_normal = analysis.perform_dge(expr_df)
	visualization.create_interactive_heatmap(expr_df, dge_results, matched_tumor, matched_normal, outputdir, 100)

	# Step 2: Survival Analysis
	print("\n--- Running Step 2: Survival Analysis ---")
	survival_data = analysis.prepare_survival_data(expr_df, clinical_df)
	visualization.create_KM_plot(survival_data, outputdir, gene_symbol=survgene)

def get_args():
	parser = argparse.ArgumentParser(description='Analyze TCGA data.')
	parser.add_argument('--outputdir',
				type=str,
				help='Output path to heatmap and survival curve',
				required=True)
	parser.add_argument('--survgene',
				type=str,
				help='Gene for which differential survival plots need to be generated',
				required=True)
	return parser

def main():
	parser = get_args()
	args = parser.parse_args()
	run_analysis(args.outputdir, args.survgene)

if __name__ == '__main__':
	main()
