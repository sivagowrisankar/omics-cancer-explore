import pandas as pd
import os,sys
import argparse
import analysis,visualization
import logging
from datetime import datetime

def run_analysis(outputdir, survgene):
	"""Main function to run the entire analysis pipeline."""

	# Load Data
	logging.info("Loading Data...")

	# Pre-run checks
	logging.info("Performing pre-run checks...")
	data_dir = 'data'
	expr_file = os.path.join(data_dir, 'TCGA.COADREAD.sampleMap_HiSeqV2.tsv')
	clinical_file = os.path.join(data_dir, 'clinical.tsv')
	for path in [data_dir, expr_file, clinical_file]:
		if not os.path.exists(path):
			logging.error(f"Error: Required file or directory not found: {path}", file=sys.stderr)
			sys.exit(1)

	try:
		# Load data
		logging.info("Logging data...")	
		expr_df = pd.read_csv(expr_file,sep='\t', index_col=0)
		clinical_df = pd.read_csv(clinical_file, sep='\t')
		clinical_df = clinical_df.drop_duplicates(subset=['cases.submitter_id'], keep='first')
		clinical_df = clinical_df.set_index('cases.submitter_id')

		# Step 1: DGE and Heatmap
		logging.info("--- Running Step 1: DGE and Heatmap--")
		dge_results, matched_tumor, matched_normal = analysis.perform_dge(expr_df)
		visualization.create_interactive_heatmap(expr_df, dge_results, matched_tumor, matched_normal, outputdir, 100)

		# Step 2: Survival Analysis
		logging.info("--- Running Step 2: Survival Analysis ---")
		survival_data = analysis.prepare_survival_data(expr_df, clinical_df)
		visualization.create_KM_plot(survival_data, outputdir, gene_symbol=survgene)
	
	except (FileNotFoundError, KeyError, ValueError) as e:
		logging.error(f"An error occurred during analysis: {e}", exc_info=True)
		sys.exit(1)

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
	date_time = datetime.now().strftime("%Y%m%d_%H%M%S")
	log_file_name = f"analysis_{date_time}.log"
	log_file_path = os.path.join(args.outputdir, log_file_name)
	
	#Ensure output directory exists
	os.makedirs('output', exist_ok = True)

	# Configure root logger
	logging.basicConfig(
		level=logging.INFO,
		format="%(asctime)s [%(levelname)s] %(message)s",
		handlers=[
			logging.FileHandler(log_file_path),
			logging.StreamHandler(sys.stdout) # Also log to console
		]
	)

	logging.info("Logger configured. Starting the analysis pipeline.")
	logging.info(f"Arguments: output_dir={args.outputdir}, gene={args.survgene}")
	run_analysis(args.outputdir, args.survgene)

if __name__ == '__main__':
	main()
