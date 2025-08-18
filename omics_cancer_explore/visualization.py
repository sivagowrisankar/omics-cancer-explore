import dash_bio as dashbio
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import logging

#Get a logger instance
logger = logging.getLogger(__name__)

def create_interactive_heatmap(expr_df, dge_results_df, tumor_samples, normal_samples, outputdir, top_n=50):
	"""Generates and saves an interactive heatmap of top diff expr genes."""
	logger.info(f"Generating heat map for {top_n} genes...")
	top_genes = dge_results_df.sort_values('fdr').head(top_n)	
	heatmap_data = expr_df.loc[top_genes.index, tumor_samples+normal_samples].copy()
	heatmap_data = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(heatmap_data.std(axis=1), axis=0)

	fig = dashbio.Clustergram(
		data=heatmap_data,
		column_labels=list(heatmap_data.columns.values),
		row_labels=list(heatmap_data.index),
		height=1800,
		width=1400,
		line_width=2
	)
	
	fig.update_layout(
		title=f'Top {top_n} Differentially expressed genes (Tumor vs. Matched Normal)',
	)

	outfilename = f"{outputdir}/diffexp_heatmap.html"

	fig.write_html(outfilename)
	logger.info(f"Interactive heatmap saved to '{outfilename}'")

def create_KM_plot(survival_df, outputdir, gene_symbol='ETV4'):
	"""Generates and saves a KM survival plot."""
	logger.info(f"Plotting survival curve for gene: {gene_symbol}...")

	if gene_symbol not in survival_df.columns:
		raise KeyError(f"Error: Gene '{gene_symbol}' not found in the expression data")

	# Stratify patients by median gene expression
	median_expr = survival_df[gene_symbol].median()
	survival_df['expression_group'] = (survival_df[gene_symbol] > median_expr).map({True: 'High', False: 'Low'})

	kmf = KaplanMeierFitter()
	fig, ax = plt.subplots(figsize=(8, 6))
	
	high_group = survival_df[survival_df['expression_group'] == 'High']
	low_group = survival_df[survival_df['expression_group'] == 'Low']
	
	kmf.fit(high_group['time_to_event'], high_group['deceased'], label='High Expression')
	kmf.plot_survival_function(ax=ax)

	kmf.fit(low_group['time_to_event'], low_group['deceased'], label='Low Expression')
	kmf.plot_survival_function(ax=ax)

	# Perform log-rank test
	results = logrank_test(
		high_group['time_to_event'], low_group['time_to_event'],
		event_observed_A=high_group['deceased'], event_observed_B=low_group['deceased']
	)

	ax.set_title(f'Survival Curves by {gene_symbol} Expression (p-value: {results.p_value:.3g})')
	ax.set_xlabel('Days')
	ax.set_ylabel('Overall Survival Probability')
	plt.tight_layout()
	outfilename = f"{outputdir}/{gene_symbol}_surv_curve.png"
	plt.savefig(outfilename)

	logger.info(f"Survival curve saved to '{outfilename}'")
