import dash_bio as dashbio
import matplotlib.pyplot as plt

def create_interactive_heatmap(expr_df, dge_results_df, tumor_samples, normal_samples, top_n=50):
	"""Generates and saves an interactive heatmap of top diff expr genes."""
	print(f"Generating heat map for {top_n} genes...")
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

	fig.write_html('../output/diffexp_heatmap.html')
	print("Interactive heatmap saved to 'output/diffexp_heatmap.html'")
