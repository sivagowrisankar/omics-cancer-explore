import os
from omics_cancer_explore import analysis, visualization

def test_dge_to_heatmap_pipeline(mock_expression_data, tmp_path):
	"""Integration test: DGE analysis through heatmap generation."""
	os.chdir(tmp_path) # Use a temporary directory for outputs
	outputdir = 'output'
	os.makedirs(outputdir, exist_ok = True)
    
	dge_results, matched_tumor, matched_normal = analysis.perform_dge(mock_expression_data)
	visualization.create_interactive_heatmap(mock_expression_data, dge_results, matched_tumor, matched_normal, outputdir)
    
	assert os.path.exists('output/diffexp_heatmap.html')
