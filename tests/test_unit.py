from omics_cancer_explore import analysis

def test_perform_dge(mock_expression_data):
	"""Unit test for the DGE analysis function."""
	dge_results, _, _ = analysis.perform_dge(mock_expression_data)
	assert 'fdr' in dge_results.columns
	assert dge_results.loc['GENE_5']['log2fc'] > 0
	assert dge_results.sort_values('fdr').index[0] == 'GENE_5'

def test_prepare_survival_data(mock_expression_data, mock_clinical_data):
	"""Unit test for survival data preparation"""
	survival_df = analysis.prepare_survival_data(mock_expression_data, mock_clinical_data)
	assert 'time_to_event' in survival_df.columns
	assert survival_df.loc['TCGA-AA-0001']['time_to_event'] == 800 
	assert survival_df.loc['TCGA-AA-0000']['time_to_event'] == 1200
