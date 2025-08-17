import pytest
import pandas as pd
import numpy as np

@pytest.fixture(scope="module")
def mock_expression_data():
	"""Creates a mock gene expression DataFrame."""
	genes = [f'GENE_{i}' for i in range(100)]
	samples = [f'TCGA-AA-{i:04d}-01A' for i in range(5)] + [f'TCGA-AA-{i:04d}-11A' for i in range(5)]
	data = np.random.rand(100,10)
	df = pd.DataFrame(data, index=genes, columns=samples)

	#Make two diff exp genes
	df.iloc[5, :5] = 5.0 #High in Tumor (01A)
	df.iloc[5, 5:] = 0.5 #Low in Normal (11A)

	return df

@pytest.fixture(scope="module")
def mock_clinical_data():
	"""Creates a mock clinical data DataFrame."""
	patients = [f'TCGA-AA-{i:04d}' for i in range(5)]
	data = {
		'demographic.vital_status' : ['Alive', 'Dead', 'Alive', 'Dead', 'Alive'],
		'demographic.days_to_death' : [np.nan, 800, np.nan, 450, np.nan],
		'diagnoses.days_to_last_follow_up': [1200, np.nan, 950, np.nan, 700]
	}
	
	return pd.DataFrame(data, index=patients)
