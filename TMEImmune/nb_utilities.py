import numpy as np 
import pandas as pd
import os, json
import gseapy as gp
import warnings
from collections import defaultdict
from sklearn.preprocessing import StandardScaler
import importlib.resources as pkg_resources
#from TMEImmune import data_processing
from data_processing import load_data


class nb_pathway():
	def __init__(self):
		#self.coef = pd.read_csv("data/coef_best_netbio.csv")
        #response = requests.get(url + "gene_sets.json")
		#self.gene_set = data_processing.load_data("gene_sets_full.json")
		self.gene_set = load_data("gene_sets_full.json")
    
	# def get_nb_coef(self):
	# 	return self.coef
    
	def reactome_geneset(self):
		return self.gene_set
	

class netbio_data:
	def __init__(self, gene, clin, response, ssgsea = None):
		self.gene = gene
		self.clin = clin[clin[response].isin(["R", "NR"])]
		self._ssgsea = ssgsea
		self.response = response
		self.common_id = list(set(gene.columns) & set(clin.index))

	def get_gene(self, gene_id):
		gene_columns = [gene_id] + self.common_id
		return self.gene[gene_columns]
	
	def get_clin(self):
		clin = self.clin.loc[self.common_id,:]
		clin_resp = clin[self.response].apply(lambda x: 1 if x == "R" else 0)
		return clin_resp
	
	def get_ssgsea(self):
		ssgsea_result = self._ssgsea
		if self._ssgsea is None:
			ssgsea_result = get_ssgsea(self.gene)
		ssgsea_col = ['pathway'] + self.common_id
		ssgsea_result = ssgsea_result[ssgsea_col]
		return ssgsea_result


## pathway expression and immunotherapy response
def parse_reactomeExpression_and_immunotherapyResponse(dataset, Prat_cancer_type='MELANOMA'):

	# edf = parse_gene_expression(dataset)
	# epdf = parse_reactome_expression(dataset)
	# pdf = parse_immunotherapy_response(dataset)
	
    # edf = data_processing.load_data('Gide/expression_mRNA.norm3.txt')
    # epdf = data_processing.load_data('Gide/pathway_expression_ssgsea.txt')
    # pdf = data_processing.load_data('Gide/patient_df.txt')

	edf = load_data('Gide/expression_mRNA.norm3.txt')
	epdf = load_data('Gide/pathway_expression_ssgsea.txt')
	pdf = load_data('Gide/patient_df.txt')
	
	# features, labels
	exp_dic, responses = defaultdict(list), []
	e_samples = []
	for sample, response in zip(pdf['Patient'].tolist(), pdf['Response'].tolist()):
		# labels
		binary_response = response
		if response == 'NR':
			binary_response = 0
		if response == 'R':
			binary_response = 1
		# features
		for e_sample in epdf.columns:
			tmp = []

			if (sample in e_sample) or (e_sample in sample):
				e_samples.append(e_sample)
				responses.append(binary_response)	

	edf = pd.DataFrame(data=edf, columns=np.append(['gene_id'], e_samples))
	edf = edf.rename(columns = {"gene_id": "genes"})
	epdf = pd.DataFrame(data=epdf, columns=np.append(['pathway'], e_samples))
	responses = np.array(responses)
	return e_samples, edf, epdf, responses


# pathway expression
# def parse_reactome_expression(dataset):
# 	epdf = pd.DataFrame()
# 	# directory
# 	data_dir = 'Gide/pathway_expression_ssgsea.txt'
# 	epdf = data_processing.load_data(data_dir)

# 	# data cleanup
# 	if len(epdf)>0:
# 		#epdf = epdf.rename(columns={'testType':'pathway'})
# 		#epdf = epdf[epdf['pathway'].str.contains('REACTOME_')]
# 		epdf = epdf.dropna()
# 	return epdf



## immunotherapy response
# def parse_immunotherapy_response(dataset):
# 	'''
# 	Input
# 	dataset : 'IMvigor210', 'Liu', 'Riaz', 'Gide', 'Prat', 'Kim', 'Auslander'
# 	'''
# 	# directory
# 	data_dir = 'Gide/patient_df.txt'
# 	pdf = data_processing.load_data(data_dir)
# 	# if ('Liu' in fldr) or ('Gide' in fldr):
# 	# 	pdf['Response'] = pdf['Response'].astype(str)
# 	# 	pdf = pdf.loc[pdf['Response'].isin(['PD', 'PR', 'CR', 'SD']),:]
# 	# 	tmp_response = []
# 	# 	for r in pdf['Response'].tolist():
# 	# 		if r in ['CR', 'PR']:
# 	# 			tmp_response.append('R')
# 	# 		if r in ['SD', 'PD']:
# 	# 			tmp_response.append('NR')
# 	# 	pdf['Response'] = tmp_response
# 	return pdf


## gene expression 
# def parse_gene_expression(dataset):
# 	# directory
# 	data_dir = 'Gide/expression_mRNA.norm3.txt'
# 	edf = data_processing.load_data(data_dir)
# 	edf = edf.dropna()
# 	return edf

# def reactome_genes():
# 	output = defaultdict(list)
# 	output_list = []
# 	f = open('data/c2.all.v7.2.symbols.gmt','r')
# 	lines = f.readlines()
# 	for line in lines:
# 		line = line.strip().split('\t')
# 		if 'REACTOME' in line[0]:
# 			reactome = line[0]
# 			output_list.append(reactome)
# 			for i in range(2, len(line)):
# 				gene = line[i]
# 				output[reactome].append(gene)
# 	f.close()
# 	return output

def expression_StandardScaler(exp_df):
	'''
	Input : expression dataframe
	'''
	col1 = exp_df.columns[0]
	tmp = StandardScaler().fit_transform(exp_df.T.values[1:])
	new_tmp = defaultdict(list)
	new_tmp[col1] = exp_df[col1].tolist()
	for s_idx, sample in enumerate(exp_df.columns[1:]):
		new_tmp[sample] = tmp[s_idx]
	output = pd.DataFrame(data=new_tmp, columns=exp_df.columns)
	return output



def get_ssgsea(df):

    # get netbio reactome pathways and corresponding genesets
	gene_set_dict = nb_pathway().reactome_geneset()
	# Remove pathways where the gene list is empty
	gene_set_dict = {pathway: genes for pathway, genes in gene_set_dict.items() if genes}

    # perform ssgsea
	ssgsea_results = gp.ssgsea(
	data=df,
    gene_sets=gene_set_dict,  # Path to the gene set file
    min_size=0,
    outdir=None,  # Avoid file output
    verbose=True,
	sample_norm_method="rank"
    )
	nes_pivot = ssgsea_results.res2d.pivot(index='Term', columns='Name', values='NES')

	nes_pivot = nes_pivot.replace(r'^\s*$', np.nan, regex=True)
	nes_pivot = nes_pivot.replace("nan", np.nan)
	nes_pivot = nes_pivot.dropna(how='all')
	nes_pivot = nes_pivot.reset_index().rename(columns={'Term': 'pathway'}) 

	return nes_pivot