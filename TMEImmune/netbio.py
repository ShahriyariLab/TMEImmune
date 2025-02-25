import numpy as np 
import pandas as pd
import time, warnings
from collections import defaultdict
import scipy.stats as stat
from sklearn.model_selection import cross_val_score, KFold, train_test_split, GridSearchCV, StratifiedKFold
import networkx as nx
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_curve, auc, accuracy_score, f1_score, precision_score, recall_score, precision_recall_curve
from sklearn.linear_model import LogisticRegression
from scipy.special import logit
#from TMEImmune import data_processing
from data_processing import load_data

warnings.filterwarnings('ignore')
warnings.filterwarnings(action='ignore', category=DeprecationWarning)
warnings.filterwarnings(action='ignore', category=FutureWarning)


exec(open('nb_utilities.py').read())

#write a class of datasets, option to perform ssgsea

## Run cross study predictions
def get_netbio(test_gene, test_clin, test_clinid, test_ssgsea = None, 
			   cohort_targets = {'PD1':['Gide']}, train_gene = None, 
			   train_geneid = "gene_id", train_clin = None, train_clinid = None, train_ssgsea = None, 
			   nGene = 200, qval = 0.01, penalty = "l2"):

	if train_gene is None:
		train_dataset = "Gide"
		train_geneid = 'genes'
		target = "PD1"
		train_samples, train_edf, train_epdf, train_responses = parse_reactomeExpression_and_immunotherapyResponse(train_dataset)

	else:
		target, train_dataset = cohort_targets.items
		train_data = netbio_data(train_gene, train_clin, train_clinid, train_ssgsea)
		train_edf, train_epdf, train_responses = train_data.get_gene(train_geneid), train_data.get_ssgsea(), train_data.get_clin()


	#reactome = data_processing.load_data("c2.all.v7.2.symbols.gmt")
	reactome = load_data("c2.all.v7.2.symbols.gmt")

	# if 'Prat' in test_dataset:
	# 	#test_samples, test_edf, test_epdf, test_responses = parse_reactomeExpression_and_immunotherapyResponse(test_dataset)
	# 	test_samples, test_edf, _, test_responses = parse_reactomeExpression_and_immunotherapyResponse(test_dataset.split('_')[0], Prat_cancer_type=test_dataset.split('_')[1])
	# 	test_epdf = get_ssgsea(test_edf)
	# 	test_pathways = ['pathway'] + test_samples
	# 	test_epdf = test_epdf[test_pathways]
	# 	test_genes = ['genes'] + test_samples
	# 	test_edf = test_edf[test_genes]
	# 	#test_epdf.to_csv("test_epdf.txt", sep = "\t", index = None)
	# else:


	test_data = netbio_data(test_gene, test_clin, test_clinid, test_ssgsea)
	test_geneid = test_gene.columns[0]
	test_edf, test_epdf, test_responses = test_data.get_gene(test_geneid), test_data.get_ssgsea(), test_data.get_clin()

	# print('\n\n#----------------------------------------')
	# print('training: %s %s'%(train_dataset, time.ctime()))
	# print('test data --> responder : %s / non-responder : %s'%(list(test_responses).count(1), list(test_responses).count(0)))

	### data cleanup: match genes and pathways between cohorts
	#print(train_edf.shape, train_epdf.shape, test_edf.shape, test_epdf.shape)
	common_genes, common_pathways = list(set(train_edf[train_geneid].tolist()) & set(test_edf[test_geneid].tolist())), list(set(train_epdf['pathway'].tolist()) & set(test_epdf['pathway'].tolist()))
	train_edf = train_edf.loc[train_edf[train_geneid].isin(common_genes),:].sort_values(by=train_geneid)
	train_epdf = train_epdf.loc[train_epdf['pathway'].isin(common_pathways),:].sort_values(by='pathway')
	test_edf = test_edf.loc[test_edf[test_geneid].isin(common_genes),:].sort_values(by=test_geneid)
	test_epdf = test_epdf.loc[test_epdf['pathway'].isin(common_pathways),:].sort_values(by = 'pathway')
	#print(train_edf.shape, train_epdf.shape, test_edf.shape, test_epdf.shape)


	### data cleanup: expression standardization
	train_edf = expression_StandardScaler(train_edf)
	train_epdf = expression_StandardScaler(train_epdf)
	test_edf = expression_StandardScaler(test_edf)
	test_epdf = expression_StandardScaler(test_epdf)

	
	### network proximal pathways
	# gene expansion by network propagation results
# biomarker genes 
	#bio_df = data_processing.load_data("Marker_summary.txt")

	biomarker_dir = 'nb_biomarker/'
	bdf = load_data('%s/%s.txt'%(biomarker_dir, target))
	#bdf = data_processing.load_data('%s/%s.txt'%(biomarker_dir, target))
	bdf = bdf.dropna(subset=['gene_id'])
	b_genes = []
	for idx, gene in enumerate(bdf.sort_values(by=['propagate_score'], ascending=False)['gene_id'].tolist()):
		if gene in train_edf[train_geneid].tolist():
			if not gene in b_genes:
				b_genes.append(gene)
			if len(set(b_genes)) >= nGene:
				break
	# LCC function enrichment
	tmp_hypergeom = defaultdict(list)
	pvalues, qvalues = [], []
	for pw in list(reactome.keys()):
		pw_genes = list(set(reactome[pw]) & set(train_edf[train_geneid].tolist()))
		M = len(train_edf[train_geneid].tolist())
		n = len(pw_genes)
		N = len(set(b_genes))
		k = len(set(pw_genes) & set(b_genes))
		p = stat.hypergeom.sf(k-1, M, n, N)
		tmp_hypergeom['pw'].append(pw)
		tmp_hypergeom['p'].append(p)
		pvalues.append(p)
	_, qvalues, _, _ = multipletests(pvalues)
	tmp_hypergeom['q'] = qvalues
	tmp_hypergeom = pd.DataFrame(tmp_hypergeom).sort_values(by=['q'])
	proximal_pathways = tmp_hypergeom.loc[tmp_hypergeom['q']<=qval,:]['pw'].tolist() ## proximal_pathways


	### Train / Test dataset merging
	train_dic = {}
	test_dic = {}
	# 1. NetBio
	train_dic['NetBio'] = train_epdf.loc[train_epdf['pathway'].isin(proximal_pathways),:]
	test_dic['NetBio'] = test_epdf.loc[test_epdf['pathway'].isin(proximal_pathways),:]
	# # 2. controls (other biomarkers)
	# for test_type, genes in zip(bio_df['Name'].tolist(), bio_df['Gene_list'].tolist()):
	# 	genes = genes.split(':')
	# 	train_dic[test_type] = train_edf.loc[train_edf['genes'].isin(genes),:]
	# 	test_dic[test_type] = test_edf.loc[test_edf['genes'].isin(genes),:]
	# 	ptws = train_dic['NetBio']['pathway']
	
	
	### Measure Prediction Performances
	#print('\tML training & predicting, %s'%time.ctime())
	#for test_type in np.append(['NetBio'], bio_df['Name'].tolist()):
	X_train, X_test = train_dic['NetBio'].T.values[1:], test_dic['NetBio'].T.values[1:]
	y_train, y_test = train_responses, test_responses
	# if X_train.shape[1] * X_train.shape[0] * len(y_train) * len(y_test) == 0:
	# 	continue

	# make predictions
	model = LogisticRegression()
	if penalty == 'l2':
		param_grid = {'penalty':['l2'], 'max_iter':[1000], 'solver':['lbfgs'], 'C':np.arange(0.1, 1, 0.1), 'class_weight':['balanced'] }
	if penalty == 'none':
		param_grid = {'penalty':['none'], 'max_iter':[1000], 'class_weight':['balanced'] }
	cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
	gcv = GridSearchCV(model, param_grid=param_grid, cv=cv, scoring='roc_auc', 
			n_jobs=5).fit(X_train, y_train) #cv=5
	# model = gcv.best_estimator_
	# best_coef = model.coef_
	# best_intercept = model.intercept_
	pred_status = gcv.best_estimator_.predict(X_test)
	pred_proba = gcv.best_estimator_.predict_proba(X_test)[:,1]
	logit_pred = logit(np.clip(pred_proba, 1e-6, 1 - 1e-6))

		#### measure performance
		# AUC (prediction probability)
	# fpr_proba, tpr_proba, _ = roc_curve(y_test, pred_proba, pos_label=1)
	# AUC_proba = auc(fpr_proba, tpr_proba)
	# 	# AUPRC
	# precision, recall, _ = precision_recall_curve(y_test, pred_proba, pos_label=1)
	# AUPRC = auc(recall, precision)
	# expected_AUPRC = list(y_test).count(1)/len(y_test)

		# final results
	#if test_type == 'NetBio':
	# print('\n\t%s, %s, train: %s, test: %s, %s'%("NetBio", "logistic regression", train_dataset, test_dataset, time.ctime()))

	# for metric, score, expected_score in zip(['AUC_proba', 'AUPRC'], [AUC_proba, AUPRC], [0.5, expected_AUPRC]):
	# 	print('\t%s, %s -- %s : %s (random expectation=%s)'%("NetBio", "logistic regression", metric, score, expected_score))


	#print('Finished, %s'%time.ctime())

	return logit_pred
