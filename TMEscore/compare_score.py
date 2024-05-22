import pandas as pd
import statsmodels.api as sm
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
import io



def coxph(method, df):

    # p value
    model = sm.PHReg(df['time'], df[method], status=df['delta'])
    result = model.fit()
    p_value = result.pvalues[0]

    # c index

    estimator = CoxPHSurvivalAnalysis()
    estimator.fit(df[[method]], df[['delta', 'time']].to_records(index=False))
    c_index_value = concordance_index_censored(df['delta'], df['time'], estimator.predict(df[[method]]))[0]
    return p_value, c_index_value


def compare(score_df, survival_df):
    score_df.columns[0] = "id"
    survival_df.columns = ["id", "time", "vital"] # vital boolean
    total = score_df.merge(survival_df, on="id")

    total = total.dropna(subset = ['time','vital'])
    total = total[total.time > 0]

    total['time'] = total['time'].astype(float)
    total['delta'] = total['vital'] == 'dead'

    #method = ['ESTIMATE_stromal', 'ESTIMATE_immune', 'istme_stromal', 'istme_immune', 'sia_score']
    method = score_df.columns[1:]
    output = {}
    for i in method:
        p, c = coxph(i, total)
        output[i] = [p, c]

    output = pd.DataFrame(output)
    output.index = ['p-value', 'c-index']
    output = output.round(decimals=4)

    return output


def compare_tcga(cancer_type):
    contents = """
cancer,score,c_index,p_value
OV,ESTIMATE_stromal,0.5344,0.1784
LUAD,ESTIMATE_stromal,0.4843,0.9341
LUSC,ESTIMATE_stromal,0.5445,0.0481
BLCA,ESTIMATE_stromal,0.5461,0.2639
BRCA,ESTIMATE_stromal,0.5164,0.3605
KIRC,ESTIMATE_stromal,0.5153,0.6446
COADREAD,ESTIMATE_stromal,0.5685,0.0294
COAD,ESTIMATE_stromal,0.5528,0.0722
HNSC,ESTIMATE_stromal,0.5279,0.0837
SKCM,ESTIMATE_stromal,0.6053,0.0103
OV,ESTIMATE_immune,0.5107,0.3222
LUAD,ESTIMATE_immune,0.5402,0.5603
LUSC,ESTIMATE_immune,0.5216,0.3451
BLCA,ESTIMATE_immune,0.5115,0.35
BRCA,ESTIMATE_immune,0.5408,0.1299
KIRC,ESTIMATE_immune,0.528,0.1118
COADREAD,ESTIMATE_immune,0.5666,0.1392
COAD,ESTIMATE_immune,0.5641,0.187
HNSC,ESTIMATE_immune,0.5634,0.0533
SKCM,ESTIMATE_immune,0.7038,0
OV,istme_stromal,0.54,0.4497
LUAD,istme_stromal,0.5267,0.1709
LUSC,istme_stromal,0.5045,0.8856
BLCA,istme_stromal,0.5409,0.4373
BRCA,istme_stromal,0.6028,0.3227
KIRC,istme_stromal,0.6146,0
COADREAD,istme_stromal,0.5661,0.0393
COAD,istme_stromal,0.5613,0.0666
HNSC,istme_stromal,0.5054,0.2583
SKCM,istme_stromal,0.6646,1.00E-04
OV,istme_immune,0.507,0.5685
LUAD,istme_immune,0.5432,0.5321
LUSC,istme_immune,0.5074,0.6698
BLCA,istme_immune,0.5177,0.3059
BRCA,istme_immune,0.5486,0.145
KIRC,istme_immune,0.5018,0.8564
COADREAD,istme_immune,0.5649,0.2274
COAD,istme_immune,0.5639,0.2526
HNSC,istme_immune,0.5542,0.0946
SKCM,istme_immune,0.7226,0
OV,sia_score,0.5269,0.6665
LUAD,sia_score,0.5115,0.399
LUSC,sia_score,0.4372,0.4057
BLCA,sia_score,0.5054,0.9495
BRCA,sia_score,0.6059,0.1178
KIRC,sia_score,0.4867,0.1982
COADREAD,sia_score,0.5926,0.0277
COAD,sia_score,0.5912,0.0781
HNSC,sia_score,0.5634,0.6349
SKCM,sia_score,0.6459,0.52
    """
    tcga_cox = pd.read_csv(io.StringIO(contents))
    tcga_subset = tcga_cox[tcga_cox['cancer'] == cancer_type]
    if tcga_subset[tcga_subset.p_value <= 0.1].shape[0] > 0:
        subset = tcga_subset[tcga_subset.p_value <= 0.1]
    else:
        subset = tcga_subset
    max_index = subset['c_index'].idxmax()
    op_score = subset.loc[max_index][1]

    return tcga_subset, op_score

