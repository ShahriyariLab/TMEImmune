import pandas as pd
from TMEImmune import data_processing, TME_score, optimal
import os


# run this line in the terminal to test the container
# docker run --rm -v $(pwd):/app tmeimmune_image python /app/docker_test.py

def main():
    # Read clinical data
    clin = pd.read_csv("example_clin.csv", index_col=0)
    
    # Process gene expression data using normalization
    df = data_processing.normalization(
        path="example_gene.csv", 
        method='TMM', 
        batch=clin, 
        batch_col="CANCER"
    )
    
    # Compute TME scores using multiple methods
    score = TME_score.get_score(
        df, 
        method=['ESTIMATE', 'ISTME', 'NetBio', 'SIA'], 
        clin=clin, 
        test_clinid="response"
    )
    
    # Evaluate performance with the given metrics and score names
    outcome = optimal.get_performance(
        score, 
        metric=['ICI', 'survival'], 
        score_name=['EST_stromal', 'EST_immune', 'IS_immune', 'IS_stromal', 'NetBio', 'SIA'], 
        ICI_col='response', 
        surv_col=['time', 'delta'], 
        df_clin=clin
    )
    
    km_path = "output/test_km.pdf"
    roc_path = "output/test_roc.pdf"
    outcome[0][0].savefig(roc_path) 
    outcome[1][0].savefig(km_path) 

if __name__ == "__main__":
    main()
