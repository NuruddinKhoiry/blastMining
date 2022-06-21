#!/usr/bin/env python3

import pandas as pd

def summary_df(DF, sample_name):
    DF['Taxa'] = DF[['Kingdom','Phylum','Class','Order','Family','Genus', 'Species']].apply(';'.join, axis=1)
    DT = DF.drop(columns=['Kingdom','Phylum','Class','Order','Family','Genus', 'Species'])
    summary = pd.DataFrame(DT.groupby(['Taxa'])['Taxa'].count())
    summary.columns = [str(sample_name)]
    summary_DF = summary.reset_index(level=0)
    return(summary_DF)