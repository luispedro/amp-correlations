from jug import TaskGenerator, iteratetask
from preproc import filter_columns

filter_columns = TaskGenerator(filter_columns)

@TaskGenerator
def run_corrs(amp_name, motus_name, mode):
    import pandas as pd
    from scipy import stats
    import corr
    data = pd.read_table(amp_name, index_col=0)
    motus = pd.read_table(motus_name, index_col=0)

    mv = motus.T.values.copy()

    if mode == 'pearsonr':
        return pd.DataFrame(corr.pearsonr_pcorr(data.values, motus.values.T), index=data.index, columns=motus.columns)
    else:
        raise ValueError("?")


amp_name, motus_name = iteratetask(filter_columns(), 2)
#run_corrs(amp_name, motus_name, 'spearmanr')
p = run_corrs(amp_name, motus_name, 'pearsonr')
