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
    elif mode == 'spearmanr':
        return pd.DataFrame(corr.spearman_pcorr(data.values, motus.values.T), index=data.index, columns=motus.columns)
    else:
        raise ValueError("?")
@TaskGenerator
def summarize_correlations(p):
    import pandas as pd

    origins = pd.read_table('preproc/AMP_origin.tsv.gz', index_col=0)
    origin_id = origins.squeeze().map(lambda c : int(c[len('specI_v3_Cluster'):], 10))

    predictions = {}
    for k,v in p.iterrows():
        predictions[k] = {
                'prediction': v.idxmax(),
                'prediction_id': int(v.idxmax().split('_')[-1][:-1], 10),
                'r': v.max(),
                }

    predictions = pd.DataFrame(predictions).T
    predictions['origin_id'] = origin_id.reindex(predictions.index)
    predictions['correct' ] = predictions.eval('prediction_id == origin_id')
    return predictions



amp_name, motus_name = iteratetask(filter_columns(), 2)
p = run_corrs(amp_name, motus_name, 'spearmanr')
s0 = summarize_correlations(p)

p = run_corrs(amp_name, motus_name, 'pearsonr')
s1 = summarize_correlations(p)
