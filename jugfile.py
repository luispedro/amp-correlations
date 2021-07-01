from jug import TaskGenerator, iteratetask
from preproc import filter_columns, filter_human_gut

filter_columns = TaskGenerator(filter_columns)
filter_human_gut = TaskGenerator(filter_human_gut)

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
    taxo = pd.read_table('./data/db_mOTU_taxonomy_ref-mOTUs.tsv', index_col=0)
    taxo.rename(index=lambda ix: int(ix.split('_')[-1], 10), inplace=True)
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
    predictions['prediction_genus'] = taxo['genus'].reindex(predictions['prediction_id']).values
    predictions['origin_genus'] = taxo['genus'].reindex(predictions['origin_id']).values

    predictions['correct' ] = predictions.eval('prediction_id == origin_id')
    predictions['correct_genus' ] = predictions.eval('prediction_genus == origin_genus')

    return predictions

@TaskGenerator
def results_q(s):
    return {
            'nr': len(s),
            'correct': s['correct'].mean(),
            'correct_genus': s['correct_genus'].mean(),
            }


@TaskGenerator
def save_to_tsv(df, oname):
    df.to_csv(oname, sep='\t')
    return oname

amp_name, motus_name = iteratetask(filter_columns(), 2)

final = {}
for min_samples in [20, 30, 45, 60, 100, 120, 150, 200, 250, 500]:
    hg_amp_name, hg_motus_name = iteratetask(filter_human_gut(amp_name, motus_name, min_number_samples=min_samples), 2)

    p = run_corrs(hg_amp_name, hg_motus_name, 'spearmanr')
    s0 = summarize_correlations(p)
    save_to_tsv(s0, f'outputs/spearmanr-hg-results_min={min_samples}.tsv.xz')
    final[min_samples, 'spearmanr'] = results_q(s0)

    p = run_corrs(hg_amp_name, hg_motus_name, 'pearsonr')
    s1 = summarize_correlations(p)
    save_to_tsv(s1, f'outputs/pearsonr-hg-results_min={min_samples}.tsv.xz')
    final[min_samples, 'pearsonr'] = results_q(s1)
