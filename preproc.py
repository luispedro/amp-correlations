
def filter_columns():
    from os import makedirs
    import pandas as pd
    makedirs('preproc', exist_ok=True)

    def ex(d):
        if ',' in d: return None
        if d.startswith('specI_v3_'): return int(d[len('specI_v3_'):], 10)
        return int(d[len('ref_mOTU_v25_'):], 10)

    origin = pd.read_table('./data/amp_based_tax_mo.tsv', index_col=0)
    def valid_specI(x):
        return len([t for t in x.split(', ') if t != 'unknown']) == 1
    def fix_specI(x):
        ts = [t for t in x.split(', ') if t != 'unknown']
        if len(ts) > 1:
            return None
        return ts[0]
    origin = origin.loc[origin.specI.map(valid_specI)]
    origin['specI'] = origin.specI.map(fix_specI)
    pd.DataFrame({'origin': origin['specI']}).to_csv('preproc/AMP_origin.tsv.gz', index_label='AMP', sep='\t')

    interesting = set(origin.index)

    data = []
    for ch in  pd.read_table('data/abundances/amp_abundances_matrix.tsv.gz', index_col=0, chunksize=200):
        ch = ch.loc[ch.index.map(interesting.__contains__)]
        ch.fillna(0., inplace=True)
        data.append(ch)
    data = pd.concat(data)
    data = data.T[data.any()].T

    motus = pd.read_table('data/freeze.v2.motusv2_5.mg3.insertcount.tsv.gz', index_col=0)
    motus = motus[motus.any(1)]
    motus = (motus.T/motus.sum(1)).T

    samples = set(data.columns) &  set(motus.index)

    data = data[samples]

    motus = motus.loc[samples]
    mused = motus.any()
    motus = motus.T[mused].T

    amp_name = 'preproc/AMP-abundance.tsv.gz'
    motus_name = 'preproc/mOTUs.tsv.gz'

    data.to_csv(amp_name, sep='\t')
    motus.to_csv(motus_name, sep='\t')
    return amp_name, motus_name



def filter_number_samples(amp_name, motus_name, min_number_samples=30):
    import pandas as pd
    motus = pd.read_table(motus_name, index_col=0, na_filter=False)

    motus = motus.T.loc[((motus > 0).sum() >= min_number_samples)].T

    amps = pd.read_table(amp_name, index_col=0, na_filter=False)
    amps = amps[motus.index]
    amps = amps.loc[(amps > 0).sum(1) >= min_number_samples]

    omotus = f'preproc/filter_motus_min={min_number_samples}.tsv.gz'
    oamps = f'preproc/filter_amps_min={min_number_samples}.tsv.gz'

    amps.to_csv(oamps, sep='\t')
    motus.to_csv(omotus, sep='\t')
    return oamps, omotus

def filter_human_gut(amp_name, motus_name, min_number_samples=30):
    import pandas as pd
    meta = pd.read_table('data/metadata.tsv', index_col=0)

    motus = pd.read_table(motus_name, index_col=0, na_filter=False)
    meta = meta.reindex(motus.index).query('host_tax_id == 9606')
    meta = meta.loc[meta.microontology.map(lambda c: 'host-associated:animal host:digestive tract:in' in c)]
    motus = motus.loc[meta.index]

    motus = motus.T.loc[((motus > 0).sum() >= min_number_samples)].T

    amps = pd.read_table(amp_name, index_col=0, na_filter=False)
    amps = amps[motus.index]
    amps = amps.loc[(amps > 0).sum(1) >= min_number_samples]

    omotus = f'preproc/filter_human_gut_motus_min={min_number_samples}.tsv.gz'
    oamps = f'preproc/filter_human_gut_amps_min={min_number_samples}.tsv.gz'

    amps.to_csv(oamps, sep='\t')
    motus.to_csv(omotus, sep='\t')
    return oamps, omotus

