
def filter_columns():
    import pandas as pd
    ampsphere = pd.read_table('./data/AMPSphere_v.2021-03.species.tsv.gz', index_col=0)
    matched = {}
    unmatched = {}
    for k,vs in ampsphere.groupby('AMP accession').groups.items():
        vs = ampsphere.loc[vs]
        cs = vs['specI cluster'].value_counts()
        if len(cs) == 1:
            [sp] = cs.index
            matched[k] = sp
        else:
            unmatched[k] = cs
    origin = pd.Series(matched)
    pd.DataFrame({'origin': origin}).to_csv('preproc/AMP_origin.tsv.gz', index_label='AMP', sep='\t')

    interesting = set(origin.index)

    data = []
    for ch in  pd.read_table('data/abundances/amp_abundances_matrix.tsv.gz', index_col=0, chunksize=200):
        ch = ch.loc[ch.index.map(interesting.__contains__)]
        ch.fillna(0., inplace=True)
        print(len(ch))
        data.append(ch)
    data = pd.concat(data)
    data = data.T[data.any()].T

    motus = pd.read_table('data/freeze.v2.motusv2_5.mg3.insertcount.tsv.gz', index_col=0)
    motus.shape

    motus.shape
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


