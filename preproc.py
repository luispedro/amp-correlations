
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



def filter_human_gut(amp_name, motus_name, min_number_samples=30):
    import pandas as pd
    meta = pd.read_table('data/metadata.tsv', index_col=0)

    motus = pd.read_table(motus_name, index_col=0)
    meta = meta.reindex(motus.index).query('host_tax_id == 9606')
    meta = meta.loc[meta.microontology.map(lambda c: 'host-associated:animal host:digestive tract:in' in c)]
    motus = motus.loc[meta.index]

    motus = motus.T.loc[((motus > 0).sum() >= min_number_samples)].T

    amps = pd.read_table(amp_name, index_col=0)
    amps = amps[motus.index]
    amps = amps.loc[(amps > 0).sum(1) >= min_number_samples]

    omotus = 'preproc/filter_human_gut_motus.tsv.gz'
    oamps = 'preproc/filter_human_gut_amsp.tsv.gz'

    amps.to_csv(oamps, sep='\t')
    motus.to_csv(omotus, sep='\t')
    return oamps, omotus

