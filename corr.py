import pandas as pd
import numpy as np
from scipy import stats

def spearman_nz(x, y):
    v = (x > 0) | (y > 0)
    x = x[v]
    y = y[v]
    # We can also do `stats.spearmanr(x,y)[0]`, but that is >2x slower as it
    # computes p-value &c which we discard
    return np.corrcoef(stats.rankdata(x),stats.rankdata(y))[0,1]

def spearman_nzz(x, y):
    v = (x > 0) & (y > 0)
    x = x[v]
    if len(x) <= 8:
        return -2.
    y = y[v]
    # We can also do `stats.spearmanr(x,y)[0]`, but that is >2x slower as it
    # computes p-value &c which we discard
    return np.corrcoef(stats.rankdata(x),stats.rankdata(y))[0,1]

def one_step_spearman_nz_pcorr(X, Y):
    #print("spearman_nz_pcorr")
    #print(X.shape)
    #print(Y.shape)
    X = X.copy()
    Y = Y.copy()
    X[X<=0] = np.nan
    Y[Y<=0] = np.nan
    full = stats.spearmanr(X, Y, axis=1, nan_policy='omit')[0]
    return full[:len(X), len(X):]

def spearman_nzz_pcorr(X, Y):
    print("spearman_nzz_pcorr")
    print(X.shape)
    print(Y.shape)
    res = []
    for i in range(len(X)):
        if i % 10 == 0:
            print(f'iter: {i}')
        res.append([spearman_nzz(X[i], Y[j]) for j in range(len(Y))])
    return np.array(res)
def spearman_nz_pcorr(X, Y):
    print("spearman_nz_pcorr")
    print(X.shape)
    print(Y.shape)
    res = []
    for i in range(len(X)):
        if i % 10 == 0:
            print(f'iter: {i}')
        res.append([spearman_nz(X[i], Y[j]) for j in range(len(Y))])
    return np.array(res)
    #return np.array([
    #        [spearman_nz(X[i], Y[j]) for j in range(len(Y))]
    #            for i in range(len(X))])



def _norm(X):
    XC = X.copy()
    XC -= np.atleast_2d(XC.mean(axis=1)).T
    C2 = np.sqrt((XC**2).sum(axis=1))
    XN = (XC.T / C2).T
    return XN

def pearsonr_pcorr(X, Y):
    return np.dot(_norm(X), _norm(Y).T)
def spearman_pcorr(X, Y):
    from scipy import stats
    return pearsonr_pcorr(stats.rankdata(X, axis=1), stats.rankdata(Y, axis=1))


def test_corr():
    def _slow_pearsonr(X, Y):
        from scipy import stats
        return np.array([
            [stats.pearsonr(X[i], Y[j])[0] for j in range(len(Y))]
                for i in range(len(X))])
    def _slow_spearmanr(X, Y):
        from scipy import stats
        return np.array([
            [stats.spearmanr(X[i], Y[j])[0] for j in range(len(Y))]
                for i in range(len(X))])
    N0 = 120
    N1 = 80
    M = 70

    X = np.random.rand(N0, M)
    Y = np.random.rand(N1, M)

    assert _norm(X).shape == X.shape

    P = _slow_pearsonr(X, Y)
    P2 = pearsonr_pcorr(X, Y)
    assert np.allclose(P, P2)

    S = _slow_spearmanr(X, Y)
    S2 = spearman_pcorr(X, Y)
    assert np.allclose(S, S2)

