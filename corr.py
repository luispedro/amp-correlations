import pandas as pd
import numpy as np

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

