"""Point-probe helpers on Elmer VTU tet meshes (barycentric interpolation)."""
import numpy as np
from scipy.spatial import cKDTree
class Probe:
    def __init__(self, X, tets):
        self.X, self.t = X, tets
        P = X[tets]; self.T0 = P[:, 0]
        M = np.stack([P[:, 1]-P[:, 0], P[:, 2]-P[:, 0], P[:, 3]-P[:, 0]], axis=2)
        self.Minv = np.linalg.inv(M); self.tree = cKDTree(P.mean(1))
    def locate(self, pts, k=24):
        _, cand = self.tree.query(pts, k=k)
        best = -np.ones(len(pts), int); lam = np.zeros((len(pts), 4))
        for j in range(k):
            c = cand[:, j]; l = np.einsum('nij,nj->ni', self.Minv[c], pts - self.T0[c])
            L = np.c_[1 - l.sum(1), l]; ok = (L.min(1) > -1e-9) & (best < 0)
            best[ok] = c[ok]; lam[ok] = L[ok]
        return best, lam
    def interp(self, field, pts):
        e, L = self.locate(pts); ok = e >= 0
        out = np.full((len(pts),) + field.shape[1:], np.nan)
        out[ok] = np.einsum('nk,nk...->n...', L[ok], field[self.t[e[ok]]])
        return out, ok
