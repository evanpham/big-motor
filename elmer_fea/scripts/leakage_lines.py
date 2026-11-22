"""Trace B field lines that start in the middle core's back iron and show only those that do NOT pass
through the middle rotor tooth (leakage), with the iron drawn as outlines.
Seeds are placed in the back-iron cross-section (z = 0) with density proportional to |Bz|, so every
line carries the same flux and the share of drawn lines equals the leakage share of the back-iron flux.
usage: python3 -I leakage_lines.py <run dir, e.g. runs/phase_X_L1/00.0deg> <out.png> [n_seeds=600]"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readvtu, probe
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

run, out = sys.argv[1], sys.argv[2]; nseed = int(sys.argv[3]) if len(sys.argv) > 3 else 600
v = readvtu.read_results(f"{run}/results")
X = v["points"]; T = v["types"]
tet = np.array([c for c, t in zip(v["cells"], T) if t == 10]); g = v["cell"]["GeometryIds"][T == 10]
B = v["point"]["magnetic flux density e"]
P = probe.Probe(X, tet)
rng = np.random.default_rng(1)

# ---- seeds: back-iron cross-section of the middle core at z = 0, density ~ |Bz|
pc = probe.Probe(X, tet[g == 1])
u = np.linspace(70e-3, 100e-3, 300); w = np.linspace(-12e-3, 12e-3, 240); U, W = np.meshgrid(u, w)
pts = np.c_[U.ravel(), W.ravel(), np.zeros(U.size)]
f, ok = pc.interp(B, pts); bz = np.where(ok, np.abs(np.nan_to_num(f[:, 2])), 0)
pts, bz = pts[bz > 0], bz[bz > 0]
pick = rng.choice(len(pts), nseed, replace=False, p=bz / bz.sum())
dA = (u[1]-u[0]) * (w[1]-w[0])
seeds = pts[pick] + np.c_[rng.uniform(-.5, .5, (nseed, 2)) * [u[1]-u[0], w[1]-w[0]], np.zeros(nseed)]
print(f"back-iron flux sampled: {bz.sum()*dA*1e6:.1f} uWb, {nseed} lines -> {bz.sum()*dA/nseed*1e9:.1f} nWb per line")

def field(p):
    e, L = P.locate(p); okk = e >= 0
    b = np.zeros_like(p); b[okk] = np.einsum('nk,nkj->nj', L[okk], B[tet[e[okk]]])
    body = np.zeros(len(p), int); body[okk] = g[e[okk]]
    return b, body, okk

def trace(sign, h=0.15e-3, nmax=1600):
    """RK2 with unit-direction steps; returns paths [n, steps, 3] (nan after exit) and bodies visited"""
    p = seeds.copy(); alive = np.ones(nseed, bool); global closed
    path = np.full((nseed, nmax + 1, 3), np.nan); path[:, 0] = p
    seen = np.zeros((nseed, 7), bool); seen[:, 1] = True
    bod = np.zeros((nseed, nmax + 1), int); bod[:, 0] = 1
    for k in range(nmax):
        a = np.where(alive)[0]
        if len(a) == 0: break
        b1, _, ok1 = field(p[a]); d1 = sign * b1 / (np.linalg.norm(b1, axis=1, keepdims=True) + 1e-30)
        b2, _, ok2 = field(p[a] + 0.5 * h * d1); d2 = sign * b2 / (np.linalg.norm(b2, axis=1, keepdims=True) + 1e-30)
        q = p[a] + h * d2
        _, body, ok3 = field(q)
        good = ok1 & ok2 & ok3
        p[a[good]] = q[good]; path[a[good], k + 1] = q[good]
        seen[a[good], body[good]] = True; bod[a[good], k + 1] = body[good]
        # stop: left the mesh, or closed the loop (back at the start after a long way round)
        back = good & (k > 150) & (np.linalg.norm(q - seeds[a], axis=1) < 0.6e-3)
        alive[a[~good | back]] = False
        closed[a[back]] = True
    return path, seen, bod

closed = np.zeros(nseed, bool)
pf, sf, bf = trace(+1); pb, sb, bb = trace(-1)
seen = sf | sb
# useful = the line runs along the tooth through its mid-plane (same test as the 84 % section flux)
def mid(path, bod):
    z = path[:, :, 2]; cross = (np.sign(z[:, 1:]) != np.sign(z[:, :-1])) & ~np.isnan(z[:, 1:])
    return (cross & (bod[:, 1:] == 3)).any(1)
useful = mid(pf, bf) | mid(pb, bb)
partial = ~useful & seen[:, 3]
to_neigh = ~useful & ~partial & (seen[:, 5] | seen[:, 6])
local = ~useful & ~partial & ~to_neigh
lost = 1 - useful.mean()
print(f"along the middle rotor tooth (useful): {useful.mean()*100:.1f}%   leakage: {lost*100:.1f}%")
print(f"  touches the tooth then short-cuts through air {partial.mean()*100:.1f}%, to neighbour poles {to_neigh.mean()*100:.1f}%, "
      f"never reaches the tooth {local.mean()*100:.1f}%   (lines that ended without closing: {(~closed).sum()})")

# ---- iron outlines: boundary faces of each iron body (a face belonging to one tet of that body)
_, nid = np.unique(np.round(X / 1e-7).astype(np.int64), axis=0, return_inverse=True); nid = nid.ravel()
def skin(ids):
    tt = tet[np.isin(g, ids)]
    tn = nid[tt]
    F = np.concatenate([tt[:, [0, 1, 2]], tt[:, [0, 1, 3]], tt[:, [0, 2, 3]], tt[:, [1, 2, 3]]])
    Fn = np.concatenate([tn[:, [0, 1, 2]], tn[:, [0, 1, 3]], tn[:, [0, 2, 3]], tn[:, [1, 2, 3]]])
    Fs = np.sort(Fn, 1); _, inv, cnt = np.unique(Fs, axis=0, return_inverse=True, return_counts=True)
    return X[F[cnt[inv.ravel()] == 1]]
bodies = {"middle core": (1, "#7f7f7f", 0.10), "middle rotor tooth": (3, "#1f77b4", 0.12),
          "other cores": (5, "#7f7f7f", 0.04), "other teeth": (6, "#1f77b4", 0.04)}
skins = {k: skin([i]) for k, (i, _, _) in bodies.items()}

views = [("3D view", 22, -58), ("looking along the motor axis (z)", 90, -90), ("looking tangentially (side of the C)", 0, -90), ("looking radially inward", 0, 0)]
fig = plt.figure(figsize=(16, 13), constrained_layout=True)
for i, (title, el, az) in enumerate(views):
    ax = fig.add_subplot(2, 2, i + 1, projection="3d")
    for k, (bid, col, al) in bodies.items():
        ax.add_collection3d(Poly3DCollection(skins[k] * 1e3, facecolor=col, edgecolor="none", alpha=al))
    for mask, col, lab in ((local, "#d62728", f"never reaches the tooth ({local.mean()*100:.1f}% of back-iron flux)"),
                           (partial, "#9467bd", f"touches the tooth, short-cuts through air ({partial.mean()*100:.1f}%)"),
                           (to_neigh, "#ff7f0e", f"to neighbour poles ({to_neigh.mean()*100:.1f}%)")):
        for j, n in enumerate(np.where(mask)[0]):
            line = np.vstack([pb[n][::-1], pf[n][1:]]) * 1e3
            line = line[~np.isnan(line[:, 0])]
            ax.plot(line[:, 0], line[:, 1], line[:, 2], color=col, lw=0.7, alpha=0.8, label=lab if j == 0 else None)
    ax.set_xlim(70, 125); ax.set_ylim(-27.5, 27.5); ax.set_zlim(-45, 45); ax.set_box_aspect((55, 55, 90))
    ax.view_init(el, az); ax.set_proj_type("ortho"); ax.set_title(title)
    ax.set_xlabel("radial x [mm]"); ax.set_ylabel("tangential y [mm]"); ax.set_zlabel("axial z [mm]")
    if i == 0: ax.legend(loc="upper left", fontsize=9)
fig.suptitle(f"{os.path.basename(os.path.dirname(os.path.abspath(run)))}  {os.path.basename(run)}: back-iron flux that does NOT pass through "
             f"the middle rotor tooth ({lost*100:.1f}% of it; each line = equal flux; {useful.sum()} useful lines hidden)", fontsize=12)
fig.savefig(out, dpi=110)
json.dump(dict(useful=float(useful.mean()), partial=float(partial.mean()), to_neighbours=float(to_neigh.mean()), local=float(local.mean()), n=nseed),
          open(os.path.splitext(out)[0] + ".json", "w"))
np.savez_compressed(os.path.splitext(out)[0] + "_lines.npz", pf=pf, pb=pb, useful=useful, partial=partial, to_neigh=to_neigh, local=local)
print("wrote", out)
