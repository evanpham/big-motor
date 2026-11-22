"""Minimal reader for Elmer's VTU output (appended raw, UInt32 headers)."""
import re, numpy as np
def read(fn):
    d = open(fn, "rb").read()
    tag = b'<AppendedData encoding="raw">'
    i = d.find(tag); head = d[:i].decode(); base = d.find(b"_", i) + 1
    np_ = int(re.search(r'NumberOfPoints="(\d+)"', head).group(1))
    out = {"point": {}, "cell": {}}
    sect = None
    for line in head.splitlines():
        if "<PointData>" in line: sect = "point"
        elif "<CellData>" in line: sect = "cell"
        elif "<Points>" in line: sect = "points"
        elif "<Cells>" in line: sect = "cells"
        m = re.search(r'<DataArray type="(\w+)"(?: Name="([^"]*)")?(?: NumberOfComponents="(\d+)")? format="appended" offset="(\d+)"', line)
        if not m: continue
        typ, name, nc, off = m.group(1), m.group(2) or "points", int(m.group(3) or 1), int(m.group(4))
        dt = {"Float64": "<f8", "Int32": "<i4", "Float32": "<f4"}[typ]
        n = np.frombuffer(d, "<u4", 1, base + off)[0]
        a = np.frombuffer(d, dt, n // np.dtype(dt).itemsize, base + off + 4)
        if nc > 1: a = a.reshape(-1, nc)
        key = sect if sect in ("point", "cell") else sect
        if key in ("point", "cell"): out[key][name] = a
        else: out[name if sect == "cells" else "points"] = a
    off, con, typ = out["offsets"], out["connectivity"], out["types"]
    out["cells"] = [con[(off[k-1] if k else 0):off[k]] for k in range(len(off))]
    out["types"] = typ
    return out

def read_results(path):
    """Read a results folder (serial .vtu or parallel .pvtu + pieces) and merge into one mesh."""
    import glob, os
    if os.path.isfile(path): files = [path]
    else:
        pv = sorted(glob.glob(os.path.join(path, "case_t*.pvtu")))
        if pv:
            txt = open(pv[-1]).read()
            files = [os.path.join(path, f) for f in re.findall(r'<Piece Source="([^"]+)"', txt)]
        else:
            files = sorted(glob.glob(os.path.join(path, "case_t*.vtu")))[-1:]
    parts = [read(f) for f in files]
    if len(parts) == 1: return parts[0]
    out = {"point": {}, "cell": {}}
    shift = np.cumsum([0] + [len(p["points"]) for p in parts[:-1]])
    out["points"] = np.vstack([p["points"].reshape(-1, 3) for p in parts])
    out["connectivity"] = np.concatenate([p["connectivity"] + s for p, s in zip(parts, shift)])
    cshift = np.cumsum([0] + [len(p["connectivity"]) for p in parts[:-1]])
    out["offsets"] = np.concatenate([p["offsets"] + s for p, s in zip(parts, cshift)])
    out["types"] = np.concatenate([p["types"] for p in parts])
    for where in ("point", "cell"):
        for k in parts[0][where]:
            out[where][k] = np.concatenate([p[where][k] for p in parts])
    off, con = out["offsets"], out["connectivity"]
    out["cells"] = [con[(off[k-1] if k else 0):off[k]] for k in range(len(off))]
    return out
