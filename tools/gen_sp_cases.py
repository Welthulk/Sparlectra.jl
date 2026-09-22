#!/usr/bin/env python3
"""Synthetic MATPOWER cases shipped with Sparlectra.jl.

    python3 tools/gen_sp_cases.py            # writes all three into data/mpower
    python3 tools/gen_sp_cases.py 300 out.m  # one case, explicit path

Three cases, all self-built (no external data, no external license):

- sp_case9    a 9-bus system in the shape of the classic 3-machine test
              system (three generator buses behind step-up transformers
              feeding a 230 kV ring with three loads), hand-built; the
              file carries the DC power-flow reference values this script
              computes with its own DC solver, which the DC test uses as
              an independent anchor
- sp_case300  a 300-bus operated grid: three voltage levels, three
              areas, an EHV backbone, generated
- sp_case1354 a 1354-bus grid of the same construction, the size probe
              for the sparse code paths (FD coloring, N-1 batches, timings)

The generated cases are validated in this script with a Newton-Raphson
power flow (vectorized Jacobian) including generator Q-limit enforcement,
an island check, and an N-1 sweep (full for the small cases, a sample for
the large one).
"""
import math
import os
import random
import sys

import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.sparse.linalg import spsolve

BASE = 100.0

# ------------------------------------------------------------ AC power flow --
class Case:
  def __init__(self, name, buses, gens, branches, bus_name, slack, comment):
    self.name = name
    self.bus = buses            # rows of the MATPOWER bus matrix
    self.gen = gens             # rows of the MATPOWER gen matrix
    self.branch = branches      # rows of the MATPOWER branch matrix
    self.bus_name = bus_name    # list of names, one per bus row
    self.slack = slack          # bus number
    self.comment = comment      # list of header comment lines (filled later)
    self.N = len(buses)

def ybus(case, status_override=None):
  N = case.N
  Y = lil_matrix((N, N), dtype=complex)
  for k, row in enumerate(case.branch):
    f, t, r, x, b = int(row[0]), int(row[1]), row[2], row[3], row[4]
    ratio, angle = row[8], row[9]
    status = row[10] if status_override is None else status_override[k]
    if not status:
      continue
    ys = 1.0 / complex(r, x)
    tap = (ratio if ratio else 1.0) * np.exp(1j * math.radians(angle))
    i, j = f - 1, t - 1
    Y[i, i] += (ys + 1j * b / 2) / (abs(tap) ** 2)
    Y[j, j] += ys + 1j * b / 2
    Y[i, j] += -ys / np.conj(tap)
    Y[j, i] += -ys / tap
  for row in case.bus:
    i = int(row[0]) - 1
    Y[i, i] += complex(row[4], row[5]) / BASE
  return Y.tocsr()

def is_connected(case, status_override=None):
  N = case.N
  A = lil_matrix((N, N))
  for k, row in enumerate(case.branch):
    status = row[10] if status_override is None else status_override[k]
    if status:
      A[int(row[0]) - 1, int(row[1]) - 1] = 1
      A[int(row[1]) - 1, int(row[0]) - 1] = 1
  n, _ = connected_components(A.tocsr(), directed=False)
  return n == 1

def newton(case, enforce_q=True, max_iter=30, tol=1e-8, warm=None, status_override=None):
  N = case.N
  Y = ybus(case, status_override)
  Pd = np.array([r[2] for r in case.bus]) / BASE
  Qd = np.array([r[3] for r in case.bus]) / BASE
  Pg = np.zeros(N); Qg = np.zeros(N)
  vset = {}; qmax = {}; qmin = {}
  for g in case.gen:
    if not g[7]:
      continue
    i = int(g[0]) - 1
    Pg[i] += g[1] / BASE
    vset[i] = g[5]; qmax[i] = qmax.get(i, 0.0) + g[3] / BASE; qmin[i] = qmin.get(i, 0.0) + g[4] / BASE
  slack = case.slack - 1
  pv = set(vset) - {slack}
  pq = set(range(N)) - pv - {slack}
  V = np.ones(N); th = np.zeros(N)
  for i, v in vset.items():
    V[i] = v
  if warm is not None:
    V = warm[0].copy(); th = warm[1].copy()
    for i, v in vset.items():
      V[i] = v
  clamped = {}
  Pspec = Pg - Pd
  it = 0
  for outer in range(15):
    Qspec = Qg - Qd
    for it in range(max_iter):
      Vc = V * np.exp(1j * th)
      Ibus = Y @ Vc
      S = Vc * np.conj(Ibus)
      pvpq = np.array(sorted(pv | pq)); pql = np.array(sorted(pq))
      dP = Pspec[pvpq] - S.real[pvpq]
      dQ = Qspec[pql] - S.imag[pql] if len(pql) else np.zeros(0)
      mis = np.concatenate([dP, dQ])
      if np.max(np.abs(mis)) < tol:
        break
      # dS/dVa, dS/dVm (dense; N <= 1400)
      Yd = Y.toarray()
      diagV = np.diag(Vc); diagI = np.diag(Ibus); diagVn = np.diag(Vc / np.abs(Vc))
      dS_dVm = diagV @ np.conj(Yd @ diagVn) + np.conj(diagI) @ diagVn
      dS_dVa = 1j * diagV @ np.conj(diagI - Yd @ diagV)
      J11 = dS_dVa[np.ix_(pvpq, pvpq)].real
      J12 = dS_dVm[np.ix_(pvpq, pql)].real
      J21 = dS_dVa[np.ix_(pql, pvpq)].imag
      J22 = dS_dVm[np.ix_(pql, pql)].imag
      J = np.block([[J11, J12], [J21, J22]])
      dx = np.linalg.solve(J, mis)
      th[pvpq] += dx[:len(pvpq)]
      V[pql] += dx[len(pvpq):]
      if np.max(np.abs(dx)) > 50:
        return None
    else:
      return None
    if not enforce_q:
      break
    Vc = V * np.exp(1j * th)
    S = Vc * np.conj(Y @ Vc)
    Qgen = S.imag + Qd
    switched = False
    for i in sorted(pv):
      if Qgen[i] > qmax[i] + 1e-9:
        pv.discard(i); pq.add(i); Qg[i] = qmax[i]; clamped[i] = "max"; switched = True
      elif Qgen[i] < qmin[i] - 1e-9:
        pv.discard(i); pq.add(i); Qg[i] = qmin[i]; clamped[i] = "min"; switched = True
    if not switched:
      break
  Vc = V * np.exp(1j * th)
  S = Vc * np.conj(Y @ Vc)
  loads = np.zeros(len(case.branch))
  for k, row in enumerate(case.branch):
    status = row[10] if status_override is None else status_override[k]
    if not status:
      continue
    f, t, r, x, b, rateA = int(row[0]), int(row[1]), row[2], row[3], row[4], row[5]
    ratio, angle = row[8], row[9]
    ys = 1.0 / complex(r, x)
    tap = (ratio if ratio else 1.0) * np.exp(1j * math.radians(angle))
    i, j = f - 1, t - 1
    If = (Vc[i] / tap - Vc[j]) * ys / np.conj(tap) + Vc[i] / (abs(tap) ** 2) * 1j * b / 2
    It = (Vc[j] - Vc[i] / tap) * ys + Vc[j] * 1j * b / 2
    loads[k] = 100.0 * max(abs(Vc[i] * np.conj(If)), abs(Vc[j] * np.conj(It))) * BASE / rateA
  return dict(V=V, th=th, S=S, iters=it + 1, clamped=clamped, loads=loads,
              pslack=(S.real[slack] + Pd[slack]) * BASE, Qgen=(S.imag + Qd) * BASE)

def dc_powerflow(case):
  """DC power flow the MATPOWER way (B' from x and tap, shift as Pfinj)."""
  N = case.N
  B = lil_matrix((N, N)); Pfinj = np.zeros(N)
  for row in case.branch:
    if not row[10]:
      continue
    f, t, x, ratio, angle = int(row[0]) - 1, int(row[1]) - 1, row[3], row[8], row[9]
    tap = ratio if ratio else 1.0
    bft = 1.0 / (x * tap)
    B[f, f] += bft; B[t, t] += bft; B[f, t] -= bft; B[t, f] -= bft
    pshift = -bft * math.radians(angle)
    Pfinj[f] += pshift; Pfinj[t] -= pshift
  P = np.zeros(N)
  for g in case.gen:
    if g[7]:
      P[int(g[0]) - 1] += g[1] / BASE
  for row in case.bus:
    P[int(row[0]) - 1] -= row[2] / BASE + row[4] / BASE
  P -= Pfinj
  slack = case.slack - 1
  keep = [i for i in range(N) if i != slack]
  Bk = B.tocsr()[keep][:, keep]
  th = np.zeros(N)
  th[keep] = spsolve(Bk.tocsc(), P[keep])
  flows = []
  for row in case.branch:
    f, t, x, ratio, angle = int(row[0]) - 1, int(row[1]) - 1, row[3], row[8], row[9]
    tap = ratio if ratio else 1.0
    flows.append((th[f] - th[t] - math.radians(angle)) / (x * tap) * BASE)
  pslack = float((B.tocsr()[slack] @ th)[0] + Pfinj[slack]) * BASE + case.bus[slack][2] + case.bus[slack][4]
  return th, np.array(flows), pslack

# ------------------------------------------------------------- generator ----
PREFIXES = ["Alt", "Neu", "Hoch", "Nieder", "Ober", "Unter", "Wald", "Berg", "Tal", "Feld", "Bach",
            "See", "Stein", "Hain", "Lind", "Eich", "Birk", "Buch", "Wiesen", "Moor", "Sand", "Kies",
            "Rot", "Weiss", "Gruen", "Blau", "Gold", "Silber", "Nord", "Sued", "Ost", "West", "Kalt",
            "Warm", "Klein", "Gross", "Frei", "Fern", "Nah", "Lang", "Kurz", "Breit", "Tief", "Flach",
            "Ahorn", "Erlen", "Espen", "Ulmen", "Tannen", "Fichten", "Kiefern", "Hasel", "Holler",
            "Dorn", "Ginster", "Heide", "Wacholder", "Schilf", "Binsen", "Klee", "Flachs", "Hanf"]
SUFFIXES = ["dorf", "heim", "hausen", "feld", "berg", "bach", "tal", "hof", "brueck", "furt",
            "kirchen", "stadt", "au", "roda", "walde", "hagen", "wang", "loh", "ried", "eck",
            "born", "moos", "stein", "burg", "weiler", "leben", "ingen", "beuren", "ach", "see"]

def make_names(N, rng):
  names = []; seen = set()
  while len(names) < N:
    n = rng.choice(PREFIXES) + rng.choice(SUFFIXES)
    if len(names) >= len(PREFIXES) * len(SUFFIXES) * 0.9:
      n += str(len(names))
    if n in seen:
      continue
    seen.add(n); names.append(n)
  return names

def build_synthetic(name, N, n_areas, n_ehv, n_trafo, n_branches, n_gen, n_load, total_load_mw,
                    seed, hv_kv=110.0, ehv_kv=380.0, ehv_share_mw=0.55, n_radial=None, n_clamp=None):
  rng = random.Random(seed)
  buses = list(range(1, N + 1))
  # areas: contiguous index ranges
  bounds = [round(i * N / n_areas) for i in range(n_areas + 1)]
  area_of = {}
  areas = {}
  for a in range(n_areas):
    areas[a + 1] = list(range(bounds[a] + 1, bounds[a + 1] + 1))
    for b in areas[a + 1]:
      area_of[b] = a + 1
  # EHV buses spread over the index range
  ehv = sorted(set(int(round((k + 0.5) * N / n_ehv)) for k in range(n_ehv)))
  while len(ehv) < n_ehv:
    c = rng.randint(1, N)
    c in ehv or ehv.append(c)
  ehv = sorted(ehv); ehv_set = set(ehv)
  hv = [b for b in buses if b not in ehv_set]
  # layout: areas on a horizontal band, EHV on a corridor above
  side = math.sqrt(N) * 6.0
  centres = {a: ((a - 1) * side * 1.2, 0.0) for a in areas}
  pos = {}
  for a, bl in areas.items():
    cx, cy = centres[a]
    for b in bl:
      if b not in ehv_set:
        pos[b] = (cx + rng.uniform(-side / 2, side / 2), cy + rng.uniform(-side / 2.4, side / 2.4))
  xs = np.linspace(-side / 2, (n_areas - 1) * side * 1.2 + side / 2, len(ehv))
  for x, b in zip(xs, ehv):
    pos[b] = (float(x), side * 0.7 + rng.uniform(-side * 0.1, side * 0.1))

  def dist(u, v):
    return math.hypot(pos[u][0] - pos[v][0], pos[u][1] - pos[v][1])

  edges = {}   # (min,max) -> kind
  def add(u, v, kind):
    key = (min(u, v), max(u, v))
    if u == v or key in edges:
      return False
    edges[key] = kind
    return True

  # HV grid per area: Delaunay candidates, MST first, then shortest chords
  delaunay_cands = {}
  for a, bl in areas.items():
    hvb = [b for b in bl if b not in ehv_set]
    pts = np.array([pos[b] for b in hvb])
    tri = Delaunay(pts)
    cand = set()
    for s in tri.simplices:
      for i in range(3):
        u, v = hvb[s[i]], hvb[s[(i + 1) % 3]]
        cand.add((min(u, v), max(u, v)))
    delaunay_cands[a] = sorted(cand, key=lambda e: dist(*e))
    # MST over the candidates
    idx = {b: i for i, b in enumerate(hvb)}
    W = lil_matrix((len(hvb), len(hvb)))
    for u, v in cand:
      W[idx[u], idx[v]] = dist(u, v)
    T = minimum_spanning_tree(W.tocsr()).tocoo()
    for i, j in zip(T.row, T.col):
      add(hvb[i], hvb[j], "L")
  # EHV backbone chain plus chords
  for i in range(len(ehv) - 1):
    add(ehv[i], ehv[i + 1], "E")
  for i in range(0, len(ehv) - 2, 3):
    add(ehv[i], ehv[i + 2], "E")
  # transformers: chosen EHV buses to the nearest HV bus of their area
  trafo_ehv = [ehv[int(round(k * (len(ehv) - 1) / max(n_trafo - 1, 1)))] for k in range(n_trafo)]
  trafo_ehv = sorted(set(trafo_ehv))
  k = 0
  while len(trafo_ehv) < n_trafo:
    e = ehv[k % len(ehv)]; k += 1
    e in trafo_ehv or trafo_ehv.append(e)
  trafo_hv = set()   # one transformer per HV bus, so no bus collects several injections
  for e in trafo_ehv:
    a = area_of[e]
    for _, b in sorted((dist(e, b), b) for b in areas[a] if b not in ehv_set and b not in trafo_hv):
      if add(e, b, "T"):
        trafo_hv.add(b)
        break
  # inter-area HV ties: nearest pairs across neighbouring areas
  for a in range(1, n_areas):
    left = [b for b in areas[a] if b not in ehv_set]
    right = [b for b in areas[a + 1] if b not in ehv_set]
    pairs = sorted(((dist(u, v), u, v) for u in left for v in right), key=lambda p: p[0])
    n_ties = max(2, N // 60)
    added = 0
    for _, u, v in pairs:
      if added >= n_ties:
        break
      if add(u, v, "L"):
        added += 1
  # top up with shortest Delaunay chords to the target count, bridges last
  cands_all = sorted((e for a in delaunay_cands for e in delaunay_cands[a]), key=lambda e: dist(*e))
  ci = 0
  while len(edges) < n_branches and ci < len(cands_all):
    u, v = cands_all[ci]; ci += 1
    add(u, v, "L")
  assert len(edges) >= n_branches, f"{name}: only {len(edges)} branches possible, {n_branches} wanted"
  # transformer HV buses need at least three outlets
  deg = {b: 0 for b in buses}
  for (u, v) in edges:
    deg[u] += 1; deg[v] += 1
  for b in sorted(trafo_hv):
    for _, c in sorted((dist(b, c), c) for c in hv if c != b):
      if deg[b] >= 4:
        break
      if add(b, c, "L"):
        deg[b] += 1; deg[c] += 1
  # radial stubs: leave n_radial HV leaves, give every other bus degree >= 2
  deg = {b: 0 for b in buses}
  for (u, v) in edges:
    deg[u] += 1; deg[v] += 1
  leaves = [b for b in hv if deg[b] == 1]
  n_radial = len(leaves) // 4 if n_radial is None else n_radial
  keep_radial = set(leaves[:n_radial])
  for b in leaves:
    if b in keep_radial:
      continue
    for _, c in sorted((dist(b, c), c) for c in hv if c != b and c not in keep_radial):
      if add(b, c, "L"):
        break
  # trim to the exact count: longest HV chords whose ends keep degree >= 3
  # and whose removal creates no bridge
  def bridges_count(edge_keys):
    idx = {b: i for i, b in enumerate(buses)}
    A = lil_matrix((N, N))
    for (u, v) in edge_keys:
      A[idx[u], idx[v]] = 1; A[idx[v], idx[u]] = 1
    A = A.tocsr()
    # Tarjan on an adjacency list
    adj = [[] for _ in range(N)]
    for (u, v) in edge_keys:
      adj[idx[u]].append(idx[v]); adj[idx[v]].append(idx[u])
    disc = [-1] * N; low = [0] * N; cnt = 0; timer = 0
    sys.setrecursionlimit(20000)
    def dfs(u, p):
      nonlocal cnt, timer
      disc[u] = low[u] = timer; timer += 1
      skip_parent_once = True
      for v in adj[u]:
        if v == p and skip_parent_once:
          skip_parent_once = False
          continue
        if disc[v] != -1:
          low[u] = min(low[u], disc[v])
        else:
          dfs(v, u)
          low[u] = min(low[u], low[v])
          if low[v] > disc[u]:
            cnt += 1
    for s in range(N):
      if disc[s] == -1:
        dfs(s, -1)
    return cnt
  base_bridges = bridges_count(list(edges))
  while len(edges) > n_branches:
    deg = {b: 0 for b in buses}
    for (u, v) in edges:
      deg[u] += 1; deg[v] += 1
    cands = sorted(((dist(u, v), (u, v)) for (u, v), kd in edges.items()
                    if kd == "L" and deg[u] >= 3 and deg[v] >= 3 and u not in trafo_hv and v not in trafo_hv), reverse=True)
    removed = False
    for _, key in cands[:300]:
      trial = [e for e in edges if e != key]
      A = lil_matrix((N, N))
      for (u, v) in trial:
        A[u - 1, v - 1] = 1; A[v - 1, u - 1] = 1
      if connected_components(A.tocsr(), directed=False)[0] == 1 and bridges_count(trial) <= base_bridges:
        del edges[key]; removed = True
        break
    assert removed, f"{name}: cannot trim to {n_branches} without islanding"
  assert len(edges) == n_branches

  # ------------------------------------------------ electrical parameters --
  zb_hv = hv_kv ** 2 / BASE
  zb_ehv = ehv_kv ** 2 / BASE
  scale_km = 40.0 / side          # map layout units to kilometres
  branch_rows = []
  kinds = []
  for (u, v), kind in edges.items():
    if kind == "L":
      km = 4.0 + dist(u, v) * scale_km * 1.2
      r, x, b = 0.12 * km / zb_hv, 0.39 * km / zb_hv, 2.9e-6 * km * zb_hv
      rate = 200.0
    elif kind == "E":
      km = 30.0 + dist(u, v) * scale_km * 2.5
      r, x, b = 0.03 * km / zb_ehv, 0.30 * km / zb_ehv, 4.0e-6 * km * zb_ehv
      rate = 1500.0
    else:
      r, x, b = 0.0, 0.12 * BASE / 600.0, 0.0
      rate = 600.0
    ratio = round(rng.choice([0.96, 0.97, 0.98, 0.985, 1.0, 1.0]), 3) if kind == "T" else 0.0
    f, t = (u, v) if kind != "T" else ((u, v) if u in ehv_set else (v, u))
    branch_rows.append([f, t, round(r, 6), round(x, 6), round(b, 5), rate, rate, rate, ratio, 0.0, 1, -360.0, 360.0])
    kinds.append(kind)

  # ------------------------------------------------------- generation ------
  gen_buses = list(ehv)
  hv_gen = [b for b in hv[:: max(1, len(hv) // max(n_gen - len(ehv), 1))]][: n_gen - len(ehv)]
  gen_buses += hv_gen
  gen_buses = sorted(set(gen_buses))[:n_gen]
  slack = ehv[len(ehv) // 2]
  # loads
  load_buses = sorted(rng.sample([b for b in hv if b not in set(hv_gen[::3])], n_load))
  raw = [min(rng.lognormvariate(3.0, 0.6), 95.0) for _ in load_buses]
  sc = total_load_mw / sum(raw)
  pd = {b: round(w * sc) for b, w in zip(load_buses, raw)}
  qd = {b: round(pd[b] * rng.uniform(0.25, 0.40)) for b in load_buses}
  # plants: EHV plants carry ehv_share, HV plants the rest, balanced per area
  plants = {}
  area_load = {a: sum(pd.get(b, 0) for b in bl) for a, bl in areas.items()}
  for a, bl in areas.items():
    target = area_load[a] * 1.03
    ehv_a = [b for b in gen_buses if b in ehv_set and area_of[b] == a]
    hv_a = [b for b in gen_buses if b not in ehv_set and area_of[b] == a]
    e_share = target * ehv_share_mw if ehv_a else 0.0
    for b in ehv_a:
      plants[b] = round(e_share / len(ehv_a))
    h_share = target - e_share
    weights = [rng.uniform(0.3, 1.0) for _ in hv_a]
    for b, w in zip(hv_a, weights):
      plants[b] = round(h_share * w / sum(weights))
  # small HV units become synchronous condensers
  small = [b for b, p in plants.items() if b not in ehv_set and p < 25]
  moved = sum(plants[b] for b in small)
  for b in small:
    plants[b] = 0
  big = [b for b in plants if plants[b] > 0 and b != slack]
  if big:
    scale_up = 1.0 + moved / sum(plants[b] for b in big)
    for b in big:
      plants[b] = round(plants[b] * scale_up)
  plants[slack] = 0
  gen_rows = []
  for b in gen_buses:
    p = float(plants.get(b, 0))
    if b == slack:
      pmax = max(300.0, round(area_load[area_of[b]] * 0.5 / 10) * 10)
    else:
      pmax = max(50.0, round(p * 1.3 / 10.0) * 10.0)
    vg = rng.choice([1.0, 1.005, 1.01, 1.015, 1.02]) if (b in ehv_set or p > 0) else rng.choice([0.99, 0.995, 1.0, 1.005])
    gen_rows.append([b, p, 0.0, 50.0, -50.0, vg, BASE, 1, pmax, 0.0])
  # shunts: capacitor banks at a few heavy-load HV buses
  shunts = {}
  for b in sorted(load_buses, key=lambda b: -pd[b])[: max(3, N // 12)]:
    shunts[b] = (0.0, float(rng.choice([15, 20, 25, 30])))
  names = make_names(N, rng)
  bus_rows = []
  bus_name = []
  for b in buses:
    gs, bs = shunts.get(b, (0.0, 0.0))
    btype = 3 if b == slack else (2 if b in gen_buses else 1)
    kv = ehv_kv if b in ehv_set else hv_kv
    bus_rows.append([b, btype, float(pd.get(b, 0)), float(qd.get(b, 0)), gs, bs, area_of[b], 1.0, 0.0, kv, area_of[b], 1.06, 0.94])
    bus_name.append(f"{names[b - 1]}_{int(kv)}")
  case = Case(name, bus_rows, gen_rows, branch_rows, bus_name, slack, [])
  case.kinds = kinds
  case.pos = pos
  case.hv = hv
  case.protected = set(trafo_hv)
  case.zb_hv = zb_hv
  case.scale_km = scale_km
  case.rng = rng
  case.n_clamp = max(4, n_gen // 8) if n_clamp is None else n_clamp
  return case

def tune_and_validate(case, n1_full=True, n1_sample=0):
  """Q-limits from the unlimited solve (some deliberately below the need),
  ratings from the base flow, then the checks. Returns the summary dict."""
  rng = case.rng
  assert is_connected(case), f"{case.name}: not connected"
  res0 = newton(case, enforce_q=False)
  assert res0 is not None, f"{case.name}: unlimited base case did not converge"
  Qg0 = res0["Qgen"]
  gen_by_bus = {int(g[0]): g for g in case.gen}
  cand = [b for b in gen_by_bus if b != case.slack]
  rng.shuffle(cand)
  clamp = set(cand[: case.n_clamp])
  def r5(v):
    return 5.0 * round(v / 5.0)
  for b, g in gen_by_bus.items():
    if b == case.slack:
      g[3], g[4] = 9999.0, -9999.0
      continue
    q0 = Qg0[b - 1]
    if b in clamp:
      if q0 >= 0:
        g[3] = max(10.0, r5(0.75 * q0)); g[4] = -max(10.0, r5(0.4 * abs(q0)))
      else:
        g[4] = min(-10.0, r5(0.75 * q0)); g[3] = max(10.0, r5(0.4 * abs(q0)))
    else:
      g[3] = max(20.0, r5(1.3 * max(q0, 0.0) + 20.0), r5(0.3 * g[8]) if g[1] > 0 else 0.0)
      g[4] = min(-20.0, -r5(1.3 * max(-q0, 0.0) + 20.0))
  res = newton(case)
  assert res is not None, f"{case.name}: base case did not converge with limits"
  classes = {"L": [150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0],
             "E": [900.0, 1200.0, 1500.0, 2000.0, 2500.0], "T": [400.0, 600.0, 800.0, 1000.0, 1500.0]}
  for k, row in enumerate(case.branch):
    flow = res["loads"][k] / 100.0 * row[5]
    cls = classes[case.kinds[k]]
    rate = next((c for c in cls if c >= 1.3 * flow), cls[-1])
    row[5] = row[6] = row[7] = rate
  # heavy corridors must be meshed: a bus on a branch rated >= 800 MVA needs
  # at least three outlets, otherwise the outage of one corridor link is a
  # voltage collapse. Add the nearest chords, drop as many light long chords.
  for _pass in range(3):
    deg = {}
    for row in case.branch:
      deg[int(row[0])] = deg.get(int(row[0]), 0) + 1
      deg[int(row[1])] = deg.get(int(row[1]), 0) + 1
    weak = sorted({int(b) for k, row in enumerate(case.branch) if row[5] >= 800.0 and case.kinds[k] == "L"
                   for b in (row[0], row[1]) if deg[int(b)] < 3})
    if not weak:
      break
    existing = {(min(int(r[0]), int(r[1])), max(int(r[0]), int(r[1]))) for r in case.branch}
    added = 0
    for b in weak:
      for _, c in sorted((math.hypot(case.pos[b][0] - case.pos[c][0], case.pos[b][1] - case.pos[c][1]), c) for c in case.hv if c != b):
        key = (min(b, c), max(b, c))
        if deg[b] >= 3:
          break
        if key in existing:
          continue
        km = 4.0 + math.hypot(case.pos[b][0] - case.pos[c][0], case.pos[b][1] - case.pos[c][1]) * case.scale_km * 1.2
        r, x, bb = 0.12 * km / case.zb_hv, 0.39 * km / case.zb_hv, 2.9e-6 * km * case.zb_hv
        case.branch.append([b, c, round(r, 6), round(x, 6), round(bb, 5), 200.0, 200.0, 200.0, 0.0, 0.0, 1, -360.0, 360.0])
        case.kinds.append("L")
        existing.add(key); deg[b] += 1; deg[c] = deg.get(c, 0) + 1; added += 1
    # drop the same number of light, long, non-bridge chords between well-connected buses
    to_remove = []
    st = [r[10] for r in case.branch]
    for k in sorted(range(len(case.branch)), key=lambda k: -case.branch[k][3]):
      if len(to_remove) >= added:
        break
      row = case.branch[k]
      f, t = int(row[0]), int(row[1])
      if case.kinds[k] != "L" or row[5] > 200.0 or deg[f] < 4 or deg[t] < 4 or f in case.protected or t in case.protected:
        continue
      st[k] = 0
      if not is_connected(case, st):
        st[k] = 1
        continue
      deg[f] -= 1; deg[t] -= 1
      to_remove.append(k)
    for k in sorted(to_remove, reverse=True):
      del case.branch[k]; del case.kinds[k]
    res = newton(case)
    assert res is not None, f"{case.name}: base case did not converge after corridor meshing"
    for k, row in enumerate(case.branch):
      flow = res["loads"][k] / 100.0 * row[5]
      cls = classes[case.kinds[k]]
      row[5] = row[6] = row[7] = next((c for c in cls if c >= 1.3 * flow), cls[-1])
  res = newton(case)
  V = res["V"]
  summary = dict(iters=res["iters"], vmin=float(V.min()), vmax=float(V.max()), loading=float(res["loads"].max()),
                 pslack=float(res["pslack"]), clamped=len(res["clamped"]), n_gen=len(case.gen), n_branch=len(case.branch),
                 load_mw=sum(r[2] for r in case.bus), load_mvar=sum(r[3] for r in case.bus))
  assert 0.94 <= V.min() and V.max() <= 1.06, f"{case.name}: voltage band {V.min():.3f}..{V.max():.3f}"
  assert res["loads"].max() < 90.0, f"{case.name}: base overload {res['loads'].max():.0f} %"
  # N-1 over the branches
  ks = list(range(len(case.branch)))
  if not n1_full:
    rng.shuffle(ks); ks = sorted(ks[:n1_sample])
  n_island = n_fail = 0
  for k in ks:
    st = [r[10] for r in case.branch]; st[k] = 0
    if not is_connected(case, st):
      n_island += 1; continue
    if newton(case, warm=(res["V"], res["th"]), status_override=st) is None:
      n_fail += 1
  summary.update(n1_checked=len(ks), n1_island=n_island, n1_fail=n_fail)
  return summary, res

# ---------------------------------------------------------------- writer ---
def fmt(v):
  if isinstance(v, float):
    s = f"{v:.6g}"
    return s if ("." in s or "e" in s or s in ("inf", "nan")) else s + ".0"
  return str(v)

def write_case(case, path, header):
  out = [f"function mpc = {case.name}"]
  out += [("% " + line).rstrip() for line in header]
  out += ["", "%% MATPOWER Case Format : Version 2", "mpc.version = '2';", "",
          "%%-----  Power Flow Data  -----%%", "%% system MVA base", "mpc.baseMVA = 100;", "",
          "%% bus data", "%\tbus_i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin", "mpc.bus = ["]
  out += ["\t" + "\t".join(fmt(v) for v in r) + ";" for r in case.bus]
  out += ["];", "", "%% generator data", "%\tbus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin", "mpc.gen = ["]
  out += ["\t" + "\t".join(fmt(v) for v in r) + ";" for r in case.gen]
  out += ["];", "", "%% branch data",
          "%\tfbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\tangmin\tangmax", "mpc.branch = ["]
  out += ["\t" + "\t".join(fmt(v) for v in r) + ";" for r in case.branch]
  out += ["];", "", "%% bus names (imported with matpower_import.apply_bus_names)", "mpc.bus_name = {"]
  out += [f"\t'{n}';" for n in case.bus_name]
  out += ["};", ""]
  with open(path, "w", newline="\n") as fh:
    fh.write("\n".join(out))

# ------------------------------------------------------------- sp_case9 ----
def build_case9():
  # three machines behind step-up transformers, a 230 kV ring with three loads
  kv = {1: 16.5, 2: 18.0, 3: 13.8, 4: 230.0, 5: 230.0, 6: 230.0, 7: 230.0, 8: 230.0, 9: 230.0}
  names = ["Kraftwerk_Nord", "Kraftwerk_Ost", "Kraftwerk_Sued", "Ringpunkt_Nord", "Stadtwerk_West",
           "Industriepark", "Ringpunkt_Ost", "Umspannwerk_Mitte", "Ringpunkt_Sued"]
  bus = []
  loads = {5: (90.0, 30.0), 6: (100.0, 35.0), 8: (125.0, 50.0)}
  for b in range(1, 10):
    pd, qd = loads.get(b, (0.0, 0.0))
    btype = 3 if b == 1 else (2 if b in (2, 3) else 1)
    bus.append([b, btype, pd, qd, 0.0, 0.0, 1, 1.0, 0.0, kv[b], 1, 1.1, 0.9])
  gen = [[1, 0.0, 0.0, 300.0, -300.0, 1.04, BASE, 1, 250.0, 10.0],
         [2, 150.0, 0.0, 120.0, -60.0, 1.025, BASE, 1, 300.0, 10.0],
         [3, 90.0, 0.0, 80.0, -40.0, 1.025, BASE, 1, 270.0, 10.0]]
  # transformers with off-nominal taps so the DC reference covers taps
  branch = [
    [1, 4, 0.0, 0.0600, 0.0, 250.0, 250.0, 250.0, 1.0, 0.0, 1, -360.0, 360.0],
    [2, 7, 0.0, 0.0650, 0.0, 300.0, 300.0, 300.0, 0.98, 0.0, 1, -360.0, 360.0],
    [3, 9, 0.0, 0.0600, 0.0, 300.0, 300.0, 300.0, 1.025, 0.0, 1, -360.0, 360.0],
    [4, 5, 0.0100, 0.0850, 0.176, 250.0, 250.0, 250.0, 0.0, 0.0, 1, -360.0, 360.0],
    [4, 6, 0.0170, 0.0920, 0.158, 250.0, 250.0, 250.0, 0.0, 0.0, 1, -360.0, 360.0],
    [5, 7, 0.0320, 0.1610, 0.306, 250.0, 250.0, 250.0, 0.0, 0.0, 1, -360.0, 360.0],
    [6, 9, 0.0390, 0.1700, 0.358, 150.0, 150.0, 150.0, 0.0, 0.0, 1, -360.0, 360.0],
    [7, 8, 0.0085, 0.0720, 0.149, 250.0, 250.0, 250.0, 0.0, 0.0, 1, -360.0, 360.0],
    [8, 9, 0.0119, 0.1008, 0.209, 150.0, 150.0, 150.0, 0.0, 0.0, 1, -360.0, 360.0],
  ]
  case = Case("sp_case9", bus, gen, branch, [f"{n}_{int(kv[i+1]) if kv[i+1] >= 100 else kv[i+1]}".replace(".", "p") for i, n in enumerate(names)], 1, [])
  case.kinds = ["T", "T", "T", "L", "L", "L", "L", "L", "L"]
  return case

def case9():
  case = build_case9()
  res = newton(case)
  assert res is not None and 0.9 <= res["V"].min() and res["V"].max() <= 1.1
  th, flows, pslack = dc_powerflow(case)
  header = [
    "SP_CASE9  Nine-bus system in the shape of the classic three-machine test",
    "system: three generator buses behind step-up transformers (two of them on",
    "off-nominal taps), a 230 kV ring with three loads. Self-built, no external",
    "data, no external license. Generated by tools/gen_sp_cases.py.",
    "",
    f"AC base case (Newton-Raphson, {res['iters']} iterations): Vmin {res['V'].min():.4f} pu at bus",
    f"{int(res['V'].argmin()) + 1}, Vmax {res['V'].max():.4f} pu, slack P {res['pslack']:.3f} MW, highest branch",
    f"loading {res['loads'].max():.1f} %.",
    "",
    "DC power-flow reference (independent implementation in the generator,",
    "MATPOWER conventions: B' from x and tap, Pfinj from the shift, slack bus 1).",
    "Used by test/test_dc_powerflow.jl as the DC anchor. Bus angles in degrees:",
  ]
  header.append("  " + ", ".join(f"{math.degrees(t):.6f}" for t in th))
  header.append("From-side active flows in MW, branch order as below:")
  header.append("  " + ", ".join(f"({int(r[0])},{int(r[1])}) {f:.4f}" for r, f in zip(case.branch, flows)))
  header.append(f"Slack generation {pslack:.4f} MW.")
  return case, header, dict(th=th, flows=flows, pslack=pslack, res=res)

# ---------------------------------------------------------------- main -----
def make_300():
  return build_synthetic("sp_case300", N=300, n_areas=3, n_ehv=24, n_trafo=14, n_branches=411,
                         n_gen=69, n_load=201, total_load_mw=9500.0, seed=300, hv_kv=220.0, ehv_kv=380.0, ehv_share_mw=0.45)

def make_1354():
  return build_synthetic("sp_case1354", N=1354, n_areas=6, n_ehv=90, n_trafo=54, n_branches=2040,
                         n_gen=260, n_load=900, total_load_mw=30000.0, seed=1354, hv_kv=220.0, ehv_kv=380.0)

def header_for(case, s):
  return [
    f"{case.name.upper()}  Synthetic {case.N}-bus operated grid, generated by tools/gen_sp_cases.py",
    f"(seed in the script). Two voltage levels ({int(case.bus[0][9])} kV and the backbone at",
    f"{int(max(r[9] for r in case.bus))} kV), {len(set(r[6] for r in case.bus))} areas, {s['n_branch']} branches, {s['n_gen']} generator buses",
    f"({sum(1 for g in case.gen if g[1] == 0.0 and int(g[0]) != case.slack)} synchronous condensers), {sum(1 for r in case.bus if r[2] > 0)} loads with",
    f"{s['load_mw']:.0f} MW / {s['load_mvar']:.0f} Mvar. Self-built, no external data, no external license.",
    "",
    f"Validated in the generator: base case {s['iters']} NR iterations, Vmin {s['vmin']:.4f} pu, Vmax",
    f"{s['vmax']:.4f} pu, highest branch loading {s['loading']:.1f} %, {s['clamped']} generators at a reactive",
    f"limit, slack (bus {case.slack}) {s['pslack']:.1f} MW. N-1 over {s['n1_checked']} of {s['n_branch']} branches:",
    f"{s['n1_island']} outages island a bus, {s['n1_fail']} do not converge. Ratings are sized per",
    "branch to the base-case flow with 30 percent headroom in the usual classes.",
  ]

if __name__ == "__main__":
  root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  outdir = os.path.join(root, "data", "mpower")
  which = sys.argv[1] if len(sys.argv) > 1 else "all"
  explicit = sys.argv[2] if len(sys.argv) > 2 else None
  if which in ("all", "9"):
    case, header, dc = case9()
    path = explicit or os.path.join(outdir, "sp_case9.m")
    write_case(case, path, header)
    print("sp_case9:", f"AC {dc['res']['iters']} it, Vmin {dc['res']['V'].min():.4f}, slack {dc['res']['pslack']:.2f} MW;",
          "DC angles", [round(math.degrees(t), 4) for t in dc["th"]], "slack", round(dc["pslack"], 4), "->", path)
  if which in ("all", "300"):
    case = make_300()
    s, _ = tune_and_validate(case, n1_full=True)
    path = explicit or os.path.join(outdir, "sp_case300.m")
    write_case(case, path, header_for(case, s))
    print("sp_case300:", s, "->", path)
  if which in ("all", "1354"):
    case = make_1354()
    s, _ = tune_and_validate(case, n1_full=False, n1_sample=60)
    path = explicit or os.path.join(outdir, "sp_case1354.m")
    write_case(case, path, header_for(case, s))
    print("sp_case1354:", s, "->", path)