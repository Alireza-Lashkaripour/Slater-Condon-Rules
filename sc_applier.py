import re
import csv
from itertools import combinations
from collections import Counter, defaultdict

def parse_line(line):
    m = re.match(r"\s*D\s*[_-]?\s*([A-Za-z0-9]+)\s*=\s*\((.*)\)\s*$", line)
    if not m:
        return None
    label = str(m.group(1))
    parts = [p.strip() for p in m.group(2).split(',') if p.strip()]
    out = []
    for p in parts:
        mm = re.match(r"f_\{([^}]+)\}\s*([abAB])\s*$", p)
        if not mm: return None
        out.append((mm.group(1).strip(), mm.group(2).lower()))
    return label, out

def load_determinants(path):
    labels, dets = [], []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if not s: continue
            r = parse_line(s)
            if r is None: continue
            labels.append(r[0]); dets.append(r[1])
    return labels, dets

def format_coeff(c):
    if c == 0.5:  return "1/2"
    if c == -0.5: return "-1/2"
    if c == 1.0:  return "+1"
    if c == -1.0: return "-1"
    return ("%+g" % c)

def diagonal_terms(det):
    terms = []
    for i in range(len(det)):
        p, sp = det[i]
        terms.append(("E", [p], [p], [(sp, sp)], 1.0))
    for i, j in combinations(range(len(det)), 2):
        pi, si = det[i]; pj, sj = det[j]
        terms.append(("J", [pi, pj], [pi, pj], [(si, si), (sj, sj)], 0.5))
        terms.append(("K", [pi, pj], [pi, pj], [(si, sj), (sj, si)], -0.5))
    return terms

def single_excitation_terms_orbs(ref, hole, part, phase):
    terms = []
    pi, si = hole; pj, sj = part
    terms.append(("E", [pi], [pj], [(si, sj)], float(phase)))
    for pk, sk in ref:
        if (pk, sk) == (pi, si): continue
        terms.append(("J", [pi, pk], [pj, pk], [(si, sj), (sk, sk)],  float(phase)))
        terms.append(("K", [pi, pk], [pj, pk], [(si, sk), (sk, sj)], -float(phase)))
    return terms

def double_excitation_terms_orbs(holes, parts_ordered, phase):
    (p1, s1), (p2, s2) = holes
    (qA, tA), (qB, tB) = parts_ordered
    out = []
    if (s1 == tA) and (s2 == tB):
        out.append(("J", [p1, p2], [qA, qB], [(s1, tA), (s2, tB)], float(phase)))
    if (s1 == tB) and (s2 == tA):
        out.append(("K", [p1, p2], [qA, qB], [(s1, tB), (s2, tA)], -float(phase)))
    return out

def align_exact_transpositions(ref, det):
    det = list(det)
    swaps = 0
    phase = 1
    n = min(len(ref), len(det))
    for i in range(n-1, -1, -1):
        if det[i] == ref[i]:
            continue
        want = ref[i]
        found = None
        for j in range(n):
            if j == i:
                continue
            if det[j] == want:
                found = j
                break
        if found is None:
            continue
        det[i], det[found] = det[found], det[i]
        swaps += 1
        phase *= -1
    return det, swaps, phase

def _display_target_with_spin_relabel(ref, ket):
    from collections import defaultdict, Counter
    avail = defaultdict(Counter)
    for sp, s in ket:
        avail[sp][s] += 1
    want = defaultdict(Counter)
    for sp, s in ref:
        want[sp][s] += 1
    assigned = defaultdict(list)
    for sp in avail:
        supply = avail[sp].copy()
        for s in ("a", "b"):
            use = min(supply[s], want[sp][s])
            if use > 0:
                assigned[sp].extend([s]*use)
                supply[s] -= use
                want[sp][s] -= use
        for s in ("a", "b"):
            if supply[s] > 0:
                assigned[sp].extend([s]*supply[s])
                supply[s] = 0
    next_pick = {sp: Counter(assigned[sp]) for sp in assigned}
    ket_display = []
    for sp, _s in ket:
        for s_try in ("a", "b"):
            if next_pick.get(sp, Counter())[s_try] > 0:
                ket_display.append((sp, s_try))
                next_pick[sp][s_try] -= 1
                break
        else:
            ket_display.append((sp, _s))
    return ket_display

def _pair_score(h, p):
    return (100 if h[0]==p[0] else 0) + (10 if h[1]==p[1] else 0)

def _choose_parts_order_for_degree2(ref, ket, holes, parts):
    ket_disp = _display_target_with_spin_relabel(ref, ket)
    h1, h2 = holes
    pA, pB = parts
    s1 = _pair_score(h1, pA) + _pair_score(h2, pB)
    s2 = _pair_score(h1, pB) + _pair_score(h2, pA)
    if s2 > s1:
        return [pB, pA]
    return [pA, pB]

def map_summary_for_degree2(holes, parts_ordered):
    pairs = [(holes[0], parts_ordered[0]), (holes[1], parts_ordered[1])]
    spatial_all = all(h[0]==p[0] for h,p in pairs)
    spins_all   = all(h[1]==p[1] for h,p in pairs)
    if spatial_all and spins_all:
        tag = "direct"
    elif spatial_all and not spins_all:
        tag = "direct (spin-flip)"
    else:
        tag = "mixed"
    return " ".join([
        f"{tag}:",
        " ; ".join([f"f_{{{h[0]}}} {h[1]} -> f_{{{p[0]}}} {p[1]}" for h,p in pairs])
    ])

def analyze_pair(Li, Di, Lj, Dj):
    DjA, swaps, phase = align_exact_transpositions(Di, Dj)
    Ai = Counter(Di)
    Bj = Counter(DjA)
    holes, parts = [], []
    for orb, c in (Ai - Bj).items():
        holes.extend([orb]*c)
    for orb, c in (Bj - Ai).items():
        parts.extend([orb]*c)
    if len(holes) != len(parts):
        return {"pair": (Li, Lj), "degree": -1, "phase": 0, "swaps": swaps, "terms": [], "maps": ""}
    deg = len(holes)
    if   deg == 0:
        terms = diagonal_terms(Di)
        maps = ""
    elif deg == 1:
        terms = single_excitation_terms_orbs(Di, holes[0], parts[0], phase)
        maps = ""
    elif deg == 2:
        parts_ordered = _choose_parts_order_for_degree2(Di, DjA, holes, parts)
        terms = double_excitation_terms_orbs(holes, parts_ordered, phase)
        maps = map_summary_for_degree2(holes, parts_ordered)
    else:
        terms = []
        maps = ""
    return {"pair": (Li, Lj), "degree": deg, "phase": phase, "swaps": swaps, "terms": terms, "maps": maps}

def filter_terms(terms):
    kept = []
    for tag, L, R, spins, coeff in terms:
        if any(a != b for (a, b) in spins):
            continue
        kept.append((tag, L, R, spins, coeff))
    return kept

def format_term(tag, L, R, spins, coeff):
    c = format_coeff(coeff)
    if tag == "E":
        sL, tR = spins[0]
        return f"{c} h_{{{L[0]} {sL},{R[0]} {tR}}}"
    sL0, tR0 = spins[0] if spins else ("", "")
    sL1, tR1 = spins[1] if len(spins) > 1 else ("", "")
    if tag == "J":
        left  = f"{L[0]} {sL0} {L[1]} {sL1}"
        right = f"{R[0]} {tR0} {R[1]} {tR1}"
        return f"{c} J_{{{left}, {right}}}"
    if tag == "K":
        left  = f"{L[0]} {sL0} {L[1]} {sL1}"
        right = f"{R[1]} {tR0} {R[0]} {tR1}"
        return f"{c} K_{{{left}, {right}}}"
    return ""

def main():
    in_path = input("Enter path to determinants file (.txt or .csv): ").strip()
    out_path = input("Enter output csv filename: ").strip()
    labels, dets = load_determinants(in_path)
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(["pair","degree","phase","phase_changes","exprs","maps"])
        n=len(dets)
        for i in range(n):
            for j in range(n):
                res = analyze_pair(labels[i], dets[i], labels[j], dets[j])
                if res["degree"] == -1:
                    continue
                kept = filter_terms(res["terms"])
                exprs = " + ".join(format_term(t, L, R, S, c) for (t,L,R,S,c) in kept)
                w.writerow([f"<D_{labels[i]}|H|D_{labels[j]}>", res["degree"], res["phase"], res["swaps"], exprs, res.get("maps","")])

if __name__=="__main__":
    main()
