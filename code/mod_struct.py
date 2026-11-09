"""
Analysis 25b - structural consequence of virulence-panel variants.

Pipeline per variant (gene Rv number + protein change p.XnnnY):
  1. map CDS residue number -> AlphaFold model residue number (validated offset)
  2. confirm the wild-type residue in the model matches the variant's WT residue
  3. pLDDT model confidence at that residue
  4. secondary structure  (Kabsch & Sander DSSP H-bond energy, implemented here)
  5. relative solvent accessibility (Shrake-Rupley, implemented here)
  6. distance to nearest UniProt active site / binding site / domain
  7. FoldX 5.1 RepairPDB + BuildModel  ->  ddG (kcal/mol)

Everything except FoldX is dependency-free: no gemmi, no Biopython, no DSSP binary.
"""
import os, gzip, json, math, subprocess, shutil, tempfile, re

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AF_DIR   = os.path.join(P, "ref", "af_pdb")
FOLDX    = os.path.join(P, "src", "foldx", "foldx_20261231")
FOLDX_MOL= os.path.join(P, "src", "foldx", "molecules")

_cache = {}
def _ref(name):
    if name not in _cache:
        _cache[name] = json.load(open(os.path.join(P, "ref", name)))
    return _cache[name]

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
AA1 = {v: k for k, v in AA3.items()}

# Tien et al. 2013, theoretical maximum solvent accessibility (A^2)
MAXASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,
          'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}

VDW = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'H':1.20}


# ---------------------------------------------------------------- PDB reading
def load_model(acc):
    """Parse an AlphaFold PDB. Returns dict with per-residue atoms and pLDDT."""
    path = os.path.join(AF_DIR, "AF-%s-F1-model_v6.pdb.gz" % acc)
    if not os.path.exists(path):
        return None
    res = {}
    order = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            name = line[12:16].strip()
            rn   = line[17:20].strip()
            num  = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            b    = float(line[60:66])
            el   = (line[76:78].strip() or name[0]).upper()
            if num not in res:
                res[num] = {"aa": AA3.get(rn, 'X'), "atoms": [], "plddt": b}
                order.append(num)
            res[num]["atoms"].append((name, x, y, z, el))
    return {"res": res, "order": order, "acc": acc, "path": path}


# ------------------------------------------------- secondary structure (DSSP)
def _hbond_energy(don, acc):
    """Kabsch & Sander electrostatic H-bond energy (kcal/mol)."""
    try:
        N, H = don["N"], don["H"]
        C, O = acc["C"], acc["O"]
    except KeyError:
        return 0.0
    def d(a, b):
        return math.dist(a, b)
    rON, rCH, rOH, rCN = d(O, N), d(C, H), d(O, H), d(C, N)
    if min(rON, rCH, rOH, rCN) < 0.5:
        return 0.0
    return 0.084 * (1/rON + 1/rCH - 1/rOH - 1/rCN) * 332.0


def secondary_structure(model):
    """DSSP H/E/C assignment from backbone geometry."""
    order = model["order"]
    bb = {}
    for i, num in enumerate(order):
        a = {n: (x, y, z) for n, x, y, z, e in model["res"][num]["atoms"]}
        if not {"N", "CA", "C", "O"} <= set(a):
            continue
        rec = {"N": a["N"], "CA": a["CA"], "C": a["C"], "O": a["O"]}
        # amide H placed along the bisector, standard DSSP approximation
        if i > 0:
            pn = order[i-1]
            pa = {n: (x, y, z) for n, x, y, z, e in model["res"][pn]["atoms"]}
            if "C" in pa and "O" in pa and model["res"][num]["aa"] != 'P':
                cx, cy, cz = pa["C"]; ox, oy, oz = pa["O"]
                dx, dy, dz = cx-ox, cy-oy, cz-oz
                L = math.sqrt(dx*dx+dy*dy+dz*dz) or 1.0
                nx, ny, nz = rec["N"]
                rec["H"] = (nx+dx/L, ny+dy/L, nz+dz/L)
        bb[num] = rec

    nums = [n for n in order if n in bb]
    pos  = {n: bb[n]["CA"] for n in nums}
    hb = set()
    for i, ni in enumerate(nums):
        for nj in nums:
            if abs(ni - nj) < 2:
                continue
            if math.dist(pos[ni], pos[nj]) > 9.0:
                continue
            if _hbond_energy(bb[ni], bb[nj]) < -0.5:
                hb.add((ni, nj))          # ni donates to nj

    ss = {n: 'C' for n in order}
    # alpha helix: i -> i-4 pattern on consecutive residues
    for n in nums:
        if (n+4 in bb) and ((n+4, n) in hb) and (n+5 in bb) and ((n+5, n+1) in hb):
            for k in range(n, n+5):
                if k in ss:
                    ss[k] = 'H'
    # beta: reciprocal or bridging H-bonds between distant residues
    for (a, b) in hb:
        if abs(a - b) < 5:
            continue
        if (b, a) in hb or (a+2, b-2) in hb or (b+2, a-2) in hb:
            for k in (a, b):
                if ss.get(k) == 'C':
                    ss[k] = 'E'
    return ss


# -------------------------------------------- solvent accessibility (Shrake-Rupley)
def _sphere(n=92):
    pts, off = [], 2.0/n
    inc = math.pi * (3.0 - math.sqrt(5.0))
    for k in range(n):
        y = k*off - 1.0 + off/2.0
        r = math.sqrt(max(0.0, 1.0 - y*y))
        phi = k*inc
        pts.append((math.cos(phi)*r, y, math.sin(phi)*r))
    return pts
_SPH = _sphere(92)


def solvent_accessibility(model, targets, probe=1.4):
    """Absolute and relative SASA for the requested residue numbers."""
    atoms = []
    for num in model["order"]:
        for name, x, y, z, el in model["res"][num]["atoms"]:
            atoms.append((num, x, y, z, VDW.get(el, 1.70)))
    out = {}
    for tnum in targets:
        if tnum not in model["res"]:
            continue
        total = 0.0
        for (num, x, y, z, r) in atoms:
            if num != tnum:
                continue
            R = r + probe
            near = [(x2, y2, z2, r2 + probe) for (n2, x2, y2, z2, r2) in atoms
                    if n2 != num or True
                    if (x2-x)**2 + (y2-y)**2 + (z2-z)**2 < (R + r2 + probe)**2
                    and not (x2 == x and y2 == y and z2 == z)]
            acc_pts = 0
            for (sx, sy, sz) in _SPH:
                px, py, pz = x + R*sx, y + R*sy, z + R*sz
                buried = False
                for (x2, y2, z2, R2) in near:
                    if (px-x2)**2 + (py-y2)**2 + (pz-z2)**2 < R2*R2:
                        buried = True
                        break
                if not buried:
                    acc_pts += 1
            total += 4.0*math.pi*R*R * acc_pts/len(_SPH)
        aa = model["res"][tnum]["aa"]
        rsa = 100.0*total/MAXASA.get(aa, 200)
        out[tnum] = {"sasa": round(total, 1), "rsa": round(rsa, 1),
                     "burial": "buried" if rsa < 20 else ("partially exposed" if rsa < 50 else "exposed")}
    return out


# ------------------------------------------------- functional site annotation
def _parse_sites(field, tag):
    """Pull residue positions out of a UniProt feature-table string."""
    out = []
    for m in re.finditer(tag + r'\s+(\d+)(?:\.\.(\d+))?', field or ""):
        a = int(m.group(1)); b = int(m.group(2) or m.group(1))
        out.extend(range(a, b+1))
    return sorted(set(out))


def functional_sites(acc):
    meta = _ref("af_index.json")["acc2meta"].get(acc, {})
    return {"active_site":  _parse_sites(meta.get("act", ""),  "ACT_SITE"),
            "binding_site": _parse_sites(meta.get("bind", ""), "BINDING"),
            "transmembrane":_parse_sites(meta.get("tm", ""),   "TRANSMEM"),
            "signal_peptide":_parse_sites(meta.get("sig", ""), "SIGNAL"),
            "domains": meta.get("dom", "")}


def nearest_site(model, resnum, sites):
    """Minimum CA-CA distance from resnum to each class of functional site."""
    def ca(n):
        if n not in model["res"]:
            return None
        for name, x, y, z, e in model["res"][n]["atoms"]:
            if name == "CA":
                return (x, y, z)
        return None
    p = ca(resnum)
    out = {}
    if p is None:
        return out
    for cls in ("active_site", "binding_site"):
        best = None
        for s in sites.get(cls, []):
            q = ca(s)
            if q is None:
                continue
            d = math.dist(p, q)
            if best is None or d < best[0]:
                best = (d, s)
        if best:
            out[cls] = {"distance_A": round(best[0], 1), "nearest_residue": best[1],
                        "contact": best[0] <= 8.0}
    return out


# ----------------------------------------------------------------- FoldX ddG
def foldx_ddg(model, muts, workdir, repair=True):
    """
    muts: list of (wt1, resnum_model, mt1)
    Returns {mutstring: ddG_kcal_per_mol}
    """
    if not os.path.exists(FOLDX):
        return {"error": "FoldX binary not available"}
    os.makedirs(workdir, exist_ok=True)
    if os.path.isdir(FOLDX_MOL) and not os.path.isdir(os.path.join(workdir, "molecules")):
        shutil.copytree(FOLDX_MOL, os.path.join(workdir, "molecules"))
    pdb = os.path.join(workdir, "model.pdb")
    with gzip.open(model["path"], 'rt') as src, open(pdb, 'w') as dst:
        for line in src:
            if line.startswith(("ATOM", "TER", "END")):
                dst.write(line)

    def run(args):
        return subprocess.run([FOLDX] + args, cwd=workdir,
                              capture_output=True, text=True, timeout=14400)

    target = "model.pdb"
    if repair:
        r = run(["--command=RepairPDB", "--pdb=model.pdb", "--output-dir=."])
        if os.path.exists(os.path.join(workdir, "model_Repair.pdb")):
            target = "model_Repair.pdb"

    listfile = os.path.join(workdir, "individual_list.txt")
    strings = ["%s%s%d%s" % (wt, "A", num, mt) for (wt, num, mt) in muts]
    with open(listfile, 'w') as fh:
        for s in strings:
            fh.write(s + ";\n")

    r = run(["--command=BuildModel", "--pdb=" + target,
             "--mutant-file=individual_list.txt", "--numberOfRuns=3", "--output-dir=."])

    ddg = {}
    fit = os.path.join(workdir, "Dif_" + target.replace(".pdb", "") + ".fxout")
    if os.path.exists(fit):
        vals = []
        for line in open(fit):
            p = line.split("\t")
            if len(p) > 1 and p[0].startswith(target.replace(".pdb", "")):
                try:
                    vals.append(float(p[1]))
                except ValueError:
                    pass
        # numberOfRuns=3 -> three consecutive rows per mutation
        per = 3
        for i, s in enumerate(strings):
            chunk = vals[i*per:(i+1)*per]
            if chunk:
                ddg[s] = {"ddG_kcal_mol": round(sum(chunk)/len(chunk), 2),
                          "runs": [round(v, 2) for v in chunk],
                          "sd": round((sum((v-sum(chunk)/len(chunk))**2 for v in chunk)/len(chunk))**0.5, 2)}
    if not ddg:
        ddg = {"error": "FoldX produced no output", "stderr": (r.stderr or "")[-400:]}
    return ddg


# ------------------------------------------------------------------ top level
def interpret(ddg, rsa, plddt, site):
    """Plain-language reading of the structural result."""
    bits = []
    if plddt is not None:
        bits.append("model confidence %s (pLDDT %.0f)" % (
            "very high" if plddt >= 90 else "confident" if plddt >= 70 else
            "low" if plddt >= 50 else "very low", plddt))
    if rsa is not None:
        bits.append("residue is %s (RSA %.0f%%)" % (rsa["burial"], rsa["rsa"]))
    for cls, lbl in (("active_site", "catalytic site"), ("binding_site", "ligand-binding site")):
        if cls in site:
            s = site[cls]
            bits.append("%.1f A from the nearest %s (residue %d)%s" % (
                s["distance_A"], lbl, s["nearest_residue"],
                " - within contact range" if s["contact"] else ""))
    if isinstance(ddg, dict) and "ddG_kcal_mol" in ddg:
        v = ddg["ddG_kcal_mol"]
        cat = ("strongly destabilising" if v > 2.0 else "destabilising" if v > 0.5 else
               "neutral" if v >= -0.5 else "stabilising")
        bits.append("predicted %s (ddG %+0.2f kcal/mol)" % (cat, v))
    return "; ".join(bits)


def analyse_variant(rv, wt, cds_resnum, mt, workdir=None, run_foldx=True):
    """
    rv          Rv number, e.g. 'Rv1908c'
    wt, mt      one-letter wild-type and mutant residues
    cds_resnum  residue number in the H37Rv CDS translation (1-based)
    """
    mp = _ref("rv2struct.json").get(rv)
    out = {"gene": rv, "variant": "p.%s%d%s" % (wt, cds_resnum, mt)}
    if not mp:
        out["status"] = "no AlphaFold model for this gene"
        return out
    acc = mp["acc"]
    model = load_model(acc)
    if model is None:
        out["status"] = "structure file missing"
        return out
    out["uniprot"] = acc
    out["mapping"] = {"method": mp["match"], "offset": mp["offset"],
                      "identity_pct": mp.get("identity"), "coverage_pct": mp.get("coverage")}

    # CDS numbering -> model numbering
    mres = cds_resnum - (mp["offset"] - 1)
    out["model_residue"] = mres
    if mres not in model["res"]:
        out["status"] = "residue %d falls outside the modelled region (model covers 1-%d)" % (
            mres, max(model["order"]))
        return out
    obs = model["res"][mres]["aa"]
    if obs != wt:
        out["status"] = ("wild-type mismatch: variant expects %s, model has %s at residue %d "
                         "- not analysed" % (wt, obs, mres))
        return out
    out["status"] = "ok"
    out["plddt"] = round(model["res"][mres]["plddt"], 1)
    ss = secondary_structure(model)
    out["secondary_structure"] = {"H": "alpha helix", "E": "beta strand", "C": "coil/loop"}[ss.get(mres, 'C')]
    out["solvent_accessibility"] = solvent_accessibility(model, [mres])[mres]
    sites = functional_sites(acc)
    out["functional_sites"] = nearest_site(model, mres, sites)
    if sites["transmembrane"]:
        out["transmembrane_residue"] = mres in sites["transmembrane"]
    ddg = {}
    if run_foldx:
        wd = workdir or tempfile.mkdtemp(prefix="foldx_%s_" % rv)
        out["foldx"] = foldx_ddg(model, [(wt, mres, mt)], wd)
        # foldx_ddg is keyed by FoldX mutation string; lift the single result out
        for k, v in out["foldx"].items():
            if isinstance(v, dict) and "ddG_kcal_mol" in v:
                ddg = v
                out["ddG_kcal_mol"] = v["ddG_kcal_mol"]
                break
    out["interpretation"] = interpret(ddg, out["solvent_accessibility"],
                                      out["plddt"], out["functional_sites"])
    return out
