#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re, io, gzip, argparse
from pathlib import Path
import pandas as pd

# ---------- safe CSV loader ----------
def safe_read_csv(path, **kwargs):
    base_kwargs = dict(low_memory=False, dtype=str)
    base_kwargs.update(kwargs)
    try:
        return pd.read_csv(path, **base_kwargs)
    except Exception:
        pass
    try:
        return pd.read_csv(path, engine="python", on_bad_lines="skip", **base_kwargs)
    except Exception:
        pass
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        txt = f.read()
    txt = txt.replace("“", '"').replace("”", '"').replace("’", "'").replace("\x00", "")
    txt = re.sub(r'"\n', '\n', txt)
    buf = io.StringIO(txt)
    return pd.read_csv(buf, engine="python", on_bad_lines="skip", **base_kwargs)

# ---------- binsummary ----------
def _normalize_columns(df):
    cols = pd.Index(df.columns.astype(str))
    cols = (cols.str.replace(r"^\ufeff","",regex=True)
                 .str.strip().str.lower()
                 .str.replace(r"[ \t\-]+","_",regex=True))
    df = df.copy(); df.columns = cols
    return df

def _pick_col(df, candidates, name):
    for c in candidates:
        if c in df.columns: return c
    raise KeyError(f"Required column not found for '{name}'. "
                   f"Available: {list(df.columns)[:20]}")

def read_binsummary(binsummary_path: Path):
    df = safe_read_csv(binsummary_path, sep="\t")
    df = _normalize_columns(df)
    sample_col = _pick_col(df, ["sample_id","sample","sampleid","sid","samp_id"], "sample_id")
    bin_col    = _pick_col(df, ["bin_id","bin","bid","binid"], "bin_id")
    comp_col   = _pick_col(df, ["completeness","complete","comp"], "completeness")
    out = df[[sample_col, bin_col, comp_col]].copy()
    out.columns = ["sample_id","bin_id","completeness"]
    out["sample_id"] = out["sample_id"].astype(str)
    out["bin_id"]    = out["bin_id"].astype(str)
    out["completeness"] = pd.to_numeric(out["completeness"], errors="coerce").fillna(0.0) / 100.0
    return out

# ---------- Roary (optional; for logging only) ----------
def parse_roary_gene_presence(roary_gpa_csv: Path):
    if not roary_gpa_csv.exists(): return None
    df = safe_read_csv(roary_gpa_csv)
    if "Gene" not in df.columns: return None
    meta = {"Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences",
            "Avg sequences per isolate","Genome Fragment","Order within Fragment",
            "Accessory Fragment","Accessory Order with Fragment","QC",
            "Min group size nuc","Max group size nuc","Avg group size nuc"}
    mag_cols = [c for c in df.columns if str(c) not in meta]
    _ = df[mag_cols].notna().astype(int)  # sanity only
    return True

# ---------- eMapper ----------
def _open_text_any(path: Path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")

def parse_emapper_annotations_relaxed(path: Path):
    out = {}
    with _open_text_any(path) as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if not parts: continue
            gene = parts[0]
            Ks = set(re.findall(r"\bK\d{5}\b", line))
            Ms = set(re.findall(r"\bM\d{5}\b", line))
            out[gene] = {"kegg": sorted(Ks), "keggmodule": sorted(Ms)}
    return out

# ---------- path helpers ----------
def make_bin_name(sample_id: str, bin_id: str) -> str:
    s = str(sample_id); b = str(bin_id)
    if re.match(r".+C\d+$", b, flags=re.IGNORECASE): return b
    if re.match(r"^[Cc]\d+$", b): return f"{s}{b.upper()}"
    if re.match(r"^\d+$", b):    return f"{s}C{b}"
    return f"{s}C{b}"

def find_emapper_path(annotation_base: Path, sample_id: str, bin_name: str):
    base = annotation_base / sample_id / f"{bin_name}.em"
    for cand in [
        base / f"{bin_name}.annotations",
        base / f"{bin_name}.annotations.gz",
        base / f"{bin_name}.emapper.annotations",
        base / f"{bin_name}.emapper.annotations.gz",
        base / f"{bin_name}.tsv",
        base / f"{bin_name}.tsv.gz",
    ]:
        if cand.exists(): return cand
    return None

def collect_bin_feature_sets_from_HaLD(bin_df: pd.DataFrame,
                                       annotation_base="/home/caefs/microbiome/projects/metagenome_analysis/annotation_HaLD"):
    bin2ko, bin2mod = {}, {}
    base = Path(annotation_base)
    hit, miss = 0, 0
    for _, row in bin_df.iterrows():
        sid = str(row["sample_id"]); bid = str(row["bin_id"])
        bname = make_bin_name(sid, bid)
        em_path = find_emapper_path(base, sid, bname)
        kos, mods = set(), set()
        if em_path is not None:
            em = parse_emapper_annotations_relaxed(em_path)
            for _, feats in em.items():
                kos.update(feats["kegg"]); mods.update(feats["keggmodule"])
            hit += 1
        else:
            print(f"[WARN] eMapper not found: {sid}/{bid} -> {bname}")
            miss += 1
        bin2ko[bname]  = kos
        bin2mod[bname] = mods
    print(f"[INFO] eMapper bins hit={hit}, miss={miss}")
    return bin2ko, bin2mod

# ---------- presence calculators ----------
def weighted_presence_per_asv(bin_df: pd.DataFrame, bin2set: dict):
    if not bin2set: return {}
    all_feats = set().union(*[s for s in bin2set.values() if s]) if len(bin2set) else set()
    if not all_feats: return {}
    weights = {}
    for _, r in bin_df.iterrows():
        bname = make_bin_name(str(r["sample_id"]), str(r["bin_id"]))
        weights[bname] = float(r["completeness"])
    denom = sum(weights.values())
    if denom <= 0: return {f: 0.0 for f in all_feats}
    out = {}
    for f in all_feats:
        num = sum(weights.get(b,0.0) for b, feats in bin2set.items() if f in feats)
        out[f] = num / denom
    return out

def binary_presence_per_asv(bin_df: pd.DataFrame, bin2set: dict):
    """Any bin has feature -> 1, else 0."""
    if not bin2set: return {}
    all_feats = set().union(*[s for s in bin2set.values() if s]) if len(bin2set) else set()
    if not all_feats: return {}
    out = {}
    for f in all_feats:
        out[f] = 1 if any(f in feats for feats in bin2set.values()) else 0
    return out

# ---------- per-ASV processing ----------
def process_one_asv(asv_dir: Path, mode="weighted", log=True):
    asv_id = asv_dir.name
    binsum = asv_dir / "binsummary.txt"
    gpa    = asv_dir / "roary" / "gene_presence_absence.csv"

    if log:
        print(f"[ASV {asv_id}] binsummary path: {binsum}")
        print(f"[ASV {asv_id}] roary GPA path: {gpa}")

    if not binsum.exists():
        if log: print(f"[ASV {asv_id}] SKIP: binsummary.txt not found")
        return asv_id, {}, {}

    try:
        bin_df = read_binsummary(binsum)
    except Exception as e:
        if log: print(f"[ASV {asv_id}] SKIP: binsummary parse failed: {e}")
        return asv_id, {}, {}

    if gpa.exists():
        try:
            _ = parse_roary_gene_presence(gpa)
            if log: print(f"[ASV {asv_id}] roary GPA loaded")
        except Exception as e:
            if log: print(f"[ASV {asv_id}] WARN: roary GPA load failed: {e}")

    bin2ko, bin2mod = collect_bin_feature_sets_from_HaLD(bin_df)

    if mode == "weighted":
        ko_row  = weighted_presence_per_asv(bin_df, bin2ko) if bin2ko else {}
        mod_row = weighted_presence_per_asv(bin_df, bin2mod) if bin2mod else {}
    else:  # binary
        ko_row  = binary_presence_per_asv(bin_df, bin2ko) if bin2ko else {}
        mod_row = binary_presence_per_asv(bin_df, bin2mod) if bin2mod else {}
    return asv_id, ko_row, mod_row

# ---------- finalize ----------
def finalize(rows_dict, out_csv, label):
    if not rows_dict:
        print(f"[INFO] No ASV processed for {label}."); return
    all_feats = sorted(set().union(*[set(d.keys()) for d in rows_dict.values() if d]))
    df = pd.DataFrame.from_dict(rows_dict, orient="index")
    for k in all_feats:
        if k not in df.columns: df[k] = 0.0
    df = df.fillna(0.0)
    df.insert(0, "ASV", df.index)
    df.to_csv(out_csv, index=False)
    print(f"Saved {label}: {out_csv}  ({df.shape[0]} ASVs x {df.shape[1]-1} features)")

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--basedir", required=True, help="dir containing per-ASV folders")
    ap.add_argument("--outprefix", default=None, help="output prefix (default: basedir name)")
    ap.add_argument("--annotation_base", default="/home/caefs/microbiome/projects/metagenome_analysis/annotation_HaLD",
                    help="root of annotation_HaLD")
    ap.add_argument("--mode", choices=["weighted","binary"], default="weighted",
                    help="weighted: completeness-weighted (0..1); binary: 1/0 presence")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    base = Path(args.basedir).resolve()
    outprefix = args.outprefix if args.outprefix else base.name
    log = (not args.quiet)

    rows_ko, rows_mod = {}, {}
    for asv_dir in sorted([p for p in base.iterdir() if p.is_dir()]):
        if not (asv_dir / "binsummary.txt").exists(): continue
        asv_id, ko_row, mod_row = process_one_asv(asv_dir, mode=args.mode, log=log)
        rows_ko[asv_id]  = ko_row if ko_row is not None else {}
        rows_mod[asv_id] = mod_row if mod_row is not None else {}

    if args.mode == "weighted":
        ko_name  = f"asv_ko_weighted_{outprefix}.csv"
        mod_name = f"asv_keggmodule_weighted_{outprefix}.csv"
    else:
        ko_name  = f"asv_ko_presence_{outprefix}.csv"
        mod_name = f"asv_keggmodule_presence_{outprefix}.csv"

    finalize(rows_ko,  ko_name,  "KO")
    finalize(rows_mod, mod_name, "Module")

if __name__ == "__main__":
    main()
