#!/bin/bash
# Recompute the per-sample pipeline fingerprint and compare it to the frozen value.
# Run before every batch. A mismatch means the cohort would split: STOP and re-run all.
P=/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2
python3 - <<'PY'
import hashlib,os,json
P="/tmp/claude-0/-home-claude/65ea7895-0684-55a5-8ec1-c52349b315bb/scratchpad/p2"
meta=json.load(open(P+"/ref/pipeline_version.json"))
h=hashlib.sha256()
for f in sorted(meta["files"]):
    p=os.path.join(P,f)
    if not os.path.exists(p): print("FAIL missing:",f); raise SystemExit(2)
    h.update(f.encode())
    with open(p,'rb') as fh:
        for c in iter(lambda: fh.read(1<<20), b''): h.update(c)
now=h.hexdigest()[:16]
if now==meta["fingerprint"]:
    print("PIPELINE FROZEN AND UNCHANGED  (%s)"%now)
else:
    print("*** PIPELINE CHANGED ***  frozen=%s  now=%s"%(meta["fingerprint"],now))
    print("Every completed sample must be re-run, or the cohort is not comparable.")
    raise SystemExit(1)
PY
