# MTB virulence WGS pipeline

Frozen analysis code (release 1.0.0, fingerprint 74d12c396d65013d) for a whole-genome sequencing study of virulence genes in 55 *Mycobacterium tuberculosis* isolates from pulmonary TB, mapped to H37Rv NC_000962.3.

- `code/` per-sample pipeline (`run_sample.sh`, `mod_*.py`, `report.py`, `check.py`) and cohort analyses (`c1_` to `c6_`)
- `ref/panel/` the frozen 946-gene tiered panel (BED, 0-based half-open)
- `env/` tool and database versions, Python packages, pipeline fingerprint
- `RESTORE.md` how to rebuild the environment and verify the panel

No sequence data are included. Licence: MIT. Cite via `CITATION.cff` or the Zenodo DOI.
