#!/usr/bin/env python3
"""
T3-E Step 1: ENCODE ATAC-seq Feasibility Check (metadata only)

Queries the ENCODE REST API for ATAC-seq experiment metadata across 10 target
cell types in both human and mouse. Assesses whether sufficient publicly available
chromatin accessibility data exists to test the hypothesis that promoter-level
chromatin accessibility conservation predicts Procrustes rigidity.

Biology: After 7 mechanistic nulls (housekeeping, TF complexity, niche adaptation,
within-type variance, inter-donor variance, expression-level confounds, PPI centrality),
chromatin architecture is the surviving hypothesis. This script determines whether
ENCODE has the data to test it.

Math: The downstream analysis will compute a Spearman correlation between chromatin
conservation scores (Jaccard similarity of open promoters) and Procrustes rigidity
residuals across n cell types. At n=7, requires |rho|>=0.786 for p<0.05. Higher n
gives more statistical power.

NO peak files are downloaded. This is metadata retrieval only.
"""

import json
import os
import time
import urllib.request
import urllib.parse
import urllib.error
from pathlib import Path
from datetime import datetime

# ── Configuration ────────────────────────────────────────────────────────────

OUTPUT_DIR = Path("output/validation/t3e_feasibility")
CACHE_DIR = OUTPUT_DIR / "encode_raw_metadata"
REPORT_PATH = OUTPUT_DIR / "encode_feasibility_report.md"

ENCODE_BASE = "https://www.encodeproject.org"
REQUEST_TIMEOUT = 30  # seconds
MAX_RETRIES = 3
BACKOFF_DELAYS = [5, 15, 45]  # seconds between retries

# Target cell types with ENCODE-compatible search terms and synonyms
# Each entry: (our_label, short_key, [list of ENCODE term names to try])
CELL_TYPES = [
    ("CD8+ T cell", "cd8t", [
        "CD8-positive, alpha-beta T cell",
        "cytotoxic T cell",
        "CD8+ T cell",
        "CD8-positive T cell",
    ]),
    ("CD4+ T cell", "cd4t", [
        "CD4-positive, alpha-beta T cell",
        "CD4+ T cell",
        "CD4-positive T cell",
        "helper T cell",
        "naive thymus-derived CD4-positive, alpha-beta T cell",
    ]),
    ("B cell", "bcell", [
        "B cell",
        "naive B cell",
        "B lymphocyte",
        "B-cell",
    ]),
    ("NK cell", "nk", [
        "natural killer cell",
        "NK cell",
    ]),
    ("Monocyte", "monocyte", [
        "monocyte",
        "CD14-positive monocyte",
        "CD14+ monocyte",
        "classical monocyte",
        "CD14-positive, CD16-negative classical monocyte",
    ]),
    ("Macrophage", "macrophage", [
        "macrophage",
        "monocyte-derived macrophage",
    ]),
    ("Endothelial cell", "endothelial", [
        "endothelial cell",
        "endothelial cell of umbilical vein",  # HUVEC — flag separately
        "lung microvascular endothelial cell",
        "pulmonary artery endothelial cell",
        "cardiac endothelial cell",
        "dermis microvascular blood vessel endothelial cell",
    ]),
    ("Hepatocyte", "hepatocyte", [
        "hepatocyte",
        "liver cell",
    ]),
    ("Neutrophil", "neutrophil", [
        "neutrophil",
        "granulocyte",
        "polymorphonuclear leukocyte",
    ]),
    ("Plasma cell", "plasma", [
        "plasma cell",
        "plasmablast",
    ]),
]

ORGANISMS = [
    ("Homo sapiens", "human"),
    ("Mus musculus", "mouse"),
]

# Cell line biosample terms to exclude
CELL_LINE_TERMS = {
    "K562", "GM12878", "Jurkat", "HEK293", "HEK293T", "HeLa",
    "A549", "MCF-7", "HepG2", "Hep G2", "THP-1", "U937",
}

# Treatment keywords that indicate stimulation
STIM_KEYWORDS = [
    "LPS", "PMA", "anti-CD3", "anti-CD28", "IFN", "interferon",
    "lipopolysaccharide", "phorbol", "ionomycin", "IL-2", "IL-4",
    "IL-6", "TNF", "activation", "stimulat", "polariz",
]

# HUVEC-related terms to flag
HUVEC_TERMS = ["HUVEC", "umbilical vein", "human umbilical"]

FAILED_QUERIES = []


def encode_request(url, description=""):
    """
    Make an HTTP request to ENCODE with timeout and retry logic.

    Returns parsed JSON on success, None on failure after all retries.
    """
    for attempt in range(MAX_RETRIES):
        try:
            headers = {"Accept": "application/json"}
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as response:
                data = json.loads(response.read().decode("utf-8"))
                return data
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError,
                ConnectionError, OSError) as e:
            delay = BACKOFF_DELAYS[attempt] if attempt < len(BACKOFF_DELAYS) else 45
            print(f"  [Retry {attempt+1}/{MAX_RETRIES}] {description}: {e}")
            if attempt < MAX_RETRIES - 1:
                print(f"  Waiting {delay}s before retry...")
                time.sleep(delay)
            else:
                print(f"  FAILED after {MAX_RETRIES} attempts: {description}")
                FAILED_QUERIES.append({
                    "description": description,
                    "url": url,
                    "error": str(e),
                    "error_type": type(e).__name__,
                })
                return None
    return None


def query_encode_experiments(organism, term_name):
    """
    Query ENCODE for ATAC-seq experiments matching organism and biosample term.

    Returns list of experiment dicts from the @graph field, or empty list.
    """
    params = {
        "type": "Experiment",
        "assay_title": "ATAC-seq",
        "status": "released",
        "replicates.library.biosample.donor.organism.scientific_name": organism,
        "biosample_ontology.term_name": term_name,
        "format": "json",
        "limit": "all",
    }
    url = f"{ENCODE_BASE}/search/?{urllib.parse.urlencode(params, safe='+')}"
    desc = f"ATAC-seq {organism} '{term_name}'"

    data = encode_request(url, desc)
    if data is None:
        return []

    graph = data.get("@graph", [])
    return graph


def extract_experiment_info(exp):
    """
    Extract relevant metadata fields from an ENCODE experiment JSON object.
    """
    accession = exp.get("accession", "unknown")

    # Biosample info
    biosample_summary = exp.get("biosample_summary", "")
    biosample_ontology = exp.get("biosample_ontology", {})
    biosample_term = biosample_ontology.get("term_name", "unknown")
    biosample_type = biosample_ontology.get("classification", "unknown")

    # Replicate count
    replicates = exp.get("replicates", [])
    bio_rep_numbers = set()
    for rep in replicates:
        brn = rep.get("biological_replicate_number")
        if brn is not None:
            bio_rep_numbers.add(brn)
    n_bio_reps = len(bio_rep_numbers)

    # Files — check for peak file types
    files = exp.get("files", [])
    has_idr_peaks = False
    has_overlap_peaks = False
    peak_file_count = 0
    for f in files:
        if isinstance(f, str):
            # Sometimes files are just accession strings; can't assess type
            continue
        output_type = f.get("output_type", "") if isinstance(f, dict) else ""
        file_status = f.get("status", "") if isinstance(f, dict) else ""
        if file_status != "released":
            continue
        if "IDR thresholded peaks" in output_type:
            has_idr_peaks = True
            peak_file_count += 1
        elif "overlap" in output_type.lower() and "peak" in output_type.lower():
            has_overlap_peaks = True
            peak_file_count += 1

    # Treatment info
    treatments = []
    # Check experiment-level treatment
    exp_treatments = exp.get("treatments", [])
    for t in exp_treatments:
        if isinstance(t, dict):
            treatments.append(t.get("treatment_term_name", str(t)))
        elif isinstance(t, str):
            treatments.append(t)

    # Check biosample description for stimulation keywords
    description = exp.get("description", "")
    is_stimulated = False
    stim_details = []
    text_to_check = f"{biosample_summary} {description} {' '.join(treatments)}".lower()
    for kw in STIM_KEYWORDS:
        if kw.lower() in text_to_check:
            is_stimulated = True
            stim_details.append(kw)

    # Check if cell line
    is_cell_line = biosample_type == "cell line"
    for cl in CELL_LINE_TERMS:
        if cl.lower() in biosample_summary.lower() or cl.lower() in biosample_term.lower():
            is_cell_line = True

    # Check if HUVEC / cell line endothelial
    is_huvec = False
    for ht in HUVEC_TERMS:
        if ht.lower() in biosample_summary.lower() or ht.lower() in biosample_term.lower():
            is_huvec = True

    # Lab and date
    lab = exp.get("lab", {})
    lab_name = lab.get("title", "unknown") if isinstance(lab, dict) else str(lab)
    date_released = exp.get("date_released", "unknown")

    return {
        "accession": accession,
        "biosample_term": biosample_term,
        "biosample_type": biosample_type,
        "biosample_summary": biosample_summary,
        "n_bio_replicates": n_bio_reps,
        "has_idr_peaks": has_idr_peaks,
        "has_overlap_peaks": has_overlap_peaks,
        "peak_file_count": peak_file_count,
        "treatments": treatments,
        "is_stimulated": is_stimulated,
        "stim_details": stim_details,
        "is_cell_line": is_cell_line,
        "is_huvec": is_huvec,
        "description": description,
        "lab": lab_name,
        "date_released": date_released,
    }


def classify_experiment(info, our_cell_type):
    """
    Apply inclusion/exclusion/flag criteria to an experiment.

    Returns (status, reasons) where status is one of:
    - "include": passes all criteria
    - "flag": passes but has concerns
    - "exclude": fails criteria
    """
    reasons = []

    # EXCLUDE: cell lines
    if info["is_cell_line"]:
        return "exclude", ["cell line"]

    # EXCLUDE: stimulated/activated
    if info["is_stimulated"]:
        return "exclude", [f"stimulated ({', '.join(info['stim_details'])})"]

    # Start with include, accumulate flags
    status = "include"

    # FLAG: single replicate
    if info["n_bio_replicates"] < 2:
        status = "flag"
        reasons.append(f"single replicate (n={info['n_bio_replicates']})")

    # FLAG: no peak files visible in metadata
    if not info["has_idr_peaks"] and not info["has_overlap_peaks"]:
        # This is common — peak files may exist but not be embedded in search results
        # Flag but don't exclude; we'll check file endpoint separately if needed
        reasons.append("no peak files visible in experiment metadata (may exist — check files endpoint)")

    # FLAG: hepatocyte (QC-dependent per project design)
    if our_cell_type == "Hepatocyte":
        reasons.append("hepatocyte — QC-dependent per project design")

    # FLAG: macrophage activation state
    if our_cell_type == "Macrophage":
        summary_lower = info["biosample_summary"].lower()
        term_lower = info["biosample_term"].lower()
        if "m0" not in summary_lower and "m0" not in term_lower:
            if "m1" in summary_lower or "m2" in summary_lower:
                reasons.append(f"non-M0 macrophage (check: {info['biosample_summary'][:100]})")
            else:
                reasons.append(f"macrophage activation state unclear (check: {info['biosample_summary'][:100]})")

    # FLAG: HUVEC / endothelial cell line
    if info["is_huvec"]:
        status = "flag"
        reasons.append("HUVEC / cell line endothelial — report separately")

    # Check biosample_type for primary cell
    if info["biosample_type"] not in ("primary cell", "tissue"):
        if info["biosample_type"] == "in vitro differentiated cells":
            status = "flag"
            reasons.append(f"in vitro differentiated (biosample_type={info['biosample_type']})")
        elif info["biosample_type"] == "cell line":
            return "exclude", ["cell line (biosample_type)"]
        else:
            status = "flag"
            reasons.append(f"unusual biosample_type: {info['biosample_type']}")

    if status == "include" and reasons:
        status = "flag"

    return status, reasons


def query_cell_type(our_label, short_key, term_names, organism_scientific, organism_key):
    """
    Query ENCODE for a single cell type + organism combination, trying multiple
    term names. Returns dict with results and the successful term.
    """
    cache_file = CACHE_DIR / f"{short_key}_{organism_key}.json"

    # Check cache
    if cache_file.exists():
        print(f"  [CACHED] {our_label} / {organism_key}")
        with open(cache_file) as f:
            return json.load(f)

    print(f"  Querying {our_label} / {organism_key}...")

    all_experiments = []
    successful_terms = []
    tried_terms = []

    for term in term_names:
        tried_terms.append(term)
        exps = query_encode_experiments(organism_scientific, term)
        if exps:
            # Deduplicate by accession
            seen = {e["accession"] for e in all_experiments if "accession" in e}
            for exp in exps:
                if exp.get("accession") not in seen:
                    all_experiments.append(exp)
                    seen.add(exp.get("accession"))
            successful_terms.append(term)
            print(f"    '{term}' → {len(exps)} experiment(s)")
        else:
            print(f"    '{term}' → 0 experiments")
        # Small delay to be polite to ENCODE API
        time.sleep(0.5)

    # Extract and classify
    experiments_info = []
    for exp in all_experiments:
        info = extract_experiment_info(exp)
        status, reasons = classify_experiment(info, our_label)
        info["classification"] = status
        info["classification_reasons"] = reasons
        experiments_info.append(info)

    result = {
        "our_label": our_label,
        "organism": organism_key,
        "organism_scientific": organism_scientific,
        "tried_terms": tried_terms,
        "successful_terms": successful_terms,
        "total_experiments": len(experiments_info),
        "experiments": experiments_info,
        "query_timestamp": datetime.now().isoformat(),
    }

    # Save to cache
    with open(cache_file, "w") as f:
        json.dump(result, f, indent=2)

    return result


def query_encode_files_for_experiment(accession):
    """
    Query the ENCODE files endpoint for peak files associated with an experiment.
    Uses the experiment detail page which includes file info.
    """
    url = f"{ENCODE_BASE}/experiments/{accession}/?format=json"
    desc = f"Files for {accession}"
    data = encode_request(url, desc)
    if data is None:
        return {"has_idr": False, "has_overlap": False, "files": []}

    files = data.get("files", [])
    peak_files = []
    has_idr = False
    has_overlap = False

    for f in files:
        if isinstance(f, str):
            # It's a path reference, not embedded — skip
            continue
        if not isinstance(f, dict):
            continue
        output_type = f.get("output_type", "")
        file_format = f.get("file_format", "")
        status = f.get("status", "")

        if status != "released":
            continue

        if "IDR thresholded peaks" in output_type:
            has_idr = True
            peak_files.append({
                "accession": f.get("accession", ""),
                "output_type": output_type,
                "file_format": file_format,
                "assembly": f.get("assembly", ""),
            })
        elif "overlap" in output_type.lower() and "peak" in output_type.lower():
            has_overlap = True
            peak_files.append({
                "accession": f.get("accession", ""),
                "output_type": output_type,
                "file_format": file_format,
                "assembly": f.get("assembly", ""),
            })

    return {"has_idr": has_idr, "has_overlap": has_overlap, "files": peak_files}


def check_geo_metadata(accession, description):
    """
    Check GEO for dataset metadata using E-utilities API.
    Returns basic metadata dict or None.
    """
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gds&term={accession}&retmode=json"
    )
    data = encode_request(url, f"GEO search {accession}")
    if data is None:
        return None

    result = data.get("esearchresult", {})
    count = int(result.get("count", 0))
    id_list = result.get("idlist", [])

    geo_info = {
        "accession": accession,
        "description": description,
        "found": count > 0,
        "gds_ids": id_list,
    }

    # If found, get summary
    if id_list:
        summary_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gds&id={','.join(id_list[:3])}&retmode=json"
        )
        summary_data = encode_request(summary_url, f"GEO summary {accession}")
        if summary_data and "result" in summary_data:
            summaries = []
            for gid in id_list[:3]:
                s = summary_data["result"].get(gid, {})
                if s:
                    summaries.append({
                        "title": s.get("title", ""),
                        "summary": s.get("summary", "")[:500],
                        "taxon": s.get("taxon", ""),
                        "gdstype": s.get("gdstype", ""),
                        "n_samples": s.get("n_samples", ""),
                        "platform": s.get("gpl", ""),
                    })
            geo_info["summaries"] = summaries

    return geo_info


def generate_report(all_results, tier2_results):
    """
    Generate the structured feasibility report as Markdown.
    """
    lines = []
    lines.append("# T3-E ENCODE ATAC-seq Feasibility Report")
    lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("\n## Purpose")
    lines.append("Assess whether ENCODE has sufficient ATAC-seq data in matched")
    lines.append("human/mouse primary cell types to test the hypothesis that chromatin")
    lines.append("accessibility conservation at identity-gene promoters predicts")
    lines.append("Procrustes rigidity across cell types.")
    lines.append("")

    # ── Section 1: Per-Cell-Type Coverage Table ──────────────────────────────
    lines.append("## Section 1: Per-Cell-Type Coverage Table\n")

    for our_label, short_key, _ in CELL_TYPES:
        lines.append(f"### {our_label}\n")

        for _, org_key in ORGANISMS:
            key = f"{short_key}_{org_key}"
            result = all_results.get(key)

            if result is None:
                lines.append(f"**{org_key.capitalize()}:** QUERY FAILED — no data retrieved\n")
                continue

            exps = result.get("experiments", [])
            included = [e for e in exps if e["classification"] == "include"]
            flagged = [e for e in exps if e["classification"] == "flag"]
            excluded = [e for e in exps if e["classification"] == "exclude"]

            lines.append(f"**{org_key.capitalize()}:**")
            lines.append(f"- Experiments found: {len(exps)} total")
            lines.append(f"- Passing inclusion criteria: **{len(included)}**")
            lines.append(f"- Flagged (concerns): {len(flagged)}")
            lines.append(f"- Excluded: {len(excluded)}")

            if result.get("successful_terms"):
                lines.append(f"- ENCODE term(s) used: {', '.join(repr(t) for t in result['successful_terms'])}")
            else:
                lines.append(f"- ENCODE terms tried (none returned results): {', '.join(repr(t) for t in result.get('tried_terms', []))}")

            # List included experiments
            if included:
                lines.append(f"- Included experiments:")
                for e in included:
                    peak_info = ""
                    if e["has_idr_peaks"]:
                        peak_info = "IDR peaks available"
                    elif e["has_overlap_peaks"]:
                        peak_info = "overlap peaks only"
                    else:
                        peak_info = "peak files not visible in search metadata"
                    lines.append(f"  - **{e['accession']}**: {e['biosample_term']} "
                                f"({e['biosample_type']}), "
                                f"{e['n_bio_replicates']} bio reps, "
                                f"{peak_info}, "
                                f"lab: {e['lab']}")

            # List flagged experiments
            if flagged:
                lines.append(f"- Flagged experiments:")
                for e in flagged:
                    reasons_str = "; ".join(e["classification_reasons"])
                    lines.append(f"  - **{e['accession']}**: {e['biosample_term']} — "
                                f"FLAGS: {reasons_str}")

            # List excluded experiments
            if excluded:
                lines.append(f"- Excluded experiments:")
                for e in excluded:
                    reasons_str = "; ".join(e["classification_reasons"])
                    lines.append(f"  - {e['accession']}: {e['biosample_term']} — "
                                f"REASON: {reasons_str}")

            lines.append("")

    # ── Section 2: Recommended Mapping ───────────────────────────────────────
    lines.append("## Section 2: Recommended Cell Type → ENCODE Experiment Mapping\n")

    types_with_both = []
    types_missing = []

    for our_label, short_key, _ in CELL_TYPES:
        human_key = f"{short_key}_human"
        mouse_key = f"{short_key}_mouse"
        human_result = all_results.get(human_key)
        mouse_result = all_results.get(mouse_key)

        human_included = []
        mouse_included = []
        human_flagged = []
        mouse_flagged = []

        if human_result:
            human_included = [e for e in human_result.get("experiments", []) if e["classification"] == "include"]
            human_flagged = [e for e in human_result.get("experiments", []) if e["classification"] == "flag"]
        if mouse_result:
            mouse_included = [e for e in mouse_result.get("experiments", []) if e["classification"] == "include"]
            mouse_flagged = [e for e in mouse_result.get("experiments", []) if e["classification"] == "flag"]

        # Consider flagged experiments as available (flagged != excluded)
        human_available = human_included + human_flagged
        mouse_available = mouse_included + mouse_flagged

        has_both = len(human_available) > 0 and len(mouse_available) > 0

        if has_both:
            types_with_both.append(our_label)
            lines.append(f"### {our_label} — DATA AVAILABLE IN BOTH SPECIES\n")

            # Recommend best human experiment
            # Prefer included over flagged, then most replicates
            def rank_exp(e):
                return (
                    0 if e["classification"] == "include" else 1,
                    -e["n_bio_replicates"],
                    0 if e["has_idr_peaks"] else 1,
                )
            human_sorted = sorted(human_available, key=rank_exp)
            mouse_sorted = sorted(mouse_available, key=rank_exp)

            rec_h = human_sorted[0]
            rec_m = mouse_sorted[0]

            lines.append(f"**Recommended human:** {rec_h['accession']} — "
                        f"{rec_h['biosample_term']}, {rec_h['n_bio_replicates']} bio reps"
                        f"{' (IDR peaks)' if rec_h['has_idr_peaks'] else ''}")
            if len(human_sorted) > 1:
                alts = [f"{e['accession']} ({e['n_bio_replicates']} reps)" for e in human_sorted[1:3]]
                lines.append(f"  Alternatives: {', '.join(alts)}")

            lines.append(f"**Recommended mouse:** {rec_m['accession']} — "
                        f"{rec_m['biosample_term']}, {rec_m['n_bio_replicates']} bio reps"
                        f"{' (IDR peaks)' if rec_m['has_idr_peaks'] else ''}")
            if len(mouse_sorted) > 1:
                alts = [f"{e['accession']} ({e['n_bio_replicates']} reps)" for e in mouse_sorted[1:3]]
                lines.append(f"  Alternatives: {', '.join(alts)}")

            # Special flags for monocyte/macrophage
            if our_label == "Monocyte":
                lines.append("\n**Note (monocyte):** Should map to CD14+ monocyte primary cells. "
                            "Verify ENCODE biosample matches CD14-positive monocyte ontology.")
            elif our_label == "Macrophage":
                lines.append("\n**Note (macrophage):** Flag if monocyte-derived. Check activation "
                            "state — M0 preferred. Any M1/M2 polarized samples should be noted.")

            if rec_h.get("classification_reasons") or rec_m.get("classification_reasons"):
                lines.append("\n**Flags on recommended experiments:**")
                if rec_h.get("classification_reasons"):
                    lines.append(f"  Human: {'; '.join(rec_h['classification_reasons'])}")
                if rec_m.get("classification_reasons"):
                    lines.append(f"  Mouse: {'; '.join(rec_m['classification_reasons'])}")

            lines.append("")
        else:
            missing_species = []
            if not human_available:
                missing_species.append("human")
            if not mouse_available:
                missing_species.append("mouse")
            types_missing.append((our_label, missing_species))
            lines.append(f"### {our_label} — MISSING DATA ({', '.join(missing_species)})\n")
            if human_available:
                lines.append(f"  Human has {len(human_available)} experiment(s), but mouse has none.")
            elif mouse_available:
                lines.append(f"  Mouse has {len(mouse_available)} experiment(s), but human has none.")
            else:
                lines.append(f"  No qualifying experiments in either species.")
            lines.append("")

    # ── Section 3: Power Assessment ──────────────────────────────────────────
    lines.append("## Section 3: Power Assessment\n")

    n_available = len(types_with_both)
    lines.append(f"Cell types with valid ENCODE data in BOTH human and mouse: **n = {n_available}**\n")

    if types_with_both:
        lines.append(f"Available types: {', '.join(types_with_both)}\n")

    if types_missing:
        lines.append(f"Missing types: {', '.join(f'{t} (no {s})' for t, ss in types_missing for s in ss)}\n")

    thresholds = {
        5: 0.900,
        6: 0.829,
        7: 0.786,
        8: 0.738,
        9: 0.683,
        10: 0.648,
    }

    if n_available in thresholds:
        rho_needed = thresholds[n_available]
        lines.append(f"At n={n_available}: Spearman requires |ρ| ≥ {rho_needed} for p < 0.05 (two-tailed)\n")
    elif n_available > 10:
        lines.append(f"At n={n_available}: Spearman threshold < 0.648 (well-powered)\n")
    elif n_available < 5:
        lines.append(f"At n={n_available}: UNDERPOWERED — Spearman test not meaningful below n=5\n")

    # Check for failed queries and present scenarios
    if FAILED_QUERIES:
        n_failed_types = len(set(q["description"].split("'")[1] for q in FAILED_QUERIES if "'" in q["description"]))
        n_optimistic = n_available + n_failed_types
        lines.append(f"\n### Partial failure scenarios:")
        lines.append(f"- **Pessimistic** (failed queries have no data): n = {n_available}")
        if n_available in thresholds:
            lines.append(f"  → requires |ρ| ≥ {thresholds[n_available]}")
        lines.append(f"- **Optimistic** (failed queries match average): n = {n_optimistic}")
        if n_optimistic in thresholds:
            lines.append(f"  → requires |ρ| ≥ {thresholds[n_optimistic]}")
        elif n_optimistic > 10:
            lines.append(f"  → threshold < 0.648 (well-powered)")

    lines.append("")

    # Pre-registered thresholds from task description
    lines.append("### Pre-registered T3-E thresholds:")
    lines.append("- ρ ≥ 0.50: POSITIVE — chromatin conservation predicts rigidity")
    lines.append("- ρ < 0.35: triggers 8th null closure")
    lines.append("")
    if n_available >= 7:
        lines.append(f"**Assessment:** At n={n_available}, the pre-registered positive "
                     f"threshold (ρ≥0.50) is BELOW the statistical significance threshold "
                     f"(|ρ|≥{thresholds.get(n_available, '?')}). This means a positive "
                     f"result (ρ≥0.50) could still be non-significant at p<0.05. The "
                     f"null closure threshold (ρ<0.35) is testable — any ρ in that range "
                     f"would clearly fail significance.")
    elif n_available >= 5:
        lines.append(f"**Assessment:** At n={n_available}, MARGINAL statistical power. "
                     f"The Spearman significance threshold (|ρ|≥{thresholds.get(n_available, '?')}) "
                     f"is very high. Only extremely strong correlations would reach significance. "
                     f"Consider supplementing with Tier 2 sources to increase n.")
    else:
        lines.append(f"**Assessment:** At n={n_available}, INSUFFICIENT statistical power. "
                     f"Cannot run a meaningful Spearman test. Must supplement with Tier 2 "
                     f"sources or reconsider analysis design.")

    lines.append("")

    # ── Section 4: Tier 2 Source Summary ─────────────────────────────────────
    lines.append("## Section 4: Tier 2 Source Summary\n")

    if tier2_results.get("calderon"):
        geo = tier2_results["calderon"]
        lines.append("### Calderon et al. 2019 (GSE118189) — Human Immune ATAC-seq\n")
        if geo.get("found"):
            lines.append(f"- GEO record FOUND (GDS IDs: {', '.join(geo.get('gds_ids', []))})")
            if geo.get("summaries"):
                for s in geo["summaries"]:
                    lines.append(f"- Title: {s.get('title', 'N/A')}")
                    lines.append(f"- Taxon: {s.get('taxon', 'N/A')}")
                    lines.append(f"- Samples: {s.get('n_samples', 'N/A')}")
                    lines.append(f"- Summary: {s.get('summary', 'N/A')[:300]}")
            lines.append("\nExpected cell types available (from publication):")
            lines.append("- CD8+ T, CD4+ T, B cell, NK cell, monocyte, neutrophil (and others)")
            lines.append("- Processed peak files typically available on GEO")
            lines.append("- Human only — no mouse equivalent in this dataset")
        else:
            lines.append("- GEO record NOT FOUND or unavailable")
            lines.append("- Manual check required at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189")
        lines.append("")

    if tier2_results.get("immgen"):
        geo = tier2_results["immgen"]
        lines.append("### ImmGen ATAC (Yoshida et al. 2019, GSE131651) — Mouse Immune ATAC-seq\n")
        if geo.get("found"):
            lines.append(f"- GEO record FOUND (GDS IDs: {', '.join(geo.get('gds_ids', []))})")
            if geo.get("summaries"):
                for s in geo["summaries"]:
                    lines.append(f"- Title: {s.get('title', 'N/A')}")
                    lines.append(f"- Taxon: {s.get('taxon', 'N/A')}")
                    lines.append(f"- Samples: {s.get('n_samples', 'N/A')}")
                    lines.append(f"- Summary: {s.get('summary', 'N/A')[:300]}")
            lines.append("\nExpected cell types available (from publication):")
            lines.append("- CD8+ T, CD4+ T, B cell, NK cell, macrophage, neutrophil, and many others")
            lines.append("- Mouse only — pairs with Calderon for cross-species coverage")
            lines.append("- Processed peak files typically available on GEO")
        else:
            lines.append("- GEO record NOT FOUND or unavailable")
            lines.append("- Manual check required at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131651")
        lines.append("")

    lines.append("### Feasibility of combining Tier 1 + Tier 2\n")
    lines.append("Combining ENCODE (Tier 1) with Calderon/ImmGen (Tier 2) is feasible but ")
    lines.append("introduces potential systematic differences:")
    lines.append("- Different library preparation protocols (ENCODE standard vs lab-specific)")
    lines.append("- Different peak calling pipelines (ENCODE uniform pipeline vs published)")
    lines.append("- Different quality thresholds and filtering criteria")
    lines.append("")
    lines.append("**Mitigation:** For the Jaccard/overlap analysis planned in T3-E, re-call ")
    lines.append("peaks from aligned BAM files using a uniform pipeline would be ideal but ")
    lines.append("requires downloading raw data. Alternative: use published peak files and ")
    lines.append("assess batch effects by comparing overlapping cell types (e.g., CD8+ T in ")
    lines.append("both ENCODE and Calderon).")
    lines.append("")

    # ── Section 5: Blockers and Decisions Required ───────────────────────────
    lines.append("## Section 5: Blockers and Decisions Required\n")

    lines.append("### 5.1 Term name mappings resolved\n")
    for our_label, short_key, _ in CELL_TYPES:
        for _, org_key in ORGANISMS:
            key = f"{short_key}_{org_key}"
            result = all_results.get(key)
            if result and result.get("successful_terms"):
                lines.append(f"- {our_label} ({org_key}): matched via "
                            f"{', '.join(repr(t) for t in result['successful_terms'])}")
            elif result:
                lines.append(f"- {our_label} ({org_key}): NO MATCH FOUND "
                            f"(tried: {', '.join(repr(t) for t in result.get('tried_terms', []))})")
    lines.append("")

    lines.append("### 5.2 Decisions requiring human confirmation\n")

    decision_items = []

    # Flag any mapping ambiguities
    for our_label, short_key, _ in CELL_TYPES:
        for _, org_key in ORGANISMS:
            key = f"{short_key}_{org_key}"
            result = all_results.get(key)
            if result:
                exps = result.get("experiments", [])
                flagged = [e for e in exps if e["classification"] == "flag"]
                for e in flagged:
                    for reason in e.get("classification_reasons", []):
                        if "activation state" in reason or "HUVEC" in reason or "in vitro" in reason:
                            decision_items.append(
                                f"- **{our_label} ({org_key}) — {e['accession']}**: {reason}"
                            )

    # Flag protocol differences between human and mouse
    for our_label, short_key, _ in CELL_TYPES:
        h_result = all_results.get(f"{short_key}_human")
        m_result = all_results.get(f"{short_key}_mouse")
        h_included = [e for e in (h_result or {}).get("experiments", []) if e["classification"] in ("include", "flag")]
        m_included = [e for e in (m_result or {}).get("experiments", []) if e["classification"] in ("include", "flag")]

        if h_included and m_included:
            h_types = set(e["biosample_type"] for e in h_included)
            m_types = set(e["biosample_type"] for e in m_included)
            if h_types != m_types:
                decision_items.append(
                    f"- **{our_label}**: human biosample types {h_types} vs "
                    f"mouse biosample types {m_types} — protocol mismatch?"
                )

    if decision_items:
        for item in decision_items:
            lines.append(item)
    else:
        lines.append("- No ambiguous mappings requiring human confirmation.")

    lines.append("")

    # Failed queries
    lines.append("### 5.3 Failed queries\n")
    if FAILED_QUERIES:
        for fq in FAILED_QUERIES:
            lines.append(f"- **{fq['description']}**: {fq['error_type']} — {fq['error']}")
            lines.append(f"  URL: {fq['url'][:200]}")
    else:
        lines.append("- No queries failed. All cell type / organism combinations retrieved successfully.")
    lines.append("")

    # Overall recommendation
    lines.append("### 5.4 Overall recommendation\n")
    if n_available >= 7:
        lines.append(f"With n={n_available} cell types having matched human/mouse ENCODE data, ")
        lines.append(f"T3-E is **feasible** for a Spearman correlation test. Proceed to Step 2 ")
        lines.append(f"(peak file download) after advisor confirms cell type selections.")
    elif n_available >= 5:
        lines.append(f"With n={n_available} cell types from ENCODE alone, T3-E is **marginal**. ")
        lines.append(f"Consider supplementing with Tier 2 sources (Calderon + ImmGen) to increase n. ")
        lines.append(f"Advisor should assess whether the statistical power is acceptable.")
    else:
        lines.append(f"With n={n_available} cell types from ENCODE, T3-E is **underpowered**. ")
        lines.append(f"ENCODE alone is insufficient. Tier 2 sources are required. Advisor should ")
        lines.append(f"assess whether a mixed-source approach is methodologically acceptable.")
    lines.append("")

    return "\n".join(lines)


def main():
    print("=" * 70)
    print("T3-E Step 1: ENCODE ATAC-seq Feasibility Check")
    print("Metadata retrieval only — no peak file downloads")
    print("=" * 70)

    # ── Query ENCODE for all cell types ──────────────────────────────────────
    all_results = {}

    for our_label, short_key, term_names in CELL_TYPES:
        print(f"\n{'─' * 50}")
        print(f"Cell type: {our_label}")
        for organism_sci, org_key in ORGANISMS:
            result = query_cell_type(our_label, short_key, term_names,
                                     organism_sci, org_key)
            all_results[f"{short_key}_{org_key}"] = result

            # Summary
            exps = result.get("experiments", [])
            included = sum(1 for e in exps if e["classification"] == "include")
            flagged = sum(1 for e in exps if e["classification"] == "flag")
            excluded = sum(1 for e in exps if e["classification"] == "exclude")
            print(f"    → {org_key}: {included} included, {flagged} flagged, "
                  f"{excluded} excluded (total {len(exps)})")

    # ── For experiments that passed but lack peak file info, do a spot-check
    #    on the experiment detail endpoint for a few key types ─────────────
    print(f"\n{'─' * 50}")
    print("Spot-checking experiment detail endpoints for peak file availability...")

    for our_label, short_key, _ in CELL_TYPES:
        for _, org_key in ORGANISMS:
            key = f"{short_key}_{org_key}"
            result = all_results.get(key)
            if not result:
                continue
            for exp_info in result.get("experiments", []):
                if exp_info["classification"] in ("include", "flag"):
                    if not exp_info["has_idr_peaks"] and not exp_info["has_overlap_peaks"]:
                        acc = exp_info["accession"]
                        print(f"  Checking {acc} ({our_label} / {org_key})...")
                        file_info = query_encode_files_for_experiment(acc)
                        if file_info["has_idr"]:
                            exp_info["has_idr_peaks"] = True
                            exp_info["peak_files_from_detail"] = file_info["files"]
                            print(f"    → IDR peaks found!")
                        elif file_info["has_overlap"]:
                            exp_info["has_overlap_peaks"] = True
                            exp_info["peak_files_from_detail"] = file_info["files"]
                            print(f"    → Overlap peaks found!")
                        else:
                            print(f"    → No peak files found in detail endpoint either")
                        time.sleep(0.5)

    # ── Check Tier 2 sources ─────────────────────────────────────────────────
    print(f"\n{'─' * 50}")
    print("Checking Tier 2 sources...")

    tier2_results = {}

    print("  Calderon et al. 2019 (GSE118189)...")
    tier2_results["calderon"] = check_geo_metadata(
        "GSE118189",
        "Calderon et al. 2019 — Human immune ATAC-seq across 19 cell types"
    )

    time.sleep(1)

    print("  ImmGen ATAC / Yoshida et al. (GSE131651)...")
    tier2_results["immgen"] = check_geo_metadata(
        "GSE131651",
        "Yoshida et al. 2019 — ImmGen mouse immune ATAC-seq"
    )

    # Save tier2 results
    with open(CACHE_DIR / "tier2_results.json", "w") as f:
        json.dump(tier2_results, f, indent=2)

    # ── Generate Report ──────────────────────────────────────────────────────
    print(f"\n{'─' * 50}")
    print("Generating feasibility report...")

    report = generate_report(all_results, tier2_results)

    with open(REPORT_PATH, "w") as f:
        f.write(report)

    print(f"\nReport saved to: {REPORT_PATH}")

    # ── Print Summary ────────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")

    n_types_with_both = 0
    for our_label, short_key, _ in CELL_TYPES:
        h = all_results.get(f"{short_key}_human", {})
        m = all_results.get(f"{short_key}_mouse", {})
        h_ok = any(e["classification"] in ("include", "flag")
                   for e in h.get("experiments", []))
        m_ok = any(e["classification"] in ("include", "flag")
                   for e in m.get("experiments", []))
        status = "BOTH" if h_ok and m_ok else "HUMAN ONLY" if h_ok else "MOUSE ONLY" if m_ok else "NEITHER"
        if h_ok and m_ok:
            n_types_with_both += 1
        print(f"  {our_label:25s} → {status}")

    print(f"\nCell types with matched human+mouse data: n = {n_types_with_both}")

    thresholds = {5: 0.900, 6: 0.829, 7: 0.786, 8: 0.738, 9: 0.683, 10: 0.648}
    if n_types_with_both in thresholds:
        print(f"Spearman significance threshold at n={n_types_with_both}: "
              f"|ρ| ≥ {thresholds[n_types_with_both]}")
    elif n_types_with_both > 10:
        print(f"Spearman significance threshold at n={n_types_with_both}: < 0.648 (well-powered)")

    if FAILED_QUERIES:
        print(f"\nWARNING: {len(FAILED_QUERIES)} queries failed after retries:")
        for fq in FAILED_QUERIES:
            print(f"  - {fq['description']}")

    print(f"\nFull report: {REPORT_PATH}")
    print(f"Raw metadata cached in: {CACHE_DIR}")


if __name__ == "__main__":
    main()
