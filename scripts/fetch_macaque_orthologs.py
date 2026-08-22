"""
Fetch macaque-human orthologs from Ensembl BioMart.

Queries the M. fascicularis gene_ensembl dataset with hsapiens_homolog_*
cross-reference attributes to produce data/macaque/biomart_macaque_human_orthologs.csv.

Reproducibility:
- Ensembl BioMart release 115 (accessed 2026-03-15 for the manuscript run).
- Re-running this script against the current BioMart release may return
  materially different ortholog assignments (Ensembl releases shift over time).
- The archived CSV at data/macaque/biomart_macaque_human_orthologs.csv is the
  canonical source for reproduction; re-run this script only if you need to
  refresh against a newer Ensembl release.

Downstream filtering: scripts/nhp_ortholog_assessment.py filters this raw
return to the 13,927-gene 1:1 macaque-human ortholog operating space reported
in S1 Text section 7 (macaque extension).
"""

from pathlib import Path

from pybiomart import Dataset


OUTPUT_PATH = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "macaque"
    / "biomart_macaque_human_orthologs.csv"
)


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    dataset = Dataset(
        name="mfascicularis_gene_ensembl",
        host="http://www.ensembl.org",
    )
    df = dataset.query(
        attributes=[
            "ensembl_gene_id",
            "external_gene_name",
            "hsapiens_homolog_associated_gene_name",
            "hsapiens_homolog_orthology_type",
        ],
    )
    df.to_csv(OUTPUT_PATH, index=False)
    print(f"Wrote {OUTPUT_PATH}: {len(df):,} rows")


if __name__ == "__main__":
    main()
