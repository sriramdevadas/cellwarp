# SUBMISSION_RECORD

What was sent to PLOS ONE, as built and as compiled. Companion to `DEPOSIT_MANIFEST.md`, which
records the Zenodo side; this file records the journal side. Every value below was measured from
the artifact named beside it, not carried forward from a plan.

Source commit: `5e77414`. `git ls-files` at that commit: 1094 (this file is added after it).

## Identifiers

- **Editorial Manager submission id: `PONE-S-26-56365`.**
- **This is not the manuscript number.** The compiled PDF's page-1 cover sheet prints
  `Manuscript Number:` with an empty value and is titled `--Manuscript Draft--`. The assigned
  manuscript number arrives later and must be added here when it does.
- **Provenance caution.** The string `PONE-S-26-56365` does not appear anywhere inside the compiled
  PDF -- not in page text, not in link annotations, and not in the raw file bytes (`PONE` and
  `56365` are both absent). Its only local evidence is the downloaded filename
  `PONE-S-26-56365.pdf`. Recorded as an identifier supplied by Editorial Manager, not as one
  verifiable from the artifact.

## The compiled PDF

Editorial Manager builds this from the uploaded items; it is not something this repository produces.

```
file          PONE-S-26-56365.pdf
md5           173564d321a57140960559d337d61d7c
bytes         5,238,533
pages         92
producer      iTextSharp 5.5.13.2 (iText Group NV, AGPL version)
built         2026-08-24 15:38:50 -04:00  =  2026-08-24 19:38:50 UTC
encrypted     no        embedded file attachments  0
```

**Items: 28** = 26 attachment stubs (pages 67-92, one per uploaded figure / SI / cover item) plus the
manuscript body and the cover letter. Each stub carries an Editorial Manager download link; the five
figure pages carry two identical link annotations each, one distinct file id per figure:

```
Fig 1  45751340     Fig 2  45751342     Fig 3  45751343
Fig 4  45751344     Fig 5  45751345
```

**Supporting information: 20 distinct captions** in the compiled body -- S1-S5 Fig (5),
S1-S7 and S9-S14 Table (13, no S8), S1-S2 Text (2). These 20 caption *items* are carried by 21
*files*: S9 Table is two CSVs. See `docs/supplementary_materials/MANIFEST.md`.

**DOIs printed in the compiled PDF**, by page:

```
10.5281/zenodo.20735611   CODE concept    pages 4, 7, 46
10.5281/zenodo.20735639   DATA concept    pages 5, 7, 46
10.5281/zenodo.22073208   basal ganglia   pages 5, 47
```

The two concept DOIs resolve to the v2 records (CODE v2 `10.5281/zenodo.22083132`, DATA v2
`10.5281/zenodo.22083465`), which is why publishing v2 required no change to submitted text. The
basal-ganglia deposit prints its **version** DOI deliberately; the reasoning is in the commit
message of `1570187` and does not need repeating here.

## Figures as the compiled proof renders them

Editorial Manager rasterises each figure for the proof. Measured from the PDF:

| page | figure | embedded raster | mode | ppi | JPEG bytes | of its own raw |
|---|---|---|---|---|---|---|
| 67 | Fig 1 | 405x720 | 8-bit RGB JPEG | 72 | 17,865 | 2.04% |
| 68 | Fig 2 | 540x645 | 8-bit RGB JPEG | 72 | 22,561 | 2.16% |
| 69 | Fig 3 | 292x720 | 8-bit RGB JPEG | 72 | 14,250 | 2.26% |
| 70 | Fig 4 | 311x719 | 8-bit RGB JPEG | 72 | 20,391 | 3.04% |
| 71 | Fig 5 | 540x626 | 8-bit RGB JPEG | 72 | 39,465 | 3.89% |

Each figure page also carries a second, incidental raster: a 9x10 4-bit indexed PNG (~130 bytes) at
108 ppi, an Editorial Manager interface glyph rather than figure content.

**These are proof renderings, not the submitted files.** The deposited TIFFs are 300 ppi and 1905-2220
px wide. The percentage column is JPEG bytes over *that 72-ppi raster's* uncompressed size; measuring
it against the deposited TIFF instead gives 0.18-0.42% and means something different. The uploaded
originals reach a reviewer through the per-figure download links above, not through the proof image.

## What was uploaded

28 files, assembled at `5e77414`. md5s are of the copies actually sent.

| # | file | bytes | md5 |
|---|---|---|---|
| 1 | `CellWarp_PLOSONE_manuscript.docx` | 72,805 | `8d5a5237ee7bae13779b76d415f00c33` |
| 2 | `Fig1.tif` | 337,886 | `ee9d4bdb1fa2a990122dc0b859661b21` |
| 3 | `Fig2.tif` | 539,834 | `9450a8d52729ef9175e2925aeab83665` |
| 4 | `Fig3.tif` | 245,284 | `ae4289ace2d685f0e004062ea1476baf` |
| 5 | `Fig4.tif` | 189,214 | `d5646cf15755e91005d6d8646db33a92` |
| 6 | `Fig5.tif` | 630,936 | `be7c883a4b8ca73677072717bd82cd09` |
| 7 | `figS1_pipeline_validation.pdf` | 637,070 | `a3e97cafb7ec3e2f5e5f962354e55236` |
| 8 | `figS2_parameter_protocol_sensitivity.pdf` | 1,118,522 | `cedd035a56ddd10b41b924a6a1118376` |
| 9 | `figS3_bootstrap_rankings.pdf` | 525,295 | `8d92a430098918dc38021414e710ebf1` |
| 10 | `figS4_matched_scale_control.pdf` | 58,243 | `d200c824c5a769e098010f1feec8dee2` |
| 11 | `figS5_markernull.pdf` | 24,130 | `96204a01c93fbff5e64dae9acd2da403` |
| 12 | `table_S1.xlsx` | 9,161 | `48f3138a7c91390265cec72f5a290ccf` |
| 13 | `table_S2.xlsx` | 8,634 | `8870b738a2d93fea58d422d63267ecb0` |
| 14 | `table_S3.csv` | 1,759 | `85481d90eb703c13f16ca916e9c4cbbd` |
| 15 | `table_S4.csv` | 394 | `22d16410ce5bde2abd1a94ce2d4785cc` |
| 16 | `table_S5.csv` | 4,338 | `b961236aad7a1d4a939b50ee49021983` |
| 17 | `Table_S6_CPC1_driver_genes.xlsx` | 59,428 | `df2b3767768645c3a71cf123ff5db048` |
| 18 | `table_S7_layer1_housekeeping_exclusion.csv` | 190 | `302c45fb2aca9b41827b5a838b006c0b` |
| 19 | `table_S9_genestd_standardization.csv` | 298 | `64525cf0007cd47b0b01d851a9c76609` |
| 20 | `table_S9_schemeB_CPC1_markers.csv` | 175 | `87dcb0b27606aafc0ef4f97f1a6c3661` |
| 21 | `table_S10_markernull.csv` | 598 | `1ca1ef09b0eb49d0020e55bb4718be46` |
| 22 | `table_S11_gene_conservation.csv` | 2,788,747 | `18ccbdd76d18a25e2c86e7157f6789b0` |
| 23 | `table_S12_software_environment.csv` | 784 | `1a208d4e44f2854a2c104bd78d039b59` |
| 24 | `table_S13_test_inventory.xlsx` | 12,545 | `b88e57aebb7af47bfd4ea63e7cdbbe39` |
| 25 | `table_S14_layer2_dimension.csv` | 3,732 | `a7b2c5ba50293f6a6863f1bd1deda394` |
| 26 | `S1_Text.txt` | 40,043 | `412e83edd33b623dceb284802bdc68ef` |
| 27 | `S2_Text.txt` | 6,703 | `ad93400db8743c0d2b2214af7dd9f135` |
| 28 | `coverletter.txt` | 3,347 | `b495455ea5d8253f2cbd88261693bb15` |

**The manuscript DOCX md5 is not a pin.** A DOCX embeds a creation timestamp, so
`build_manuscript_docx.py` is content-reproducible but never byte-reproducible: rebuilding from the
same pinned text yields the same 72,805 bytes and a different md5. The value that *is* stable is the
source text it was built from, `manuscript_combined.txt` at md5
`35b9b55c552e32afd2d9dcff445074fa`, pinned in `reproduce/MANUSCRIPT_MD5`.

`coverletter.txt` is correspondence with the editor. It is ignored by `.gitignore:196`
(`cover*letter*.txt`) and is deliberately absent from both Zenodo archives; it is listed here because
it was submitted, not because it is archived.

## The basal-ganglia record's `References` relation — checked, correct, closed

`10.5281/zenodo.22073208`'s DataCite metadata carries `References -> 10.5281/zenodo.20735612`, which
is CODE **v1's version DOI**, not the concept DOI `20735611` the manuscript cites. That looks like a
mismatch and is not one.

The basal-ganglia deposit pins CellWarp to v1 deliberately. Its `requirements.txt` reads
`cellwarp @ git+https://github.com/sriramdevadas/cellwarp.git@pcompbiol-submission-2026-06`, and its
own comment states the pin is "the exact commit that produced the paper's results (tag
`pcompbiol-submission-2026-06` == commit `4fca942` == the Zenodo code deposit
`10.5281/zenodo.20735612`)". The tag peels to `4fca942`, and `4fca942` is what CODE v1 was built
from.

So a frozen replication points at the exact version it ran against, which is what it should do. The
relation is correct and must not be "corrected" to the concept DOI: the concept DOI resolves to v2,
which is **not** what the basal-ganglia results were produced with. No change to the Zenodo record.

## Ethics statement as entered

The Editorial Manager ethics field was filled with the manuscript's own ethics statement verbatim,
the text at `manuscript_combined.txt:136` under the Methods heading at `:134`. It is reproduced in
`docs/declarations.txt` Section 7 so that the repository and the submission form carry the same
bytes. The form answer to the human-participants question was "This study did not involve human
participants, data, specimens or tissue"; no PLOS Human Participants Research Checklist was required
or attached.
