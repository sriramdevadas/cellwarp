# SUBMISSION_RECORD

What was sent to PLOS ONE, as built and as compiled. Companion to `DEPOSIT_MANIFEST.md`, which
records the Zenodo side; this file records the journal side. Every value below was measured from
the artifact named beside it, not carried forward from a plan.

Source commit: `5e77414`. `git ls-files` at that commit: 1094 (this file is added after it).

## Identifiers

- **Manuscript number: `PONE-D-26-42571`**, assigned by PLOS ONE when the submission returned from
  technical check on 2026-08-25. This is the number to quote in all correspondence.
- **Original submission id: `PONE-S-26-56365`**, the identifier the 24 August submission carried
  before a manuscript number existed. Superseded, kept because the 24 August compiled PDF and the
  first technical-check correspondence are filed under it.
- **Provenance caution, and it applies to both.** Neither identifier is verifiable from the artifact
  it names. `PONE-S-26-56365` appears nowhere inside the 24 August compiled PDF -- not in page text,
  not in link annotations, not in the raw file bytes, where `PONE` and `56365` are both absent -- and
  its only local evidence was the downloaded filename. `PONE-D-26-42571` arrived in correspondence
  and has exactly the same status: Editorial-Manager-supplied, recorded on their authority, not
  derived from anything here. The 24 August cover sheet printed `Manuscript Number:` empty and was
  titled `--Manuscript Draft--`, which is what a pre-assignment submission looks like.

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

**This table is the 24 August submission, 28 files assembled at `5e77414`.** The 25 August
resubmission changed three of its rows and is recorded in the resubmission section below: the
manuscript DOCX was rebuilt (72,741 B, from the abstract edit at `4b2d14b`), `coverletter.txt` was
replaced by `coverletter.docx`, and the completed Human Participants Research Checklist was added,
taking the File Inventory from 28 items to 29. Every other row below is still current.

md5s are of the copies actually sent.

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

**The manuscript DOCX md5 is not a pin.** `build_manuscript_docx.py` is content-reproducible but
never byte-reproducible: rebuilding from the same pinned text yields the same byte count and a
different md5. The value that *is* stable is the source text it was built from,
`manuscript_combined.txt`, pinned in `reproduce/MANUSCRIPT_MD5`.

**The mechanism, corrected.** This section previously said a DOCX "embeds a creation timestamp".
That named the wrong cause. `dcterms:created` and `dcterms:modified` are **fixed** at
`2013-12-23T23:15:00Z` in both DOCX files produced here -- a constant inherited from python-docx's
shipped `default.docx`, not a wall clock. What actually varies is the **ZIP member timestamps** in
the local file headers, stamped at save. The conclusion was right and the stated reason was not.

**How that was found, because the method transfers.** Two consecutive builds returned the same md5
and the file was nearly recorded as byte-reproducible on that evidence. A third build, after a
deliberate two-second pause, differed. **A reproducibility test that runs faster than the resolution
of the clock it is testing will confirm itself.** Any such test here needs a pause longer than one
second between builds, or it is measuring nothing.

`coverletter.txt` is correspondence with the editor. It is ignored by `.gitignore:196`
(`cover*letter*.txt`) and is deliberately absent from both Zenodo archives; it is listed here because
it was submitted, not because it is archived.

## The resubmission — technical check, 2026-08-25

There have been **two submission events**, and everything above describes the first unless it says
otherwise.

**24 August 2026.** Submitted as `PONE-S-26-56365`. Compiled PDF 92 pages, 28 items, md5
`173564d321a57140960559d337d61d7c`, 5,238,533 bytes. That artifact is now superseded.

**25 August 2026.** Passed PLOS ONE's in-house technical check with **two administrative queries and
no scientific ones**, both satisfied on resubmission. Manuscript number `PONE-D-26-42571` assigned.

### What technical check raised, and what was done

1. **Cover letter format.** *"We notice that your Cover Letter in File Inventory are in Note pad.
   Please upload Cover Letter in File Inventory in word format."* Converted to `coverletter.docx`,
   one `.txt` line to one `.docx` paragraph, with the round trip asserted in the converter so the
   text is provably unchanged: 504 words, 21 paragraphs, `extracted == source`. The `.txt` remains
   the source of record; only the uploaded `.txt` was removed from the File Inventory.
2. **Human Participants Research Checklist.** Requested, and **completed rather than declined**. An
   earlier reading held that completing it would require asserting things untrue of this study; that
   was wrong, and it was wrong because the form had not been read. The form is conditional
   throughout, and every item has a designed exempt pathway: item 1 answered N/A with the
   manuscript's ethics statement pasted into its explanation box, items 2 and 3 answered N/A because
   their conditions -- prospective recruitment, and retrospective study of medical records or
   archived samples -- are untrue of a reanalysis of published atlases. Uploaded as file type
   `Other`.

### The abstract was edited in the same window

Not requested by PLOS. The window existed because of the two items above, and the abstract carried a
defect worth the cost of using it: the PLOS Computational Biology desk rejection had objected that
the paper did not define cell type geometry **or its relevance to the broader field**, and while the
construct was defined at word 2, the practice claim the paper contradicts appeared in the cover
letter and the introduction and **not in the abstract** -- the one document an Academic Editor reads.
The thesis also sat at word 229 of 300 against the cover letter's 11%. Seven of twelve sentences were
rewritten; five are byte-identical; 300 words to 280; no result, number or claim changed. Commit
`4b2d14b`.

### The rebuilt compiled PDF

```
pages        93            (was 92; +1 for the checklist placeholder)
bytes        5,350,072
inventory    29 items      Cover Letter + Manuscript + 5 Figures + 21 SI + 1 Other
```

**Its md5 is not a stable identifier and must not be recorded as one.** Three downloads of the same
document taken on 2026-08-25 gave **5,350,072 bytes every time and two different md5s**:
`1358b39d9a8c88ef7761db5cb27e9b76` for a build stamped `D:20260825120731`, and
`766196b3dd07e3d3a6d467af01a4ebf1` for one stamped `D:20260825124135`, downloaded twice. Page text is
**identical across all three**. Editorial Manager recompiles the PDF and stamps the build time into
it, so its md5 identifies a download rather than the document. **Byte count and page count are the
stable identifiers.**

### Verified on the rebuilt build

Two independent checks agreed on every value:

- abstract form field **280 words** by `str.split()`, word-identical to the manuscript's abstract
- abstract superscripts are true `U+207B` + `U+2076` / `U+2074`, ArialUnicodeMS, rendered raised and
  reduced -- body runs 8.0 pt against 12.0 pt, raised 4.51-4.61 pt
- figures on pp. 67-71, `Fig1.tif` … `Fig5.tif` ascending, Fig 2 the full A/B/C composite
- DOIs: `20735611` pp. 4, 7, 46 · `20735639` pp. 5, 7, 46 · `22073208` pp. 5, 47
- `mint|minted|minting|pending`: zero, whole-word and as substrings
- `_S8`: zero occurrences anywhere in the file

**Superscript runs: 51 total carrying 96 characters**, of which **44 are minus-led**. Both figures
circulated and both are right; they count different things. Decomposed from the DOCX rather than
from either report: `−4` x33, `−6` x8, `−5` x2, `−14` x1 = 44 minus-led, plus `ᵀ` x5, the Frobenius
`²` x1, and the Fig 5 caption's `⁶` x1 = 51.

### What technical check did NOT query

Absence of a query is evidence, and nobody will reconstruct it later. All of the following passed
**unchanged and unremarked**:

- **the Data Availability statement**, including its volunteered disclosure that three classes of
  intermediate input are not reconstructible from the deposit -- the rolling-path UCSC/DoRothEA
  sources, the ChEMBL and DrugBank derived tables with no recorded release, and the deposited
  artifacts with no producing script
- **the `.txt` Supporting Information files**, `S1_Text.txt` and `S2_Text.txt`, uploaded as plain
  text rather than converted -- the DOC/DOCX/RTF rule reaches the manuscript, not SI
- **all five figures**, at 300 ppi RGB LZW, including the sub-8-pt text recorded as
  `AFTER_ACCEPTANCE.md` item 10
- **the ethics statement's content**, queried only as a missing checklist and never as to what it said
- **the references**

The cover letter's *format* was queried and its *text* was not.

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

The Editorial Manager ethics field is built on the manuscript's own ethics statement, the text at
`manuscript_combined.txt:136` under the Methods heading at `:134`.

**They are not the same bytes, and an earlier version of this section said they were.** Measured:

```
manuscript_combined.txt:136              495 bytes   65 words   3 sentences
docs/declarations.txt line 205           495 bytes   byte-identical to the manuscript
checklist explanation box                495 bytes   byte-identical to the manuscript
Editorial Manager ethics field           626 bytes   83 words   4 sentences
```

**Three places hold 495 and one holds 626.** The difference is purely additive:
`EM.startswith(manuscript)` is True, the three shared sentences are byte-identical at 184, 189 and
120 characters, and the remainder is one appended sentence of 131 characters including its leading
space -- "No human participants were recruited, no specimens or tissue were collected, and no
identifiable private information was accessed." That sentence **appears nowhere in this repository**,
searched on six distinct phrases from it.

So it is a scope difference rather than a substantive one, but the claim of identical bytes was
false as written. The Editorial Manager ethics screen instructs in bold that all information entered
there be included in the Methods section, and as submitted that instruction is not satisfied for the
fourth sentence. `AFTER_ACCEPTANCE.md` item 13 carries the full analysis and both resolution costs.
**Neither resolution is adopted here**; it is the coordinator's decision and carries no deadline.

The form answer to the human-participants question was "This study did not involve human
participants, data, specimens or tissue". At the 24 August submission no Human Participants Research
Checklist was attached; one was completed and uploaded at the 25 August resubmission -- see the
resubmission section below.
