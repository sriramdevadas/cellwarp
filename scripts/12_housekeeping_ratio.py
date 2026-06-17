#!/usr/bin/env python3
"""
Thread 4A: Housekeeping gene ratio vs Procrustes rigidity.

Biology:
    Housekeeping (HK) genes are constitutively expressed across tissues and cell
    types. They encode core cellular machinery (ribosomal proteins, metabolic
    enzymes, chaperones, proteasome subunits). Because HK genes serve functions
    in many cell types simultaneously, mutations altering their expression have
    pleiotropic effects — they hurt many tissues at once. This pleiotropic
    constraint means evolution cannot easily change HK gene expression without
    deleterious consequences.

    Lineage-specific genes, by contrast, are active in only one or a few cell
    types. Mutations in these genes have more limited downstream effects, so they
    face weaker purifying selection and can drift more freely between species.

    Hypothesis: cell types that derive more of their expression "budget" from
    housekeeping genes should be more evolutionarily rigid (lower Procrustes
    residual magnitude), because a larger share of their program is locked in
    by pleiotropy.

Math:
    For each cell type c with mean expression vector μ_c ∈ R^G:
        hk_ratio(c) = mean(μ_c[HK]) / mean(μ_c[all])
    where HK ⊂ {1,...,G} is the set of housekeeping gene indices.

    We then compute Spearman ρ between hk_ratio and residual magnitude across
    all 35 cell types.

Data:
    Housekeeping genes: HOUNKPE_HOUSEKEEPING_GENES from MSigDB (1,129 genes).
    Originally from Hounkpe et al., curated set of genes expressed uniformly
    across human tissues. We use this because the Eisenberg & Levanon 2013 list
    (3,804 genes) was unavailable for download (TAU server 403).

    Centroids: output/phase2/scaled_35types/centroids_{human,mouse}_35.csv
    Residuals: output/phase2/scaled_35types/procrustes_results_35.json
    Annotations: output/phase2/developmental_constraint/developmental_annotations.csv
    Ortholog map: data/phase1/orthologs_human_mouse.csv

Output:
    output/phase2/mechanistic/housekeeping/
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "phase2" / "mechanistic" / "housekeeping"
OUT.mkdir(parents=True, exist_ok=True)

CENTROIDS_H = ROOT / "output" / "phase2" / "scaled_35types" / "centroids_human_35.csv"
CENTROIDS_M = ROOT / "output" / "phase2" / "scaled_35types" / "centroids_mouse_35.csv"
RESULTS_JSON = ROOT / "output" / "phase2" / "scaled_35types" / "procrustes_results_35.json"
ANNOTATIONS = ROOT / "output" / "phase2" / "developmental_constraint" / "developmental_annotations.csv"
ORTHOLOGS = ROOT / "data" / "phase1" / "orthologs_human_mouse.csv"

# ── HOUNKPE housekeeping gene set (1,129 genes from MSigDB) ──────────────────
# Source: HOUNKPE_HOUSEKEEPING_GENES, MSigDB C2:CGP
# https://www.gsea-msigdb.org/gsea/msigdb/cards/HOUNKPE_HOUSEKEEPING_GENES
HOUNKPE_HK_GENES = {
    "AAMP", "AAR2", "AARS1", "AARSD1", "AASDHPPT", "ABCE1", "ABCF2", "ACAD9",
    "ACBD3", "ACBD6", "ACTR1A", "ACTR1B", "ACTR2", "ADH5", "ADIPOR2", "ADRM1",
    "AGAP3", "AGGF1", "AHSA1", "AIDA", "AIP", "AK2", "AK3", "AKIRIN1",
    "AKIRIN2", "AKT1", "ALDH9A1", "ALKBH5", "AMFR", "ANAPC2", "ANKFY1",
    "ANKLE2", "ANKRD10", "ANKRD40", "ANKRD52", "ANKS1A", "ANP32A", "ANP32B",
    "AP1G1", "AP2M1", "AP3M1", "AP3S2", "APEH", "APH1A", "API5", "APPL1",
    "ARAF", "ARCN1", "ARF1", "ARF4", "ARFGAP2", "ARFGAP3", "ARFGEF2",
    "ARHGAP1", "ARIH2", "ARL1", "ARL14EP", "ARL6IP1", "ARL8B", "ARMC1",
    "ARMCX3", "ARNT", "ARPC2", "ARPP19", "ASNSD1", "ASXL1", "ATF6", "ATF6B",
    "ATG101", "ATG3", "ATP5PF", "ATP6AP2", "ATP6V0E1", "ATPAF1", "B4GALT7",
    "BAD", "BAG1", "BAG6", "BANF1", "BAP1", "BAZ1B", "BCAP31", "BCCIP",
    "BCKDK", "BET1L", "BFAR", "BLOC1S2", "BLTP2", "BMI1", "BNIP3L", "BOD1",
    "BORCS7", "BRD3", "BRD4", "BRK1", "BRMS1", "BSDC1", "BTBD1", "BTBD2",
    "BTF3", "BTF3L4", "BUB3", "BUD31", "BZW1", "C11orf98", "C14orf119",
    "C1orf174", "C1orf43", "C6orf120", "C6orf89", "C9orf78", "CAB39",
    "CACTIN", "CALM1", "CALR", "CAND1", "CASC3", "CBX1", "CBY1", "CCAR2",
    "CCDC124", "CCDC22", "CCDC47", "CCDC50", "CCDC97", "CCNI", "CCNK", "CCT2",
    "CCT4", "CCZ1", "CD2BP2", "CDC123", "CDC16", "CDC23", "CDC26", "CDC34",
    "CDC37", "CDC37L1", "CDC42", "CDC42SE1", "CDC5L", "CDK4", "CDKN1B",
    "CEBPZ", "CEBPZOS", "CFDP1", "CFL1", "CGGBP1", "CHCHD1", "CHMP2A",
    "CHMP2B", "CHMP3", "CHMP4B", "CHMP6", "CHP1", "CHTOP", "CIPC", "CLNS1A",
    "CLPP", "CLPTM1", "CLPTM1L", "CMAS", "CMPK1", "CMTR1", "CNBP", "CNOT11",
    "CNPPD1", "COA3", "COA5", "COA6", "COASY", "COG3", "COG7", "COMMD7",
    "COMT", "COPA", "COPB2", "COPE", "COPG1", "COPS4", "COPS5", "COPS6",
    "COPZ1", "COQ5", "COX15", "COX17", "CREBBP", "CRK", "CRKL", "CRTC3",
    "CSNK1A1", "CSNK1D", "CSNK1G2", "CSNK2A1", "CSNK2A2", "CSNK2B", "CTBP1",
    "CTDNEP1", "CTNNB1", "CTNNBL1", "CUEDC2", "CUL1", "CUL3", "CUL4B",
    "CUTA", "CYB5B", "DAD1", "DAP3", "DARS1", "DCAF12", "DCAF15", "DCAF5",
    "DCTN2", "DCTN3", "DCTN4", "DCTPP1", "DDB1", "DDX1", "DDX17", "DDX18",
    "DDX19A", "DDX23", "DDX27", "DDX47", "DDX51", "DEAF1", "DENND6A", "DENR",
    "DERL1", "DHPS", "DHX15", "DHX29", "DHX38", "DHX40", "DIMT1", "DNAJA1",
    "DNAJA2", "DNAJA3", "DNAJB11", "DNAJB6", "DNAJC1", "DNAJC11", "DNAJC13",
    "DNAJC3", "DNAJC4", "DNAJC5", "DNAJC7", "DNAJC8", "DNAJC9", "DNLZ",
    "DNTTIP1", "DPY30", "DRG1", "DUSP11", "DYNC1LI1", "DYNLRB1", "E2F4",
    "E4F1", "EAPP", "EBAG9", "EBNA1BP2", "EDC3", "EEF1D", "EEFSEC", "EGLN2",
    "EIF1AX", "EIF1B", "EIF2AK1", "EIF2B1", "EIF2D", "EIF2S1", "EIF2S2",
    "EIF3A", "EIF3B", "EIF3D", "EIF3E", "EIF3G", "EIF3I", "EIF3K", "EIF3M",
    "EIF4A1", "EIF4A3", "EIF4B", "EIF4G2", "EIF4H", "EIF5A", "EIF5B",
    "ELAVL1", "ELK1", "ELL", "ELOF1", "EMC2", "EMC3", "EMC4", "EMC6", "EMC7",
    "EMC8", "EMD", "ENSA", "EPN1", "ERAL1", "ERGIC1", "ERGIC3", "ERH", "ERI3",
    "ERLEC1", "ERP29", "ERP44", "ETF1", "EXOC2", "EXOSC4", "FAF2", "FAM168A",
    "FAM168B", "FAM20B", "FAM32A", "FAM50A", "FAM98A", "FBL", "FBXL3",
    "FBXO28", "FBXO42", "FBXW4", "FBXW5", "FCF1", "FEM1B", "FHIP2A", "FIBP",
    "FIS1", "FIZ1", "FKBP8", "FLCN", "FRG1", "FTSJ3", "FUBP3", "FXR2",
    "FZR1", "G3BP1", "G3BP2", "G6PC3", "GABARAPL2", "GATAD1", "GATAD2B",
    "GBF1", "GDI2", "GET3", "GET4", "GFUS", "GGA2", "GID8", "GLE1", "GNB1",
    "GNL2", "GNPAT", "GNPTG", "GOLGA3", "GOLGA4", "GOLPH3", "GORASP2",
    "GOSR2", "GPBP1", "GPI", "GPKOW", "GPN2", "GPN3", "GPX4", "GSK3A",
    "GSPT1", "GSTO1", "GTF2B", "GTF2E2", "GTF2F1", "GTF2F2", "GTF3A",
    "GTF3C3", "GTF3C6", "GTPBP1", "GUK1", "H2AZ2", "H3-3B", "HAPSTR1",
    "HAT1", "HAX1", "HCFC1", "HDLBP", "HERPUD2", "HGS", "HIGD1A", "HINT1",
    "HIPK1", "HMGXB3", "HNRNPA1", "HNRNPA2B1", "HNRNPAB", "HNRNPD",
    "HNRNPDL", "HNRNPM", "HNRNPU", "HNRNPUL1", "HNRNPUL2", "HP1BP3", "HSBP1",
    "HSD17B12", "HSF1", "HSP90AB1", "HSP90B1", "HSPA4", "HSPE1", "HTATSF1",
    "HTT", "ICMT", "IER3IP1", "IGBP1", "IK", "ILF2", "ILF3", "ILKAP",
    "INO80", "INTS10", "INTS12", "IPO5", "IPO8", "IPO9", "IQSEC1", "ISCA2",
    "ISCU", "ITCH", "IWS1", "JAGN1", "JAK1", "JOSD1", "JTB", "KARS1",
    "KAT8", "KBTBD2", "KCMF1", "KCTD20", "KCTD5", "KDELR1", "KDELR2",
    "KDM1A", "KDM2A", "KDM4B", "KEAP1", "KHDRBS1", "KHSRP", "KIAA0319L",
    "KIFBP", "KLHDC10", "KPNA1", "KPNA3", "KPNA4", "KPNA6", "KRAS", "KXD1",
    "LAMP2", "LAMTOR1", "LAMTOR5", "LARP1", "LARP4B", "LARS1", "LCMT1",
    "LDB1", "LEMD2", "LENG1", "LENG8", "LMAN1", "LMAN2", "LMF2", "LONP1",
    "LONP2", "LRPAP1", "LRPPRC", "LRRC41", "LRRC59", "LSM1", "LSM12",
    "LSM14A", "LSM14B", "LSM4", "LTV1", "LUC7L3", "LYSET", "MAEA", "MAGOH",
    "MAN1A2", "MAP2K1", "MAP2K2", "MAPK1", "MAPK1IP1L", "MAPK8IP3", "MAPRE1",
    "MARCHF5", "MAU2", "MAX", "MAZ", "MBD1", "MBD2", "MCCC2", "MCMBP",
    "MCRS1", "MCTS1", "MEA1", "MECP2", "MED10", "MED13", "MED19", "MED21",
    "MED25", "MED28", "MED4", "MED9", "MEPCE", "METTL9", "MFN2", "MFSD14A",
    "MGRN1", "MIEF1", "MKRN2", "MLEC", "MLF2", "MMADHC", "MORC2", "MORF4L1",
    "MORF4L2", "MPHOSPH10", "MRPL14", "MRPL17", "MRPL18", "MRPL22", "MRPL24",
    "MRPL28", "MRPL36", "MRPL37", "MRPL38", "MRPL39", "MRPL41", "MRPL43",
    "MRPL44", "MRPL49", "MRPL51", "MRPL57", "MRPL9", "MRPS10", "MRPS14",
    "MRPS16", "MRPS17", "MRPS2", "MRPS25", "MRPS28", "MRPS30", "MRPS5",
    "MRPS9", "MRTO4", "MSL1", "MSL2", "MTCH1", "MTDH", "MTFR1L", "MTOR",
    "MTPN", "MYCBP2", "MYDGF", "NACA", "NAP1L1", "NAP1L4", "NAPA", "NAPG",
    "NARS1", "NAT10", "NCBP1", "NCBP2", "NCL", "NCOA5", "NDST1", "NDUFAF2",
    "NDUFB8", "NEDD8", "NELFB", "NFATC2IP", "NFX1", "NFYB", "NGDN", "NGRN",
    "NIFK", "NMD3", "NMT1", "NOB1", "NOL7", "NOP14", "NOP9", "NOSIP",
    "NPEPPS", "NPLOC4", "NRBP1", "NSA2", "NSMCE1", "NSMCE4A", "NT5C2",
    "NUCKS1", "NUDC", "NUDT16L1", "NUDT21", "NUDT9", "NUFIP2", "NUMA1",
    "NUP133", "NUP153", "NUS1", "NXF1", "OGT", "OS9", "OSBP", "OST4", "OSTC",
    "OTUB1", "OTUD4", "OXSR1", "P4HB", "PA2G4", "PABIR1", "PABPN1", "PACS1",
    "PACSIN2", "PAF1", "PAFAH1B1", "PAFAH1B2", "PAIP1", "PARK7", "PARN",
    "PARP1", "PATL1", "PBDC1", "PBX2", "PCBP1", "PCBP2", "PCGF3", "PCIF1",
    "PCMT1", "PCMTD2", "PCSK7", "PCYOX1", "PCYT1A", "PDAP1", "PDCD6",
    "PDCD7", "PDIA3", "PDIA6", "PDPK1", "PDS5A", "PDZD11", "PELP1", "PFDN1",
    "PFDN2", "PGAM5", "PGRMC2", "PHAF1", "PHF12", "PHF5A", "PHRF1", "PIGH",
    "PIN1", "PIP5K1C", "PITHD1", "PLAA", "PMPCA", "PMPCB", "PNO1", "PNRC2",
    "POFUT1", "POLD2", "POLE3", "POLR1D", "POLR2C", "POLR2J", "POLR2K",
    "POLR3C", "POP7", "PPIB", "PPID", "PPP1CC", "PPP1R11", "PPP1R15B",
    "PPP1R37", "PPP1R8", "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R5A", "PPP2R5C",
    "PPP2R5E", "PPP5C", "PPP6C", "PPP6R1", "PPP6R2", "PPT1", "PPTC7", "PRCC",
    "PRDM4", "PRDX5", "PREB", "PREP", "PRKAA1", "PRKAR1A", "PRKAR2A",
    "PRKCSH", "PRKRA", "PRPF18", "PRPF19", "PRPF38A", "PRPS1", "PRPSAP1",
    "PSMA2", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB4", "PSMB6",
    "PSMB7", "PSMC1", "PSMC3", "PSMC5", "PSMD1", "PSMD10", "PSMD11",
    "PSMD12", "PSMD14", "PSMD3", "PSMD4", "PSMD5", "PSMD6", "PSMD7", "PSMD8",
    "PSME3", "PSME3IP1", "PSMG2", "PSMG3", "PTEN", "PTPN11", "PTRH2",
    "PUF60", "PUM2", "PURB", "PWP1", "PWP2", "QRICH1", "R3HCC1", "RAB11A",
    "RAB11B", "RAB14", "RAB1B", "RAB2A", "RAB35", "RAB3GAP1", "RAB5A",
    "RAB6A", "RAB7A", "RAB9A", "RABGGTB", "RABL6", "RAC1", "RAD21", "RAD23A",
    "RAD23B", "RALA", "RALBP1", "RALY", "RANBP1", "RANBP10", "RANBP2",
    "RANBP9", "RANGAP1", "RARS2", "RBCK1", "RBM14", "RBM17", "RBM25", "RBM4",
    "RBM42", "RBM8A", "RBMX2", "RBX1", "RC3H2", "RER1", "RERE", "RGP1",
    "RHBDD2", "RHEB", "RHOT2", "RING1", "RMND5A", "RNF113A", "RNF139",
    "RNF167", "RNF168", "RNF181", "RNF20", "RNF216", "RNPS1", "RPA1", "RPA2",
    "RPL15", "RPL22", "RPL35A", "RPL4", "RPL7L1", "RPN1", "RPN2", "RPS19BP1",
    "RRAGA", "RRN3", "RRP36", "RRP7A", "RRP9", "RSL24D1", "RSRC2", "RTCB",
    "RTF1", "RUVBL1", "SAP18", "SAR1A", "SARNP", "SARS1", "SART1", "SART3",
    "SAV1", "SBDS", "SCAF1", "SCYL1", "SCYL2", "SDF4", "SDHAF2", "SEC13",
    "SEC24B", "SEC31A", "SEC61B", "SEC62", "SENP2", "SENP5", "SEPHS1",
    "SERBP1", "SERP1", "SETD5", "SF1", "SF3A2", "SF3A3", "SF3B1", "SF3B2",
    "SF3B4", "SF3B5", "SF3B6", "SFPQ", "SGTA", "SH2B1", "SH3GL1", "SHARPIN",
    "SKIC3", "SKP1", "SLC25A11", "SLC25A28", "SLC25A3", "SLC25A5", "SLC30A5",
    "SLC30A9", "SLC35A1", "SLC35A4", "SLC35E1", "SLC39A6", "SLC4A1AP",
    "SLIRP", "SLTM", "SLU7", "SMAD2", "SMAP1", "SMARCA5", "SMARCB1",
    "SMARCC1", "SMC1A", "SMDT1", "SMIM11", "SMIM12", "SMPD1", "SMU1",
    "SMYD5", "SNAP29", "SNAPIN", "SND1", "SNF8", "SNRNP35", "SNRNP70",
    "SNRPA", "SNRPB", "SNRPD3", "SNX1", "SNX12", "SNX19", "SNX27", "SNX3",
    "SOD1", "SON", "SP1", "SP2", "SPAG7", "SPAG9", "SPG11", "SPG7", "SPPL3",
    "SPRYD3", "SPTSSA", "SREBF2", "SRM", "SRP14", "SRP9", "SRPRB", "SRRM1",
    "SRSF1", "SRSF10", "SRSF11", "SRSF2", "SRSF3", "SRSF4", "SRSF6",
    "SRSF9", "SSBP1", "SSNA1", "SSR2", "SSR3", "SSRP1", "SSU72", "ST13",
    "STAT3", "STAU1", "STIM1", "STIP1", "STK11", "STK24", "STK38", "STRAP",
    "STRN4", "STT3B", "STUB1", "STX4", "STX5", "STX8", "SUB1", "SUGT1",
    "SUMO1", "SUPT6H", "SURF4", "SYF2", "SYNCRIP", "TAB1", "TADA2B", "TADA3",
    "TAF10", "TAOK2", "TARDBP", "TAX1BP1", "TBC1D20", "TBCA", "TBCB", "TBCD",
    "TBK1", "TBL1XR1", "TCF25", "TCP1", "TERF2IP", "TEX261", "THAP11",
    "THOC7", "THRAP3", "TIMM17A", "TIMM17B", "TIMM22", "TM9SF3", "TMBIM6",
    "TMCO1", "TMED1", "TMED10", "TMED4", "TMED9", "TMEM101", "TMEM127",
    "TMEM147", "TMEM165", "TMEM183A", "TMEM203", "TMEM219", "TMEM230",
    "TMEM245", "TMEM248", "TMEM256", "TMEM258", "TMEM259", "TMEM30A",
    "TMEM42", "TMEM50A", "TMEM60", "TMX3", "TMX4", "TNKS", "TNKS2", "TOLLIP",
    "TOMM20", "TOMM22", "TOMM70", "TOR1A", "TOR1B", "TPI1", "TPP2", "TRAM1",
    "TRAPPC3", "TRIM28", "TRIP12", "TRIP4", "TRMT112", "TRPC4AP", "TSN",
    "TSPYL1", "TSR1", "TSR3", "TTC1", "TUBG1", "TUSC2", "TUT1", "TXN2",
    "TXNDC12", "TXNDC9", "TXNL1", "U2AF1", "U2AF2", "UBA1", "UBA2", "UBAC1",
    "UBAP2L", "UBE2D2", "UBE2F", "UBE2G1", "UBE2I", "UBE2J1", "UBE2K",
    "UBE2Q1", "UBE2R2", "UBE2Z", "UBE3C", "UBFD1", "UBL4A", "UBL5", "UBL7",
    "UBP1", "UBQLN1", "UBQLN4", "UBR4", "UBXN1", "UFC1", "UFL1", "UFM1",
    "UGP2", "UROD", "USB1", "USP14", "USP22", "USP39", "USP4", "USP48",
    "UTP18", "UTP3", "VCP", "VDAC1", "VEZF1", "VPS16", "VPS25", "VPS26A",
    "VPS26B", "VPS29", "VPS36", "VPS37C", "VPS51", "VPS52", "VTI1B", "WAC",
    "WASL", "WBP4", "WDR45B", "WDR46", "WDR55", "WDR82", "WDR83OS", "WHAMM",
    "WIPI2", "WRNIP1", "XPOT", "YARS1", "YIF1A", "YIPF3", "YKT6", "YME1L1",
    "YTHDC1", "YTHDF1", "YWHAB", "YWHAE", "YWHAQ", "YY1", "ZBTB17",
    "ZC3H11A", "ZC3H15", "ZC3H18", "ZC3H3", "ZDHHC6", "ZDHHC7", "ZFP91",
    "ZFPL1", "ZFYVE21", "ZMAT2", "ZMIZ1", "ZMPSTE24", "ZMYND19", "ZNF146",
    "ZNF330", "ZNF593", "ZNF706", "ZNF777", "ZNF865", "ZRANB1", "ZRANB2",
    "ZRSR2", "ZYG11B",
}


def load_ortholog_map() -> dict[str, str]:
    """Build symbol→Ensembl mapping from cached ortholog table."""
    orth = pd.read_csv(ORTHOLOGS)
    # human_gene_name → human_ensembl_id
    return dict(zip(orth["human_gene_name"], orth["human_ensembl_id"]))


def map_hk_to_ensembl(symbol_to_ens: dict[str, str]) -> set[str]:
    """Map housekeeping gene symbols to Ensembl IDs in our gene space."""
    mapped = set()
    unmapped = []
    for sym in HOUNKPE_HK_GENES:
        ens = symbol_to_ens.get(sym)
        if ens:
            mapped.add(ens)
        else:
            unmapped.append(sym)
    return mapped, unmapped


def compute_hk_ratio(centroids: pd.DataFrame, hk_ensembl: set[str]) -> pd.Series:
    """
    Compute housekeeping expression ratio for each cell type.

    For each cell type c:
        hk_ratio(c) = mean(expr[HK genes]) / mean(expr[all genes])

    A ratio >1 means HK genes are on average more highly expressed than the
    genome-wide mean. A ratio <1 means HK genes are less expressed.
    """
    gene_cols = [c for c in centroids.columns if c.startswith("ENSG")]
    hk_cols = [c for c in gene_cols if c in hk_ensembl]

    expr_all = centroids[gene_cols].values  # (n_types, n_genes)
    expr_hk = centroids[hk_cols].values     # (n_types, n_hk)

    mean_all = expr_all.mean(axis=1)
    mean_hk = expr_hk.mean(axis=1)

    ratio = mean_hk / mean_all
    return pd.Series(ratio, index=centroids["cell_type"], name="hk_ratio")


def load_residuals() -> pd.Series:
    """Load per-cell-type Procrustes residual magnitudes."""
    with open(RESULTS_JSON) as f:
        data = json.load(f)
    residuals = {ct: info["magnitude"] for ct, info in data["residuals"].items()}
    return pd.Series(residuals, name="residual_magnitude")


def load_annotations() -> pd.DataFrame:
    """Load developmental annotations for cell category coloring."""
    return pd.read_csv(ANNOTATIONS)


def assign_category(row: pd.Series) -> str:
    """Assign broad cell category from lineage annotations."""
    lin = row.get("lineage", "")
    if lin == "hematopoietic":
        return "immune"
    elif lin == "epithelial":
        return "epithelial"
    elif lin in ("mesenchymal", "endothelial"):
        return "stromal"
    else:
        return "other"


def plot_scatter(
    hk_ratio: pd.Series,
    residuals: pd.Series,
    categories: pd.Series,
    species: str,
    rho: float,
    pval: float,
    out_path: Path,
):
    """Scatter plot: hk_ratio vs residual magnitude, colored by category."""
    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {
        "immune": "#2196F3",
        "epithelial": "#4CAF50",
        "stromal": "#FF9800",
        "other": "#9E9E9E",
    }

    # Align all series
    common = hk_ratio.index.intersection(residuals.index)
    hk = hk_ratio.loc[common]
    res = residuals.loc[common]
    cats = categories.loc[common]

    for cat, color in colors.items():
        mask = cats == cat
        if mask.sum() == 0:
            continue
        ax.scatter(
            hk.loc[mask], res.loc[mask],
            c=color, label=cat, s=60, alpha=0.8, edgecolors="white", linewidth=0.5,
        )

    # Label points
    for ct in common:
        label = ct
        # Abbreviate long names
        if len(label) > 25:
            label = label[:22] + "..."
        ax.annotate(
            label, (hk.loc[ct], res.loc[ct]),
            fontsize=5.5, alpha=0.75,
            xytext=(4, 4), textcoords="offset points",
        )

    ax.set_xlabel("Housekeeping Gene Expression Ratio\n(mean HK expr / mean all expr)", fontsize=11)
    ax.set_ylabel("Procrustes Residual Magnitude\n(evolutionary divergence)", fontsize=11)
    ax.set_title(
        f"Housekeeping Ratio vs Evolutionary Rigidity — {species}\n"
        f"Spearman ρ = {rho:.3f}, p = {pval:.4f}",
        fontsize=12,
    )
    ax.legend(title="Cell category", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    print("=" * 70)
    print("Thread 4A: Housekeeping Gene Ratio vs Procrustes Rigidity")
    print("=" * 70)

    # ── Step 1: Map HK genes to our Ensembl gene space ─────────────────
    print("\n── Step 1: Housekeeping gene mapping ──")
    sym_to_ens = load_ortholog_map()
    hk_ensembl, unmapped = map_hk_to_ensembl(sym_to_ens)
    print(f"  HOUNKPE housekeeping gene set: {len(HOUNKPE_HK_GENES)} symbols")
    print(f"  Mapped to Ensembl IDs in our gene space: {len(hk_ensembl)}")
    print(f"  Unmapped (not in 1:1 ortholog set): {len(unmapped)}")
    if unmapped:
        print(f"  First 10 unmapped: {unmapped[:10]}")

    # ── Step 2: Load centroids ─────────────────────────────────────────
    print("\n── Step 2: Loading centroids ──")
    centroids_h = pd.read_csv(CENTROIDS_H)
    centroids_m = pd.read_csv(CENTROIDS_M)
    print(f"  Human: {centroids_h.shape[0]} cell types × {centroids_h.shape[1]-1} genes")
    print(f"  Mouse: {centroids_m.shape[0]} cell types × {centroids_m.shape[1]-1} genes")

    # Verify HK genes are in centroids
    gene_cols = [c for c in centroids_h.columns if c.startswith("ENSG")]
    hk_in_centroids = hk_ensembl.intersection(gene_cols)
    print(f"  HK genes present in centroid gene space: {len(hk_in_centroids)} / {len(hk_ensembl)}")

    # ── Step 3: Compute HK ratio ───────────────────────────────────────
    print("\n── Step 3: Housekeeping ratio per cell type ──")
    hk_ratio_h = compute_hk_ratio(centroids_h, hk_ensembl)
    hk_ratio_m = compute_hk_ratio(centroids_m, hk_ensembl)

    print(f"\n  Human HK ratios (top 5 highest):")
    for ct, r in hk_ratio_h.sort_values(ascending=False).head(5).items():
        print(f"    {ct}: {r:.4f}")
    print(f"  Human HK ratios (top 5 lowest):")
    for ct, r in hk_ratio_h.sort_values().head(5).items():
        print(f"    {ct}: {r:.4f}")

    # ── Step 4: Load Procrustes residuals ──────────────────────────────
    print("\n── Step 4: Loading Procrustes residual magnitudes ──")
    residuals = load_residuals()
    print(f"  Loaded residuals for {len(residuals)} cell types")

    # ── Step 5: Correlate and visualize (Human) ────────────────────────
    print("\n── Step 5: Human — HK ratio vs residual magnitude ──")
    common = hk_ratio_h.index.intersection(residuals.index)
    rho_h, pval_h = stats.spearmanr(hk_ratio_h.loc[common], residuals.loc[common])
    print(f"  N = {len(common)} cell types")
    print(f"  Spearman ρ = {rho_h:.4f}")
    print(f"  p-value = {pval_h:.4f}")

    if pval_h < 0.05:
        direction = "NEGATIVE (rigid types have higher HK ratio)" if rho_h < 0 else "POSITIVE (diverged types have higher HK ratio)"
        print(f"  SIGNIFICANT at α=0.05. Direction: {direction}")
    else:
        print(f"  NOT significant at α=0.05.")

    # Load annotations for coloring
    annot = load_annotations()
    annot["category"] = annot.apply(assign_category, axis=1)
    cat_map = annot.set_index("cell_type")["category"]

    plot_scatter(
        hk_ratio_h, residuals, cat_map,
        "Human", rho_h, pval_h,
        OUT / "hk_ratio_vs_residual_human.png",
    )

    # ── Step 6: Mouse replication ──────────────────────────────────────
    print("\n── Step 6: Mouse — HK ratio vs residual magnitude ──")
    common_m = hk_ratio_m.index.intersection(residuals.index)
    rho_m, pval_m = stats.spearmanr(hk_ratio_m.loc[common_m], residuals.loc[common_m])
    print(f"  N = {len(common_m)} cell types")
    print(f"  Spearman ρ = {rho_m:.4f}")
    print(f"  p-value = {pval_m:.4f}")

    if pval_m < 0.05:
        direction = "NEGATIVE (rigid types have higher HK ratio)" if rho_m < 0 else "POSITIVE (diverged types have higher HK ratio)"
        print(f"  SIGNIFICANT at α=0.05. Direction: {direction}")
    else:
        print(f"  NOT significant at α=0.05.")

    plot_scatter(
        hk_ratio_m, residuals, cat_map,
        "Mouse", rho_m, pval_m,
        OUT / "hk_ratio_vs_residual_mouse.png",
    )

    # ── Cross-species HK ratio agreement ───────────────────────────────
    print("\n── Cross-species HK ratio agreement ──")
    common_both = hk_ratio_h.index.intersection(hk_ratio_m.index)
    rho_cross, pval_cross = stats.spearmanr(
        hk_ratio_h.loc[common_both], hk_ratio_m.loc[common_both]
    )
    print(f"  N = {len(common_both)} cell types")
    print(f"  Spearman ρ (human HK ratio vs mouse HK ratio) = {rho_cross:.4f}")
    print(f"  p-value = {pval_cross:.6f}")
    if rho_cross > 0.7:
        print("  STRONG agreement: HK ratios are conserved across species.")
    elif rho_cross > 0.4:
        print("  MODERATE agreement: HK ratios are partially conserved.")
    else:
        print("  WEAK agreement: HK ratios differ between species.")

    # Cross-species scatter
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        hk_ratio_h.loc[common_both], hk_ratio_m.loc[common_both],
        c="#607D8B", s=60, alpha=0.8, edgecolors="white", linewidth=0.5,
    )
    for ct in common_both:
        label = ct if len(ct) <= 25 else ct[:22] + "..."
        ax.annotate(
            label, (hk_ratio_h.loc[ct], hk_ratio_m.loc[ct]),
            fontsize=5.5, alpha=0.75, xytext=(4, 4), textcoords="offset points",
        )
    ax.set_xlabel("Human HK Ratio", fontsize=11)
    ax.set_ylabel("Mouse HK Ratio", fontsize=11)
    ax.set_title(
        f"Cross-Species Housekeeping Ratio Agreement\n"
        f"Spearman ρ = {rho_cross:.3f}, p = {pval_cross:.2e}",
        fontsize=12,
    )
    # Add identity line
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, "--", color="gray", alpha=0.5, linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    cross_path = OUT / "hk_ratio_cross_species.png"
    fig.savefig(cross_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {cross_path}")

    # ── Save results ───────────────────────────────────────────────────
    print("\n── Saving results ──")

    # Combined table
    results_df = pd.DataFrame({
        "cell_type": common,
        "hk_ratio_human": hk_ratio_h.loc[common].values,
        "hk_ratio_mouse": hk_ratio_m.loc[common].values,
        "residual_magnitude": residuals.loc[common].values,
    })
    results_df = results_df.merge(
        annot[["cell_type", "lineage", "progenitor"]],
        on="cell_type", how="left",
    )
    results_df = results_df.sort_values("residual_magnitude")
    csv_path = OUT / "hk_ratio_vs_residual.csv"
    results_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # JSON summary
    summary = {
        "analysis": "Thread 4A: Housekeeping gene ratio vs Procrustes rigidity",
        "hk_gene_source": "HOUNKPE_HOUSEKEEPING_GENES (MSigDB, 1129 symbols)",
        "hk_genes_in_ortholog_space": len(hk_ensembl),
        "hk_genes_in_centroids": len(hk_in_centroids),
        "unmapped_hk_genes": len(unmapped),
        "n_cell_types": len(common),
        "human_correlation": {
            "spearman_rho": round(rho_h, 4),
            "p_value": round(pval_h, 4),
            "significant_at_005": bool(pval_h < 0.05),
        },
        "mouse_correlation": {
            "spearman_rho": round(rho_m, 4),
            "p_value": round(pval_m, 4),
            "significant_at_005": bool(pval_m < 0.05),
        },
        "cross_species_agreement": {
            "spearman_rho": round(rho_cross, 4),
            "p_value": float(f"{pval_cross:.6e}"),
        },
    }
    json_path = OUT / "hk_ratio_results.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved: {json_path}")

    # ── Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nHousekeeping genes: {len(hk_in_centroids)} / {len(HOUNKPE_HK_GENES)} mapped to our 16,959-gene space")
    print(f"\nHuman:  Spearman ρ = {rho_h:.4f}, p = {pval_h:.4f}")
    print(f"Mouse:  Spearman ρ = {rho_m:.4f}, p = {pval_m:.4f}")
    print(f"Cross-species HK agreement: ρ = {rho_cross:.4f}, p = {pval_cross:.2e}")

    print("\n── Interpretation ──")
    if pval_h < 0.05 and rho_h < 0:
        print("HYPOTHESIS SUPPORTED (Human): Cell types with higher housekeeping gene")
        print("expression ratios have lower Procrustes residuals (are more rigid).")
        print("This is consistent with pleiotropic constraint: housekeeping genes are")
        print("under stronger purifying selection, stabilizing the transcriptomic programs")
        print("of cell types that rely more heavily on them.")
    elif pval_h < 0.05 and rho_h > 0:
        print("UNEXPECTED POSITIVE CORRELATION (Human): Cell types with higher HK ratios")
        print("are MORE diverged. This contradicts the pleiotropic constraint hypothesis.")
        print("Possible explanation: cell types with high HK ratios are less specialized,")
        print("and their non-HK genes are under weaker lineage-specific constraint.")
    else:
        print("HYPOTHESIS NOT SUPPORTED (Human): No significant correlation between")
        print("housekeeping gene ratio and evolutionary rigidity. The pleiotropic")
        print("constraint mechanism, at least measured this way, does not explain")
        print("cross-species rigidity differences.")

    if pval_m < 0.05:
        print(f"\nMouse REPLICATION: {'CONFIRMED' if np.sign(rho_m) == np.sign(rho_h) else 'CONTRADICTED'}.")
    else:
        print(f"\nMouse: Also not significant (consistent with human).")

    print(f"\nAll outputs saved to: {OUT}")
    print("=" * 70)


if __name__ == "__main__":
    main()
