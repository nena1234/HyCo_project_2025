library(clusterProfiler)
library(org.Hs.eg.db)

# Define gene symbols
top_200_up_genes_SK_for_KO <- read.table("~/Desktop/HyCo_Project_RNAseq_analysis/HyCo_RNAseq_analyzed-by-ERM/top_200_upregulated_gene_ids_SK-LMS-1.txt", header = FALSE, stringsAsFactors = FALSE)

colnames(top_200_up_genes_SK_for_KO ) <- c("Gene")

# Convert gene symbols to Entrez IDs (KEGG uses Entrez IDs)
gene_entrez_SK <- bitr(top_200_up_genes_SK_for_KO$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# Convert Entrez IDs to KEGG Orthology (KO) IDs
gene_kegg_SK <- bitr_kegg(gene_entrez_SK$ENTREZID, fromType = "ncbi-geneid", toType = "kegg", organism = "hsa")

# Merge the results
final_mapping_SK <- merge(gene_entrez_SK, gene_kegg_SK, by.x = "ENTREZID", by.y = "ncbi-geneid")

# View results
head(final_mapping)

# # SK-LMS-1 upregulated genes
# genes <- c("CA9",
#            "MME",
#            "TREM1",
#            "PLOD2",
#            "SLC2A1",
#            "IL1B",
#            "ERO1A",
#            "SRD5A3",
#            "ADM",
#            "ENO1",
#            "MMP13",
#            "SNTB1",
#            "CXCL1",
#            "SLC39A14",
#            "TMEM158",
#            "LINC01929",
#            "PGK1",
#            "IGFBP3",
#            "MGLL",
#            "NDNF",
#            "LINC02154",
#            "PDK1",
#            "C4orf3",
#            "BNIP3"
#            "SOD2"
#            "DSE"
#            "APLN"
#            "CPXM2"
#            "ENO2"
#            "ND1"
#            "GBE1"
#            "CA12"
#            "ALDOC"
#            "SPOCK1"
#            "THBS1"
#            "P4HA1"
#            "OGFRL1"
#            "MMP1"
#            "PIM1"
#            "NFKBIA"
#            "TGFA"
#            "DACT1"
#            "LOC105372338"
#            "CXCL6"
#            "C3"
#            "GAPDH"
#            "GPR176"
#            "PRDM1"
#            "EHF"
#            "TMEM45A"
#            "SEMA5A"
#            "DRAM1"
#            "BNIP3L"
#            "ABCA13"
#            "LOC124904862"
#            "SLC38A5"
#            "FHIP1A"
#            "DARS1"
#            "BMPER"
#            "ANKRD30B"
#            "PECAM1"
#            "NDRG1"
#            "TRPA1"
#            "WNT5A"
#            "CSF2"
#            "ALDOA"
#            "C1QTNF1"
#            "LOC105372497"
#            "COL8A1"
#            "PCGF5"
#            "BIRC3"
#            "CXCL8"
#            "PDGFC"
#            "LOX"
#            "NRP2"
#            "POU2F2-AS2"
#            "PFKP"
#            "BICC1
#            "CAV1
#            "TNIP1
#            "KDM3A
#            "BACH1
#            "VCAM1
#            "NEK7
#            "MXI1
#            "NOX4
#            "CCDC3
#            "BPIFB4
#            "PTX3
#            "FAM180A
#            "CTSS
#            "ZC3H12A
#            "IL6ST
#            "S100A4
#            "SGK1
#            "ICAM1
#            "TAF1D
#            "CCBE1
#            "ATP13A3
#            "SLC2A3
#            "BTG1
#            "TNFAIP8
#            "SLC4A4
#            "TENM2
#            "DNAH11
#            "SERPING1
#            "ANKRD12
#            "FKBP1A
#            "SPAG4
#            "CFL2
#            "SGIP1
#            "LINC03061
#            "FAM162A
#            "ND6
#            "MSC
#            "HIVEP2
#            "OSMR
#            "IL6
#            "SLC43A3
#            "BMAL2
#            "AK4
#            "LOC102724458
#            "SERPINE2
#            "GNG11
#            "MT1E
#            "AMPD3
#            "PRUNE2
#            "CASP1
#            "MKX
#            "LHFPL2
#            "SERPINB7
#            "CXCL3
#            "FAM227B
#            "ND2
#            "TMEM47
#            "HIF3A
#            "KITLG
#            "NFKBIZ
#            "RAB3B
#            "VLDLR
#            "SNORD3C
#            "TPI1
#            "SLC36A4
#            "MIR210HG_1
#            "QKI
#            "SCN9A
#            "ZNF292
#            "GPC6
#            "DARS1-AS1
#            "C5orf46
#            "IL1A
#            "OSGIN2
#            "PTPRH
#            "XIRP2
#            "ATP8B1
#            "CLIC4
#            "DIAPH2
#            "IL4I1
#            "PRR16
#            "INPP4B
#            "NT5DC1
#            "SH3D21
#            "ITPR1
#            "IL32
#            "RAI14
#            "CHRDL1
#            "NEGR1
#            "CAMK4
#            "RPL11
#            "ELL2
#            "NT5E
#            "NAP1L1
#            "STOM
#            "BHLHE40
#            "PLAT
#            "NOG
#            "TGFBR3
#            "LRRC15
#            "IQCD
#            "PLP2
#            "VEGFA
#            "RLF
#            "RIPK2
#            "PDGFB
#            "MSRB3-AS1
#            "GXYLT2
#            "PLOD1
#            "ARRDC3
#            "TRHDE
#            "CD47
#            "FKBP1C
#            "CHRAC1
#            "ANO6
#            "PTPN12
#            "ALCAM
#            "HIVEP1
#            "EVI2A
#            "MAP3K5
#            "NCAM1
#            "RPS13")  # Replace with your list
