# Cancer type considerations for COSMIC signature extraction

The
[`trinuc_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/trinuc_mutation_rates.md)
function uses mutational signature extraction to calculate relative
trinucleotide-specific SNV mutation rates in tumors. The
“signatures_to_remove” option allows some signatures to be excluded from
this analysis, which means each tumor will receive a weight of 0 for
these signatures, indicating that none of the tumor’s SNVs are
attributable to these signatures. This page describes the behavior of
the helper function
[`suggest_cosmic_signature_exclusions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/suggest_cosmic_signature_exclusions.md)
and the reasoning behind it.

As reported in [Alexandrov
2020](https://doi.org/10.1038/s41586-020-1943-3), and in COSMIC
v3.1/v3.2 mutational signature releases, some signatures are only
expected to appear in certain cancer types. For more reliable signature
extraction, consider excluding implausible signatures when running
[`trinuc_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/trinuc_mutation_rates.md).

There are also some signatures associated with various drug treatments
(SBS11, SBS31, SBS32, SBS35, SBS86, SBS87, SBS90), so you if you know
that your samples are treatment-naive or haven’t been exposed to the
implicated drugs, some or all of these signatures can be excluded.

Note that while some COSMIC signatures are attributed to sequencing
artifacts, you shouldn’t exclude these because cancereffectsizeR already
handles these signatures specially.

The
[`suggest_cosmic_signature_exclusions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/suggest_cosmic_signature_exclusions.md)
function will identify possible signature exclusions based on cancer
type and treatment status.

``` r
  suggest_cosmic_signature_exclusions(cancer_type = "BRCA", treatment_naive = TRUE)
  suggest_cosmic_signature_exclusions(cancer_type = "Kidney-RCC")
```

The cancer type recommendations are based on Extended Data Figure 5 of
Alexandrov 2020 and the COSMIC website (for signatures released after
the paper’s publication). The first two columns of the table below give
the labels accepted by the `cancer_type` argument.

Before excluding signatures, make sure your data set does not contain
tumors from multiple PCAWG categories. For example, TCGA HNSC (head and
neck cancer) includes oral cancers, which are listed separately here as
Oral-SCC, so excluding all signatures that do not appear in Head-SCC
(such as SBS29, tobacco chewing) would not be appropriate.

(You can access a text version of this table
[here](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/main/inst/extdata/COSMIC_v3.2_signatures_by_cancer_type.txt).)

| PCAWG | Applicable TCGA | Description | SBS signatures found |
|---:|---:|---:|---:|
| ALL | (none) | acute lymphoblastic leukemia | SBS1, SBS2, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS13, SBS18, SBS40, SBS84, SBS85, SBS86, SBS87 |
| Biliary-AdenoCA | CHOL | biliary adenocarcinoma | SBS1, SBS2, SBS3, SBS4, SBS5, SBS9, SBS10a, SBS10b, SBS10c, SBS10d, SBS12, SBS13, SBS15, SBS16, SBS17a, SBS17b, SBS18, SBS20, SBS21, SBS22, SBS24, SBS28, SBS29, SBS32, SBS38, SBS40, SBS42, SBS44, SBS86, SBS87 |
| Bladder-TCC | BLCA | bladder transitional cell carcinoma | SBS1, SBS2, SBS3, SBS5, SBS8, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS18, SBS19, SBS22, SBS29, SBS40, SBS86, SBS87, SBS88, SBS92 |
| Bone-Osteosarc | (none) | sarcoma, bone | SBS1, SBS2, SBS3, SBS5, SBS8, SBS10c, SBS10d, SBS13, SBS17a, SBS17b, SBS30, SBS40, SBS86, SBS87 |
| Bone-Other | (none) | Paper vague; presumably bone cartilaginous neoplasm osteoblastoma, bone osteofibrous dysplasia, bone neoplasm epithelioid | SBS1, SBS2, SBS5, SBS8, SBS10c, SBS10d, SBS13, SBS17a, SBS17b, SBS40, SBS86, SBS87 |
| Breast | BRCA | breast adenocarcinoma, ductal carcinoma in situ, breast lobular carcinoma | SBS1, SBS2, SBS3, SBS5, SBS8, SBS9, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS15, SBS17a, SBS17b, SBS18, SBS19, SBS21, SBS26, SBS29, SBS30, SBS34, SBS37, SBS39, SBS40, SBS41, SBS44, SBS86, SBS87, SBS90 |
| Cervix | CESC | cervical squamous cell carcinoma and adenocarcinoma | SBS1, SBS2, SBS5, SBS6, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS15, SBS17a, SBS17b, SBS18, SBS26, SBS28, SBS33, SBS40, SBS44, SBS86, SBS87 |
| CNS-GBM | GBM | CNS glioblastoma | SBS1, SBS2, SBS5, SBS10a, SBS10b, SBS10c, SBS10d, SBS11, SBS13, SBS15, SBS30, SBS37, SBS40, SBS42, SBS86, SBS87 |
| CNS-LGG | LGG | CNS lower-grade glioma | SBS1, SBS5, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS11, SBS19, SBS86, SBS87 |
| CNS-Medullo | (none) | CNS medulloblastoma | SBS1, SBS3, SBS5, SBS8, SBS10c, SBS10d, SBS18, SBS39, SBS40, SBS86, SBS87 |
| CNS-Oligo | (none) | CNS oligodenroglioma | SBS1, SBS5, SBS8, SBS10c, SBS10d, SBS40, SBS86, SBS87 |
| CNS-PiloAstro | (none) | CNS pilocytic astrocytoma | SBS1, SBS5, SBS10c, SBS10d, SBS19, SBS23, SBS40, SBS86, SBS87 |
| ColoRect-AdenoCA | COAD, READ | colorectal adenocarcinoma | SBS1, SBS2, SBS3, SBS5, SBS6, SBS9, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS15, SBS17a, SBS17b, SBS18, SBS20, SBS21, SBS26, SBS28, SBS30, SBS37, SBS40, SBS41, SBS44, SBS86, SBS87, SBS88, SBS94 |
| Eso-AdenoCA | ESCA | esophagus adenocarcinoma | SBS1, SBS2, SBS3, SBS4, SBS5, SBS6, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS13, SBS16, SBS17a, SBS17b, SBS18, SBS20, SBS26, SBS28, SBS34, SBS38, SBS40, SBS44, SBS86, SBS87, SBS90 |
| Eso-SCC | ESCA | esophagus squamous cell carcinoma | SBS1, SBS2, SBS5, SBS9, SBS10b, SBS10c, SBS10d, SBS13, SBS15, SBS17a, SBS17b, SBS18, SBS22, SBS36, SBS44, SBS86, SBS87, SBS90, SBS93 |
| Ewings | (none) | Ewing’s sarcoma | SBS1, SBS2, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS13, SBS18, SBS35, SBS40, SBS86, SBS87 |
| Eye-Melanoma | UVM | uveal melanoma | SBS1, SBS5, SBS10c, SBS10d, SBS40, SBS86, SBS87 |
| Head-SCC | HNSC | head-and-neck squamous cell carcinoma (note HNSC also includes some Oral-SCC) | SBS1, SBS2, SBS3, SBS4, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS15, SBS16, SBS17a, SBS17b, SBS18, SBS21, SBS30, SBS32, SBS33, SBS36, SBS38, SBS39, SBS40, SBS44, SBS86, SBS87, SBS88 |
| Kidney-ChRCC | (none) | kidney chromophobe renal cell carcinoma | SBS1, SBS2, SBS5, SBS10c, SBS10d, SBS13, SBS17a, SBS17b, SBS29, SBS40, SBS86, SBS87 |
| Kidney-Papillary | KIRP | papillary renal cell carcinoma | SBS1, SBS2, SBS5, SBS10c, SBS10d, SBS13, SBS86, SBS87 |
| Kidney-RCC | KIRC | kidney renal cell carcinoma | SBS1, SBS2, SBS4, SBS5, SBS10c, SBS10d, SBS12, SBS13, SBS14, SBS22, SBS26, SBS29, SBS40, SBS41, SBS86, SBS87 |
| Liver-HCC | LIHC | liver hepatocellular carcinoma | SBS1, SBS3, SBS4, SBS5, SBS6, SBS9, SBS10c, SBS10d, SBS12, SBS14, SBS16, SBS17a, SBS17b, SBS18, SBS19, SBS22, SBS23, SBS24, SBS26, SBS28, SBS29, SBS30, SBS31, SBS35, SBS37, SBS40, SBS86, SBS87 |
| Lung-AdenoCA | LUAD | lung adenocarcinoma | SBS1, SBS2, SBS3, SBS4, SBS5, SBS6, SBS9, SBS10c, SBS10d, SBS13, SBS15, SBS17a, SBS17b, SBS18, SBS28, SBS29, SBS40, SBS86, SBS87 |
| Lung-SCC | LUSC | lung squamous cell carcinoma | SBS1, SBS2, SBS3, SBS4, SBS5, SBS7a, SBS7b, SBS7c, SBS8, SBS10c, SBS10d, SBS13, SBS15, SBS18, SBS29, SBS30, SBS40, SBS44, SBS86, SBS87 |
| Lymph-BNHL | (none) | lymphoid mature B-cell lymphoma | SBS1, SBS2, SBS3, SBS5, SBS6, SBS9, SBS10c, SBS10d, SBS13, SBS17a, SBS17b, SBS34, SBS36, SBS37, SBS40, SBS42, SBS84, SBS85, SBS86, SBS87 |
| Lymph-CLL | (none) | lymphoid chronic lymphocytic leukemia | SBS1, SBS5, SBS9, SBS10c, SBS10d, SBS18, SBS40, SBS84, SBS85, SBS86, SBS87 |
| Myeloid-AML | LAML | acute myeloid leukemia | SBS1, SBS5, SBS10c, SBS10d, SBS18, SBS19, SBS31, SBS40, SBS86, SBS87 |
| Myeloid-MDS/MPN | (none) | myelodysplastic syndrome and myeloproliferative neoplasm | SBS1, SBS2, SBS5, SBS10c, SBS10d, SBS19, SBS23, SBS32, SBS86, SBS87 |
| Oral-SCC | (none) | oral squamous cell carcinoma (note some TCGA HNSC are oral) | SBS1, SBS2, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS10b, SBS10c, SBS10d, SBS13, SBS20, SBS29, SBS86, SBS87, SBS88 |
| Ovary-AdenoCA | OV | ovary adenocarcinoma | SBS1, SBS2, SBS3, SBS5, SBS7b, SBS7c, SBS7d, SBS8, SBS10c, SBS10d, SBS13, SBS18, SBS26, SBS35, SBS38, SBS39, SBS40, SBS41, SBS86, SBS87 |
| Panc-AdenoCA | PAAD | pancreatic adenocarcinoma | SBS1, SBS2, SBS3, SBS5, SBS6, SBS8, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS14, SBS15, SBS17a, SBS17b, SBS18, SBS20, SBS26, SBS28, SBS30, SBS31, SBS34, SBS37, SBS40, SBS44, SBS86, SBS87 |
| Panc-Endocrine | (none) | pancreatic neuroendocrine tumor | SBS1, SBS2, SBS3, SBS5, SBS6, SBS8, SBS9, SBS10c, SBS10d, SBS11, SBS13, SBS17a, SBS17b, SBS26, SBS30, SBS36, SBS39, SBS40, SBS86, SBS87 |
| Prost-AdenoCA | PRAD | prostate adenocarcinoma | SBS1, SBS2, SBS3, SBS5, SBS6, SBS8, SBS10c, SBS10d, SBS12, SBS13, SBS18, SBS33, SBS37, SBS39, SBS40, SBS41, SBS86, SBS87 |
| Sarcoma | SARC | sarcoma | SBS1, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS18, SBS40, SBS86, SBS87 |
| Sarcoma-bone | (none) | paper unclear; perhaps osteosarcoma and other sarcomas combined? | SBS1, SBS2, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS13, SBS18, SBS35, SBS40, SBS86, SBS87 |
| Skin-BCC | (none) | skin basal cell carcinoma | SBS1, SBS2, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS86, SBS87 |
| Skin-Melanoma | SKCM | skin melanoma | SBS1, SBS2, SBS3, SBS5, SBS7a, SBS7b, SBS7c, SBS7d, SBS9, SBS10c, SBS10d, SBS11, SBS13, SBS14, SBS17a, SBS17b, SBS31, SBS36, SBS38, SBS40, SBS86, SBS87 |
| Skin-SCC | (none) | skin squamous cell carcinoma | SBS1, SBS2, SBS5, SBS6, SBS7a, SBS7b, SBS7c, SBS7d, SBS10c, SBS10d, SBS13, SBS40, SBS86, SBS87 |
| SoftTissue-Leiomyo | (none) | leiomyosarcoma soft tissue | SBS1, SBS2, SBS5, SBS10c, SBS10d, SBS13, SBS17a, SBS17b, SBS30, SBS40, SBS86, SBS87 |
| SoftTissue-Liposarc | (none) | liposarcoma soft tissue | SBS1, SBS2, SBS3, SBS5, SBS10c, SBS10d, SBS13, SBS37, SBS40, SBS86, SBS87 |
| Stomach-AdenoCA | STAD | stomach adenocarcinoma | SBS1, SBS2, SBS3, SBS5, SBS6, SBS8, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS14, SBS15, SBS17a, SBS17b, SBS18, SBS20, SBS21, SBS26, SBS28, SBS34, SBS40, SBS41, SBS44, SBS86, SBS87, SBS93 |
| Thy-AdenoCA | THCA | thyroid low-grade adenocarcinoma | SBS1, SBS2, SBS5, SBS6, SBS10c, SBS10d, SBS13, SBS29, SBS40, SBS86, SBS87 |
| Uterus-AdenoCA | UCEC | uterus adenocarcinoma | SBS1, SBS2, SBS3, SBS5, SBS6, SBS7a, SBS7b, SBS7d, SBS10a, SBS10b, SBS10c, SBS10d, SBS13, SBS14, SBS15, SBS20, SBS21, SBS26, SBS28, SBS30, SBS40, SBS44, SBS86, SBS87 |
