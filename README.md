# The Hidden Sources of Spurious Fusion Transcripts in Plants
## Overview
This repository contains the computational workflows, scripts, test datasets, and example files used in the study:
> **The Hidden Sources of Spurious Fusion Transcripts in Plants**

Fusion transcripts were initially discovered in cancer and are increasingly reported across many organisms, including plants. However, whether most reported plant fusion transcripts represent genuine biological events or technical/analytical artifacts has remained unclear.
In this study, we systematically re-evaluated fusion transcript detection in plants using multiple sequencing platforms, including:
- **PacBio Iso-Seq (IsoSeq)**
- **Nanopore direct RNA sequencing (dRNA-seq)**
- **Nanopore PCR-cDNA sequencing (cDNAseq)**
- **Short-read RNA-seq (RNAseq)**

Our results show that many reported plant fusion transcripts arise from previously underappreciated hidden sources, including:

- **Short homologous sequence (SHS)-mediated template switching**
- Alignment artifacts
- Genome assembly errors
- Annotation inconsistencies
- Reference-genome bias

This repository provides workflows to identify, classify, and visualize fusion transcripts in plants.

---

# Associated Manuscript
**The Hidden Sources of Spurious Fusion Transcripts in Plants** (under revision)

---

# Key Contributions
## Conceptual Advances
This work proposes an important distinction:
- **Fusion detection** ≠ **True biological fusion events**

Many detected fusion signals can arise during sequencing, reverse transcription, alignment, or annotation processes.

## Methodological Advances
This repository includes workflows for:
- Long-read fusion detection using an optimized **JAFFAL**
- Short-read fusion calling using **Arriba**
- Breakpoint classification (SHS / Joined / Unknown / Misaligned)
- Fusion network analysis
- Hybrid rice allelic/non-allelic fusion screening
- IGV-ready fusion visualization
- Overlap with chromatin/RNA 3D interaction datasets
- Genome annotation comparison

---

# Repository Structure
```text
.
├── README.md
├── workflow_commands.md
│
├── test_data/
│   Example datasets for running workflows
│
├── examples/
│   Example outputs used in manuscript analyses
│
├── figures/
│   Example visualization figures
│
├── environment/
│   software environments and R sessionInfo
│
├── scripts/
│   Main analysis scripts organized by module
│   ├── 01_preprocessing/
│   ├── 02_fusion_calling/
│   ├── 03_fusion_and_breakpoints_visualization/
│   ├── 04_breakpoint_classification/
│   ├── 05_fusion_network/
│   ├── 06_fusions_calling_short_reads/
│   ├── 07_hybrid_fusions/
│   ├── 08_others/
│   └── Stats_and_plot.R
````

---

# Workflow Summary
## 01_preprocessing
Preprocessing of:
* PacBio Iso-Seq raw subreads
* Nanopore direct RNA / cDNA reads
* Read correction and alignment

---

## 02_long-read_fusion_calling
Detection of putative fusion transcripts from long-read data using an optimized version of **JAFFAL**.

---

## 03_fusion_visualization
Generate IGV-compatible files for:

* Fusion transcripts
* Fusion breakpoints

---

## 04_breakpoint_classification
Classify fusion breakpoints into four categories:
* **SHS** (short homologous sequence)
* **Joined**
* **Unknown**
* **Misaligned**

---

## 05_fusion_network
Network analysis of recurrent fusion-partner genes using **igraph**.

---

## 06_fusion_calling_short_reads
Detect fusion transcripts from short-read RNA-seq using:
* STAR
* Arriba

---

## 07_hybrid_fusions
Identify candidate allelic and non-allelic fusion events in hybrid rice using the **trans-mate pair strategy**.

---

## 08_others
Additional analyses:

* Overlap with 3D genomic interaction data
* MH63 HiFi genome assembly
* Annotation-version comparisons

---


# Test Data
Small example datasets are provided under:

```text
test_data/
```

These allow users to test the pipeline structure and reproduce representative outputs.

For full-scale analyses, raw sequencing data used in the study are publicly available at NCBI: **PRJNA1291274**

---


# Citation

If you use this repository, please cite:

```text
The Hidden Sources of Spurious Fusion Transcripts in Plants
```

(doi to be added after publication)

---

# Contact
For questions, suggestions, or bug reports, please open an Issue or contact the corresponding authors.
- Xi-Tong Zhu	E-mail: z724@qq.com
- Yidan Ouyang	E-mail: diana1983941@mail.hzau.edu.cn
- Ling-Ling Chen	E-mail: llchen@gxu.edu.cn


---

# Final Note

This repository is not only a fusion-calling workflow. It is also a framework for asking a more fundamental question:

> **When we detect a fusion transcript, how confident can we be that a genuine fusion event has actually occurred?**

We hope this resource helps future studies avoid false positives and improve rigor in transcriptomics research.
