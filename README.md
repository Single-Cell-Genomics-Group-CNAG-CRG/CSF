# Decoding Immune Response in Leptomeningeal Disease through Single-Cell Sequencing of Cerebrospinal Fluid 🔬🧬🧠

## ✨ Overview 

This repository contains the code and analysis pipeline associated with the paper: [**"Decoding Immune Response in Leptomeningeal Disease through Single-Cell Sequencing of Cerebrospinal Fluid"**](https://www.biorxiv.org/content/10.1101/2025.01.27.634744v1). 📑🧐🔍

Our study utilizes single-cell TNA and TCR sequencing, bulk TCR sequencing, and spatial transcriptomics (CosMx platform) to investigate the immune landscape of cerebrospinal fluid (CSF) in patients with leptomeningeal disease (LMD). 🧪🧫🩸 

This repository provides all necessary scripts and instructions to reproduce the analysis performed in the paper. 📝💻📊

## 📂 Repository Structure

```
├── data/                # How to get the data
├── notebooks/           # Jupyter notebooks for step-by-step analysus
│   ├── python/          # Python-based analysis notebooks
│   │   ├── environment.yml  # Conda environment for Python
│   ├── R/               # R-based analysis notebooks
│   │   ├── environment.yml  # Conda environment for R
├── README.md           # This README file
```

## 🚀 Getting Started

### ✅ Clone the repository

```bash
git clone https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/CSF.git
cd REPO_NAME
```

### ✅  Set up the environments

To set up the Python environment:
```bash
cd notebooks/python
conda env create -f environment.yml
conda activate scRNAseq_python
```

To set up the R environment:
```bash
cd notebooks/R
conda env create -f environment.yml
conda activate scRNAseq_R
```

### 📊 Usage

Run the analysis step by step using Jupyter notebooks:

1. Open Jupyter Notebook:
   ```bash
   jupyter notebook
   ```
2. Navigate to the `notebooks/python/` or `notebooks/R/` directory depending on the script

## 📝 Citation

If you use this code in your research, please cite our paper:

> [Decoding immune response in leptomeningeal disease through single-cell sequencing of cerebrospinal fluid](https://www.biorxiv.org/content/10.1101/2025.01.27.634744v1)

## 📞 Contact

For questions or issues, please [open an issue](https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/CSF/issues/new/choose) in this repository or contact [Paula Nieto](mailto:paula.nieto@cnag.eu). 

---
