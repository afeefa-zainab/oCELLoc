<div align="center">
  <img src="man/figures/logo.png" height="100px"/>
  <h1>oCELLoc</h1>
</div>

**oCELLoc** is an R package designed to predict cell types in spatial transcriptomics and scRNA-seq data. It leverages shrinkage methods to provide insights into spatial patterns of gene expression.

## Package Metadata

- **Package**: oCELLoc
- **Title**: Predicts Suitable Cell Types in Spatial Transcriptomics and scRNA-seq Data
- **Version**: 1.0.0
- **Authors**: 
  - Afeefa Zainab (`afeeffazainab@gmail.com`) - Author and Creator
  - Vladyslav Honcharuk (`vladyslav.s.honcharuk@gmail.com`) - Author
  - Alexis Vandenbon (`alexisvdb@infront.kyoto-u.ac.jp`) - Author
- **Description**: Picks the suitable cell types in spatial and scRNA-seq data using shrinkage methods. The package includes curated reference gene expression profiles for human and mouse cell types, facilitating immediate application to common spatial or scRNA datasets. Additionally, users can input custom reference data to support tissue- or experiment-specific analyses.
- **License**: MIT
- **Encoding**: UTF-8

## Installation

You can install the package from GitHub using:

```r
devtools::install_github("afeefa-zainab/oCELLoc")
```

## Citation

If you use this package in a publication, please cite:

Afeefa Zainab, Vladyslav Honcharuk, Alexis Vandenbon (2025). oCELLoc: Automated Cell Type Assignment in Transcriptomics Data Using Reference Filtering. https://doi.org/10.64898/2025.12.11.693812

