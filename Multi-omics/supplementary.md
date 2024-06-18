| Tool  | Method  | Data compatible/ benchmarked | Type of integration|Installation  | Details on usage  |Link to Github|other|
|---|---|---|---|---|---|---|---|
| SpatialGlue  | GNN  | Stereo-CITE-seq, SPOTS, 10x Visium + protein co-profiling, transcriptome-epigenome, generated data | linked data|PyPI (runs ok in conda)  |  rpy2 issues in env, all data should be in .h5ad |[link](https://github.com/JinmiaoChenLab/SpatialGlue)|returns attention weights for modalities| 
| MEFISTO  |factor analysis   | generated data, 10x Visium, no examples of real integration  | - |part of MOFA|-|[link](https://biofam.github.io/MOFA2/MEFISTO.html)|weights for factors (genes)|
| SLAT  | GNN  | aligning 2 Stereo-seq slices, 3D reconstruction from 4 Stereo-seq slices, 10x Xenium and 10x Visium alignment | cross-technology alignment, different slices  | docker, PyPI  |all data should be in .h5ad, requires manual preprocessing of the data |[link](https://github.com/gao-lab/SLAT)|notebooks with options for downstream analysis|

Table 1. Summary of discussed diagonal multi-omics integration tools

## Datasets

Data used in [STalign](https://www.nature.com/articles/s41467-023-43915-7) paper: https://www.nature.com/articles/s41467-023-43915-7#data-availability

Data used in [CAST](https://www.biorxiv.org/content/10.1101/2023.08.13.552987v1.full). Link to data doesn't work.

- SPOTS with the 10x Visium technology capturing whole transcriptomes and extracellular proteins https://doi.org/10.1038/s41587-022-01536-3, GSE198353. High-resolution images (https://figshare.com/account/home#/projects/143019)
- Stereo-CITE-seq spatial transcriptomics + proteomics (https://doi.org/10.1101/2023.04.28.538364)
- spatial transcriptomics + DVP proteomics (https://doi.org/10.1038/s41593-022-01097-3)
- Spatial-ATAC-RNA-seq (https://doi.org/10.1038/s41586-023-05795-1)
- Cite-seq, proteogenomics (https://doi.org/10.1016/j.cell.2021.12.018)
- spatial CITE-seq transcriptomics+proteomics (https://doi.org/10.1038/s41587-023-01676-0)
- Benchmark datasets for 3D mass spec imaging (=2D Mass spec imaging on adjacent sections) (https://academic.oup.com/gigascience/article/4/1/s13742-015-0059-4/2707545)
- https://doi.org/10.1038/s41467-023-43105-5 (suppl table 1, collection of publicly available datasets from different studies)
- spatial-ATAC and the spatial RNA-seq (MISAR-seq, https://doi.org/10.1038/s41592-023-01884-1)
- Mass spec imaging + spatial transcriptomics (Visium): https://www.nature.com/articles/s41587-023-01937-y (see data availability, e.g. https://data.mendeley.com/datasets/w7nw4km7xd/1, sma zip file)


## Interesting papers

- [Integration of Multiple Spatial Omics Modalities Reveals Unique Insights into Molecular Heterogeneity of Prostate Cancer](https://www.biorxiv.org/content/10.1101/2023.08.28.555056v1.full) Spatial transcriptomics and Mass spec imaging were performed on adjacent sections, and registered via their respective H&E images. The datasets are not publically available.
- [Search and Match across Spatial Omics Samples
  at Single-cell Resolution](https://www.biorxiv.org/content/10.1101/2023.08.13.552987v1.full)
- https://frontlinegenomics.com/a-guide-to-multi-omics-integration-strategies/