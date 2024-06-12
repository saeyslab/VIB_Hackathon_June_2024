| Tool  | Method  | Data compatible/ benchmarked | Type of integration|Installation  | Details on usage  |Link to Github|other|
|---|---|---|---|---|---|---|---|
| SpatialGlue  | GNN  | Stereo-CITE-seq, SPOTS, 10x Visium + protein co-profiling, transcriptome-epigenome, generated data | linked data|PyPI (runs ok in conda)  |  rpy2 issues in env, all data should be in .h5ad |[link](https://github.com/JinmiaoChenLab/SpatialGlue)|returns attention weights for modalities| 
| MEFISTO  |factor analysis   | generated data, 10x Visium, no examples of real integration  | - |part of MOFA|-|[link](https://biofam.github.io/MOFA2/MEFISTO.html)|weights for factors (genes)|
| SLAT  | GNN  | aligning 2 Stereo-seq slices, 3D reconstruction from 4 Stereo-seq slices, 10x Xenium and 10x Visium alignment | cross-technology alignment, different slices  | docker, PyPI  |all data should be in .h5ad, requires manual preprocessing of the data |[link](https://github.com/gao-lab/SLAT)|notebooks with options for downstream analysis|

Table 1. Summary of discussed diagonal multi-omics integration tools
