# Analyses for 'Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays'

## Overview
The analyses here, were performed to benchmark the use of MPAC (Malinois with Parallel Aggregated Cross-validation) for non-coding variant effect prediction. Full descriptions of each analysis can be found in the README of that directory. There are a number of other repositories associated with the manuscript that we list below.

## Related Repositories
| Repository | Description |
|---|---|
| https://github.com/sjgosai/boda2 | Code for original Malinois model (Gosai et al 2024) published at: https://www.nature.com/articles/s41586-024-08070-z |
| https://github.com/Reilly-Lab-Yale/coda_mpac | Fork of Boda2 and primary code base for MPAC variant effect predictions, includes MPAC-specific command line tools for Indel and Haplotype predictions |
| https://huggingface.co/saarantras1/MPAC | Huggingface mirror of MPAC model with checkpoints converted to safetensors |
| https://github.com/john-c-butts/MPAC | Repository for analyses related to: Benchmarking, ClinVar, COSMIC, dELS saturation mutagenesis, seqlet definition, GTEx mechanisms |
| https://github.com/Reilly-Lab-Yale/MPAC_gnomAD_and_satmut | Repository for analyses related to gnomAD and promoter saturation mutagenesis |
| https://github.com/Reilly-Lab-Yale/predictions_plusplus | Code for additional prediction tasks, including dELS saturation mutagenesis |
| https://doi.org/10.5281/zenodo.15178434 | Contains all predictions, supplementary datasets, and model weights used in 'Identifying non-coding variant effects at scale via machine learning models of cis-regulatory reporter assays' |
| https://huggingface.co/datasets/saarantras1/mpac-predictions | Huggingface mirror of predictions available at https://doi.org/10.5281/zenodo.15178434 |
