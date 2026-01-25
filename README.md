# Alternating minimization algorithm for unlabeled sensing and linked linear regression
This repository contains the MATLAB implementation of the alternating minimization algorithm proposed in the paper [Alternating Minimization algorithm for unlabeled sensing and
linked linear regression](https://arxiv.org/pdf/2211.07621).

DOI: https://doi.org/10.1016/j.sigpro.2025.109927

## Repository Structure
The repository has 6 foldes: "Fig3", "Table1", "altMinProposed", "benchmarks", "dataSets", "misc". Below, I describe each of these folders and the files they contain.
* `Fig3:` To reproduce the results in Figure 3 of the paper. For example, to reproduce Fig. 3(a), run  the `fig3a.m` file in the folder `Fig3`.  The Python code for the stochastic EM benchmark algorithm in this folder was taken from the authors' own implementation at https://github.com/abidlabs/stochastic-em-shuffled-regression/tree/master,
* `Table1:` To reproduce the results in Table 1 of the paper.  
* `altMinProposed:` Code for the implementation of the proposed Alternating Minimization algorithm/
* `benchmarks:` Code for the benchmark methods proposed in other works.
* `dataSets:` The real datasets used for the results in Table 1.
* `misc:` Commonly used function throughout the code. 


## Instructions
* For linux systems, replace ' \ ' at the beginning of the the main.m files in the `figures` folder by '/'.
* To reproduce any of the results in [1], run the corerresponding file in the `Fig3` folder. For example, to reproduce Fig. 3(a), run  the `fig3a.m` file in the folder `Fig3`. 
* To reproduce the results on real-datasets in Table 1, please first install [CVX](https://cvxr.com/cvx/).
* The `sarcos` dataset (see Table 1 of Reference) is a large dataset, with a long run-time. The results for that dataset are provided in the MATLAB data file `sarcos.mat`

## References
[1].  Alternating Minimization algorithm for unlabeled sensing and linked linear regression. [URL](https://arxiv.org/pdf/2211.07621)

[2]. Ahmed Ali Abbasi, Abiy Tasissa, and Shuchin Aeron. R-local unlabeled sensing: A novel graph matching approach for multiview unlabeled sensing under local permutations. IEEE Open Journal of Signal Processing, 2:309–317, 2021.
[URL](https://ieeexplore.ieee.org/document/9440727)


## Feedback
Please email any feedback to the author Ahmed Ali Abbasi at aabbasi1@iastate.edu.

## Citation
Please cite the following papers on Unlabeled Sensing if you use this code and if these papers are relevant to your paper. Thank you!

```
1. @ARTICLE{9440727,
  author={Abbasi, Ahmed Ali and Tasissa, Abiy and Aeron, Shuchin},
  journal={IEEE Open Journal of Signal Processing}, 
  title={R-Local Unlabeled Sensing: A Novel Graph Matching Approach for Multiview Unlabeled Sensing Under Local Permutations}, 
  year={2021},
  volume={2},
  number={},
  pages={309-317},
  keywords={Sensors;Loss measurement;Covariance matrices;Signal processing algorithms;Signal processing;Task analysis;Minimization;Graph matching;multi-view;unlabeled sensing},
  doi={10.1109/OJSP.2021.3083479}}
```

```
2. @INPROCEEDINGS{9746201,
  author={Abbasi, Ahmed Ali and Tasissa, Abiy and Aeron, Shuchin},
  booktitle={ICASSP 2022 - 2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)}, 
  title={r-Local Unlabeled Sensing: Improved Algorithm and Applications}, 
  year={2022},
  volume={},
  number={},
  pages={5593-5597},
  keywords={Geometry;Linear systems;Signal processing algorithms;Signal processing;Minimization;Mathematical models;Sensors;Unlabeled sensing;unassigned distance geometry problem;proximal alternating minimization},
  doi={10.1109/ICASSP43922.2022.9746201}}
```

```
3.@article{ABBASI2025109927,
title = {Alternating minimization algorithm for unlabeled sensing and linked linear regression},
journal = {Signal Processing},
volume = {232},
pages = {109927},
year = {2025},
issn = {0165-1684},
doi = {https://doi.org/10.1016/j.sigpro.2025.109927},
url = {https://www.sciencedirect.com/science/article/pii/S0165168425000428},
author = {Ahmed Ali Abbasi and Shuchin Aeron and Abiy Tasissa}} 
```
