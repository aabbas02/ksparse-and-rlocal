


# Alternating minimization algorithm for unlabeled sensing and linked linear regression
This repository contains the MATLAB implementation of the alternating minimization algorithm proposed in the paper [Alternating Minimization algorithm for unlabeled sensing and
linked linear regression](https://arxiv.org/pdf/2211.07621).

DOI: https://doi.org/10.1016/j.sigpro.2025.109927

## Repository Structure
The repository has 6 foldes: "Fig3", "Table1", "altMinProposed", "benchmarks", "dataSets", "misc". Below, I describe each of these folders and the files they contain.
* Fig3:  To reproduce the results in Figure 3 of the paper. For example, to reproduce Fig. 3(a), run  the `fig3a.m` file in the folder `Fig3`. 
* Table1: To reproduce the results in Table 1 of the paper.  
* altMinProposed: Code for the implementation of the proposed Alternating Minimization algorithm/
* benchmarks: Code for the benchmark methods proposed in other works.
* dataSets: The real datasets used for the results in Table 1.
* misc: Commonly used function throughout the code. 


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
