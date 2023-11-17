# ClusterCirc_SPSS

This repository contains the SPSS codes to perform ClusterCirc (cc_data and cc_simu) in SPSS
or the alternative freeware PSPP. ClusterCirc codes for exemplary data and the resulting output 
are also included.

## Description of ClusterCirc

ClusterCirc is a clustering method designed for data with circular
structure. It can be used to find item clusters with optimal circumplex
spacing as an alternative to other clustering techniques like
conventional cluster analysis.

ClusterCirc can be applied to raw data or on item loadings on two
orthogonal factors or components from principal component analysis,
exploratory or confirmatory factor analysis. If ClusterCirc is used on
raw data, principal component analysis is performed before ClusterCirc
to yield loadings on two unrotated components.

The ClusterCirc algorithm uses item loadings and translates them into
angles in a circular arrangement by trigonometric conversion.
ClusterCirc then sorts items into clusters that yield optimal circumplex
spacing.

Optimal circumplex spacing for item clusters is given if clusters are
evenly distributed across the circle (equal spacing between clusters)
and if items are clustered closely around their cluster centroid
(minimal within-cluster spacing). Spacing coefficients are computed to
assess circumplex spacing of items, clusters, and the overall data.
Range of all ClusterCirc coefficients: 0-1 (0 = perfect circumplex
spacing).

## Using ClusterCirc in SPSS

ClusterCirc can be performed in SPSS by using the syntax files from this repository.

1.  **cc_data: (SPSS_cc_data_INSERT.sps)**  
    Main function. Sorts items of your dataset into clusters with optimal
    circumplex spacing. Spacing coefficients are computed for the suggested clustering.

    Follow the instructions in the comments above the code. You need to insert the following
    parameters of your data in the code:
    - working directory/path
    - p = Desired number of clusters
    - m = Number of variables
    - n = Sample size
    - data = data file on which ClusterCirc should be performed
  
3.  **cc_simu: (SPSS_cc_simu_INSERT.sps)**  
    Can be used to assess circumplex fit of the dataset.
    The function uses the specifications of the data and creates samples
    from a population with perfect circumplex spacing of clusters.
    Results for the dataset (spacing coefficients from cc_data)
    are compared to results from cc_simu to evaluate circumplexity in the data.
    cc_simu can only be used after performing cc_data.

    Follow the instructions in the comments above the code. You need to insert the following
    parameters of your data in the code:
    - working directory/path
    - n = Sample size
    - samples = Number of samples for the simulation. Recommendation: 100-500. Needs to be inserted two times.
    - m+2 = Number of items in the data + 2. Needs to be inserted three times.
  
### Usage for exemplary data:
- SPSS_cc_data_example.sps
- SPSS_cc_simu_example.sps

You only need to insert the working directory of your device to perform ClusterCirc with the example codes.

## Citing ClusterCirc:

The manuscript that presents Cluster-Circ has been submitted to a peer-
reviewed journal. When using Cluster-Circ, please cite the preprint version 
at the current stage of the publication process: 

https://psyarxiv.com/yf37w/

(Note to reviewers: Please don't open the preprint article to ensure double-blind
peer-review).

## ClusterCirc in R

ClusterCirc can also be performed in R by installing and loading the
associated R package from github: 

https://github.com/anweide/ClusterCirc
