# Topological Data Analysis: Simplicial Complex

[![DOI](https://zenodo.org/badge/1031948864.svg)](https://doi.org/10.5281/zenodo.22163397) [![CRAN status](https://www.r-pkg.org/badges/version/SimplicialComplex)](https://cran.r-project.org/package=SimplicialComplex) <a href="https://CRAN.R-project.org/package=SimplicialComplex" target="_blank" rel="noreferrer"> <img src="https://cranlogs.r-pkg.org/badges/grand-total/SimplicialComplex" alt="CRAN downloads" width="100" height="20"/> </a>

`SimplicialComplex` is a user-friendly Topological Data Analysis (TDA) package written entirely in R. While most TDA libraries (Dionysus, PHAT, GUDHI) are developed in Python and C++, implementing simplicial complexes natively in R makes them directly compatible with the rich ecosystem of statistical methods R already offers.

## Features

- **Simplicial complexes**: Build Vietoris–Rips, Alpha, Cech, Witness, Cubical, and Flood complexes from point clouds, or define abstract simplicial complexes by hand.
- **Topological invariants**: Faces, boundary matrices, Betti numbers, and the Euler characteristic.
- **Persistent homology**: filtrations, boundary-matrix reduction, persistence pairs, persistence diagrams, and persistence landscape.
- **Statistical**: Wasserstein distance & bottleneck distance were built in the latest version, with matching plot.
- **Examples**: Full worked examples in `inst/example`.

## Playground

Try the [interactive playground](https://tf3q5u-0-0.shinyapps.io/simplicialcomplex/) to get familiar with all the concepts used in TDA.

## References

- Zomorodian, A., & Carlsson, G. (2004). Computing persistent homology. *Proceedings of the Twentieth Annual Symposium on Computational Geometry*, 347–356.
- Chazal, F., & Michel, B. (2021). An introduction to topological data analysis: Fundamental and practical aspects for data scientists. *Frontiers in Artificial Intelligence*, 4, 667963.
- Graf, F., Pellizzoni, P., Uray, M., Huber, S., & Kwitt, R. (2025). The Flood Complex: Large-scale persistent homology on millions of points. *Advances in Neural Information Processing Systems*, 38.
- Otter, N., Porter, M. A., Tillmann, U., Grindrod, P., & Harrington, H. A. (2017). A roadmap for the computation of persistent homology. EPJ data science, 6(1), 17.