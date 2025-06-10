# Hua-Lab-16S-Microplastics
Hua Lab Microplastics

Data analysis and plotting for "Microplastic Fibers in Aquatic Ecosystems: Implications for Amphibians and Host-Parasite Dynamics"

## Analysis Breakdown

This project is divided into three main analysis sections. For a detailed explanation of the methodology and results for each section, please see the `README.md` file within the respective directory:

*   **[Alpha Diversity Analysis](./AlphaAnalysis/README.md)**: Investigates the effects of treatments on within-sample diversity (richness and evenness).

*   **[Beta Diversity Analysis](./BetaAnalysis/README.md)**: Examines the differences in microbial community structure between treatment groups.

*   **[Taxonomic Composition Analysis](./TaxaAnalysis/README.md)**: Visualizes the relative abundance of different bacterial phyla across treatments.


# First time setup in R shell:
make sure you set the correct path for the working directory:
```
setwd(~/tadpole-data)
```

Then install the dependencies defined in this project (renv.lock file)

```
renv::restore()
```

If you need to add more dependenciues while you are working then do:

```
renv::isntall(c('dependency1', 'dependecy2'))
```

and those will get add to the lock file.
