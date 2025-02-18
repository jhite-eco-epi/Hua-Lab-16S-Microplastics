# Hua-Lab-16S-Microplastics
Hua Lab Microplastics

Data analysis and plotting for "Microplastic Fibers in Aquatic Ecosystems: Implications for Amphibians and Host-Parasite Dynamics"


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
