# QIIME 2 Pipeline for Simpson's Alpha Diversity Analysis with Rarefaction

This document describes the complete workflow used to compute Simpson's alpha diversity for our tadpole gut microbial communities. The pipeline includes:

1. Downloading a pre-filtered feature table artifact.
2. Exporting the feature table to a BIOM file.
3. Importing the BIOM file back into QIIME 2.
4. Rarefying the feature table to standardize sequencing depth.
5. Computing Simpson’s alpha diversity from the rarefied feature table.
6. Exporting the final diversity metrics for downstream analysis.

Each step is documented with the associated QIIME 2 commands and an explanation of its purpose.

## Prerequisites

- QIIME 2 (version 2024.10.1 or later, to activate: conda activate qiime2-amplicon-2024.10)
- Access to the pre-filtered feature table artifact from Box
- Basic command-line proficiency

## Pipeline Steps

### 1. Download the Pre-filtered Feature Table

Our pre-filtered feature table is stored as a QIIME 2 artifact on Box. The file we downloaded is located at:

20250205_M008009_report/analysis/filtered/table_filtered.qza


### 2. Export the Feature Table as a BIOM File

Extract the data from the QIIME 2 artifact:

```
qiime tools export \
  --input-path table_filtered.qza \
  --output-path exported_feature_table
```

This command creates a directory exported_feature_table containing the file feature-table.biom.
### 3. Import the BIOM File into QIIME 2

If further processing within QIIME 2 is needed, import the BIOM file:

```
qiime tools import \
  --input-path exported_feature_table/feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature_table.qza
```

### 4. Rarefy the Feature Table

Rarefaction subsamples each sample in the feature table to a common sequencing depth, which normalizes differences in sequencing effort. Choose an appropriate sampling depth (for example, 10,000 sequences per sample):

```
qiime feature-table rarefy \
  --i-table feature_table.qza \
  --p-sampling-depth 10000 \
  --o-rarefied-table rarefied_table.qza
```

This produces rarefied_table.qza, a rarefied version of your feature table.
### 5. Compute Simpson's Alpha Diversity

Using the rarefied feature table, calculate Simpson's alpha diversity:

```
qiime diversity alpha \
  --i-table rarefied_table.qza \
  --p-metric simpson \
  --o-alpha-diversity rarefied_simpson_vector.qza
```

This command computes Simpson's diversity for each sample and outputs the results as simpson_vector.qza.

```
qiime diversity alpha-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table rarefied_table.qza \
  --p-metric faith_pd \
  --o-alpha-diversity rarefied_faith_pd_vector.qza
```

### 6. Export the Simpson Diversity Results 

To export the Simpson diversity metrics to a TSV file for further analysis in R or another tool, run:

```
qiime tools export \
  --input-path rarefied_simpson_vector.qza \
  --output-path exported_simpson_rarefied
```

### Bonus - adding in Beta diversity calculations:

qiime diversity beta \
  --i-table rarefied_table.qza \
  --p-metric braycurtis \
  --o-distance-matrix rarefied_braycurtis_distance_matrix.qza

qiime tools export \
  --input-path rarefied_braycurtis_distance_matrix.qza \
  --output-path exported_braycurtis_matrix

  ### Super extra bonus - taxa bar plots:

  qiime taxa barplot \
  --i-table rarefied_table.qza \
  --i-taxonomy Taxa/taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-barplots-rarefied.qzv

qiime tools export \
  --input-path taxa-barplots-rarefied.qzv \
  --output-path exported_taxa_barplot_data
