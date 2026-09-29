# Calculate SIs at gene level under pairwise epistasis model

The genes are assumed to not overlap in any ranges (the calling function
checks for this)

## Usage

``` r
pairwise_gene_epistasis(cesa, genes, samples, conf)
```

## Arguments

- cesa:

  CESAnalysis

- genes:

  two-length vector of gene names

- samples:

  validated sample subset, as from select_samples()

- conf:

  confidence level on (0, 1)
