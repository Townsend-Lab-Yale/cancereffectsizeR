# Create an independent copy of a CESAnalysis

Used internally to "copy" CESAnalysis objects while keeping memory use
to a minimum.

## Usage

``` r
copy_cesa(cesa)
```

## Arguments

- cesa:

  CESAnalysis

## Details

The trick is to use data.table's copy function on all data.tables (and
lists of data.tables) within CESAnalysis slots. (If you just call copy
on the whole object, the data tables won't be handled in a
memory-efficient way. And if you call copy on non-data.tables, it's
actually less efficient since it forces an immediate full copy instead
of the usual copy-on-modify.)
