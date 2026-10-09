# Parallel lapply helper

Parallel lapply helper

## Usage

``` r
.lapply_cores(X, FUN, n_cores = 1L)
```

## Arguments

- X:

  List/vector to iterate over.

- FUN:

  Function applied to each element.

- n_cores:

  Number of workers (forked on unix, PSOCK otherwise).

## Value

A list.
