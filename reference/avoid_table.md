# Avoidance table

`avoid_table(B, s)[k + 1] = choose(B - k, s) / choose(B, s)`, with value
zero whenever `B - k < s`.

## Usage

``` r
avoid_table(B, s)
```

## Arguments

- B, s:

  Non-negative integers.

## Value

Numeric vector of length `B + 1`.
