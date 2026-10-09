# Methods for nnG objects

Methods for nnG objects

## Usage

``` r
# S3 method for class 'nnG'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'nnG'
print(x, ...)

# S3 method for class 'nnG'
summary(object, ...)

# S3 method for class 'nnG'
plot(x, y = NULL, ...)
```

## Arguments

- x, object:

  An `nnG` object.

- row.names, optional, ...:

  Further arguments (ignored or passed on).

- y:

  Ignored; present for compatibility with the plot generic.

## Value

`as.data.frame` returns the underlying data frame, `summary` a list of
class `summary.nnG`, and `print`, `plot` return their input invisibly.
