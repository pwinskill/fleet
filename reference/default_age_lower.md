# Default (graded) age grid

Lower edges (in years) of the age groups. Fine in infancy, where
immunity and maternal dynamics move fast, and coarse in adulthood. Band
aggregation in the outputs weights each age group by the exact fraction
of its own width that falls inside the band, so a band edge landing
inside a group (as it can on a coarse custom grid) apportions that group
between the two bands instead of handing it whole to one of them; the
default grid places edges at 2, 5, 10 and 15 years, so the usual
rendering bands fall on group boundaries exactly and every weight is 0
or 1.

## Usage

``` r
default_age_lower(max_age = 80)
```

## Arguments

- max_age:

  oldest age-group lower edge, in **years** (absorbing top group). Must
  be a single finite number `>= 20`: the grid is graded up to a 5-yearly
  section starting at 15, so there is no room for an absorbing top group
  below 20.

## Value

numeric vector of age-group lower edges in years.

## Examples

``` r
default_age_lower()
#>  [1]  0.00000000  0.08333333  0.16666667  0.25000000  0.33333333  0.41666667
#>  [7]  0.50000000  0.58333333  0.66666667  0.75000000  0.83333333  0.91666667
#> [13]  1.00000000  1.25000000  1.50000000  1.75000000  2.00000000  2.25000000
#> [19]  2.50000000  2.75000000  3.00000000  3.25000000  3.50000000  3.75000000
#> [25]  4.00000000  4.25000000  4.50000000  4.75000000  5.00000000  6.00000000
#> [31]  7.00000000  8.00000000  9.00000000 10.00000000 11.00000000 12.00000000
#> [37] 13.00000000 14.00000000 15.00000000 20.00000000 25.00000000 30.00000000
#> [43] 35.00000000 40.00000000 45.00000000 50.00000000 55.00000000 60.00000000
#> [49] 65.00000000 70.00000000 75.00000000 80.00000000
# coarser top of the grid
default_age_lower(max_age = 60)
#>  [1]  0.00000000  0.08333333  0.16666667  0.25000000  0.33333333  0.41666667
#>  [7]  0.50000000  0.58333333  0.66666667  0.75000000  0.83333333  0.91666667
#> [13]  1.00000000  1.25000000  1.50000000  1.75000000  2.00000000  2.25000000
#> [19]  2.50000000  2.75000000  3.00000000  3.25000000  3.50000000  3.75000000
#> [25]  4.00000000  4.25000000  4.50000000  4.75000000  5.00000000  6.00000000
#> [31]  7.00000000  8.00000000  9.00000000 10.00000000 11.00000000 12.00000000
#> [37] 13.00000000 14.00000000 15.00000000 20.00000000 25.00000000 30.00000000
#> [43] 35.00000000 40.00000000 45.00000000 50.00000000 55.00000000 60.00000000
```
