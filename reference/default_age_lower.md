# Default (graded) age grid

Lower edges (in years) of the age groups: fine in infancy, where
immunity and maternal dynamics move fast, and coarse in adulthood. Group
edges are pinned at the conventional reporting boundaries – 0, 1, 2, 3,
5, 7, 10, 15, 20, 30, 40, 60 years – and the groups between two anchors
are of equal width, with the budget of groups spread across the anchors
in proportion to the width in `log(1 + age)`. Each interval gets at
least one group; the remainder is shared out largest-remainder.

## Usage

``` r
default_age_lower(max_age = 80, n_group = 53L)
```

## Arguments

- max_age:

  oldest age-group lower edge, in **years** (the absorbing top group).
  Must be a single finite number `>= 20`: below that the grid has too
  few anchors left to resolve the ages over which immunity develops.

- n_group:

  total number of age groups, the absorbing top group included. The
  default, 53, is one more than the 52 the old grid used, and is what
  the error figures above were measured at. Raising it refines the whole
  grid while keeping its shape, which is what a grid-convergence check
  wants.

## Value

numeric vector of age-group lower edges in years.

## Details

Log rather than linear width, because what the grid has to resolve is
the rise of immunity with age, and that is far closer to a function of
log age than of age. (Of `log(1 + age)` rather than `log(age)`, which is
infinite over the first year. The offset is in years, so the allocation
is not scale-free.) An equal-width grid at the same cost is four times
as far from the converged solution as this one, so the grading matters
more than the count.

This replaced a grid of fixed monthly / quarterly / yearly / 5-yearly
sections, halving the discretisation error at the same number of groups:
the largest departure across age bands falls from 10.7% to 4.5% and the
rms from 5.4% to 3.2%. The old grid over-resolved infancy (12 of its 52
groups in the first year, against 7 here) and under-resolved everything
above 15.

Refining removes that error: the departure from the converged profile
halves with each doubling of `n_group`, first order. What refinement
does **not** remove is the mean-field approximation itself – one
immunity value per stratum, where the IBM holds a spread of infection
histories at the same age – so a persistent difference from the IBM is
not evidence that the grid is too coarse. `validations/age-grid/run.R`
in `fleetcheck` measures both.

Band aggregation weights each age group by the exact fraction of its own
width that falls inside the band, so a band edge landing inside a group
(as one can, away from the anchors) apportions that group between the
two bands instead of handing it whole to one of them. The same is true
of intervention age targeting. Neither is snapped to the grid.

## Examples

``` r
default_age_lower()
#>  [1]  0.0000000  0.1428571  0.2857143  0.4285714  0.5714286  0.7142857
#>  [7]  0.8571429  1.0000000  1.2000000  1.4000000  1.6000000  1.8000000
#> [13]  2.0000000  2.2500000  2.5000000  2.7500000  3.0000000  3.4000000
#> [19]  3.8000000  4.2000000  4.6000000  5.0000000  5.5000000  6.0000000
#> [25]  6.5000000  7.0000000  7.7500000  8.5000000  9.2500000 10.0000000
#> [31] 11.2500000 12.5000000 13.7500000 15.0000000 16.6666667 18.3333333
#> [37] 20.0000000 22.5000000 25.0000000 27.5000000 30.0000000 33.3333333
#> [43] 36.6666667 40.0000000 44.0000000 48.0000000 52.0000000 56.0000000
#> [49] 60.0000000 65.0000000 70.0000000 75.0000000 80.0000000
# coarser top of the grid
default_age_lower(max_age = 60)
#>  [1]  0.00000  0.12500  0.25000  0.37500  0.50000  0.62500  0.75000  0.87500
#>  [9]  1.00000  1.20000  1.40000  1.60000  1.80000  2.00000  2.25000  2.50000
#> [17]  2.75000  3.00000  3.40000  3.80000  4.20000  4.60000  5.00000  5.50000
#> [25]  6.00000  6.50000  7.00000  7.75000  8.50000  9.25000 10.00000 11.00000
#> [33] 12.00000 13.00000 14.00000 15.00000 16.66667 18.33333 20.00000 22.00000
#> [41] 24.00000 26.00000 28.00000 30.00000 32.50000 35.00000 37.50000 40.00000
#> [49] 44.00000 48.00000 52.00000 56.00000 60.00000
# twice the resolution everywhere, same shape
length(default_age_lower(n_group = 105))
#> [1] 105
```
