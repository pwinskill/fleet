# Default age grid

Lower edges, in years, of the default age groups: fine in childhood,
where immunity and maternal immunity change fast, and coarse in
adulthood.

## Usage

``` r
default_age_lower(max_age = 80, n_group = 118L)
```

## Arguments

- max_age:

  lower edge, in years, of the oldest age group, which is open-ended. A
  single finite number of at least 20: below that the grid has too few
  anchors left to resolve the ages over which immunity develops.

- n_group:

  total number of age groups, the open-ended top group included. A whole
  number of at most 1e6, and at least one more than the number of anchor
  intervals below `max_age` (14 for `max_age = 80`). Changing it refines
  or coarsens the whole grid while keeping its shape, as a
  grid-convergence check needs; run time is roughly proportional to it.

## Value

A numeric vector of age-group lower edges in years, from 0 to `max_age`.

## Details

Group edges are pinned at the conventional reporting ages, 0, 1, 2, 3,
5, 7, 10, 15, 20, 30, 40 and 60 years, so that the usual output bands
fall on group boundaries, and at 21, so that whole groups make up the
year, \[20, 21), from which the IBM draws maternal immunity. Groups
between two anchors are of equal width. The groups are shared among the
anchor intervals in proportion to each interval's width in
`log(1 + age)`, weighted 1.25 below 5 years and 0.5 from 21 years; each
interval gets at least one group, and the remainder is allocated by
largest remainder.

The width is logarithmic because what the grid has to resolve is the
rise of immunity with age, which is far closer to a function of log age
than of age. (`log(1 + age)` rather than `log(age)`, which is infinite
over the first year; the offset is in years, so the allocation is not
scale-free.) The grading matters more than the count: at 53 groups an
equal-width grid is four times as far from the converged solution as a
log-graded one (rms over the age profile), and a grid of fixed monthly,
quarterly, yearly and 5-yearly sections, which over-resolves infancy and
under-resolves everything above 15, 1.7 times as far. The grid's error
is largest in clinical and severe incidence among young children at high
transmission; each age band's error comes from the groups inside it, and
the groups above 21 years have no measurable effect on any child's
incidence. The weights move groups to where the error arises, and match
the unweighted grid's accuracy with about a quarter fewer groups.

The default, 118 groups, is the smallest grid on which every falciparum
claim in `fleetcheck` passes, not the most accurate; the deciding claim
is the clinical age profile's 3-5 year band at EIR 120, which needs 14
groups between 3 and 5 years. `fleet`'s age profile converges at first
order: its departure from the converged profile is fleet-low and halves
with each doubling of `n_group`. On the default grid the largest
departure of the EIR 20 clinical age profile is 3.2%, in the oldest
band, which the weights coarsen, and the rms 1.6%; on 53 groups they are
6.1% and 3.2%. Against the IBM median, severe incidence at EIR 120 is 2%
low on the default grid and 7% low on 53 groups.

`n_group = 53` runs over twice as fast and is adequate where a result
does not rest on the level of severe incidence, or of young children's
incidence at high transmission; `n_group = 236` halves the grid's error
for twice the cost. Refinement does not remove the mean-field
approximation itself, one immunity value per stratum where the IBM holds
a spread of infection histories at the same age, so a persistent
difference from the IBM is not evidence that the grid is too coarse.
`fleetcheck` measures both, in `validations/age-grid/run.R` and the
`age-profile-clinical` and `severe-allage-eir` claims.

Output bands and intervention age targets that do not fall on a group
edge are apportioned by the fraction of each group's width they cover;
neither is snapped to the grid.

## See also

[`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md),
which takes the grid as `age_lower`.

## Examples

``` r
default_age_lower()
#>   [1]  0.00000000  0.04545455  0.09090909  0.13636364  0.18181818  0.22727273
#>   [7]  0.27272727  0.31818182  0.36363636  0.40909091  0.45454545  0.50000000
#>  [13]  0.54545455  0.59090909  0.63636364  0.68181818  0.72727273  0.77272727
#>  [19]  0.81818182  0.86363636  0.90909091  0.95454545  1.00000000  1.07142857
#>  [25]  1.14285714  1.21428571  1.28571429  1.35714286  1.42857143  1.50000000
#>  [31]  1.57142857  1.64285714  1.71428571  1.78571429  1.85714286  1.92857143
#>  [37]  2.00000000  2.10000000  2.20000000  2.30000000  2.40000000  2.50000000
#>  [43]  2.60000000  2.70000000  2.80000000  2.90000000  3.00000000  3.14285714
#>  [49]  3.28571429  3.42857143  3.57142857  3.71428571  3.85714286  4.00000000
#>  [55]  4.14285714  4.28571429  4.42857143  4.57142857  4.71428571  4.85714286
#>  [61]  5.00000000  5.25000000  5.50000000  5.75000000  6.00000000  6.25000000
#>  [67]  6.50000000  6.75000000  7.00000000  7.33333333  7.66666667  8.00000000
#>  [73]  8.33333333  8.66666667  9.00000000  9.33333333  9.66666667 10.00000000
#>  [79] 10.50000000 11.00000000 11.50000000 12.00000000 12.50000000 13.00000000
#>  [85] 13.50000000 14.00000000 14.50000000 15.00000000 15.62500000 16.25000000
#>  [91] 16.87500000 17.50000000 18.12500000 18.75000000 19.37500000 20.00000000
#>  [97] 20.50000000 21.00000000 22.80000000 24.60000000 26.40000000 28.20000000
#> [103] 30.00000000 32.50000000 35.00000000 37.50000000 40.00000000 43.33333333
#> [109] 46.66666667 50.00000000 53.33333333 56.66666667 60.00000000 64.00000000
#> [115] 68.00000000 72.00000000 76.00000000 80.00000000
# open-ended top group from 60 years
default_age_lower(max_age = 60)
#>   [1]  0.00000000  0.04166667  0.08333333  0.12500000  0.16666667  0.20833333
#>   [7]  0.25000000  0.29166667  0.33333333  0.37500000  0.41666667  0.45833333
#>  [13]  0.50000000  0.54166667  0.58333333  0.62500000  0.66666667  0.70833333
#>  [19]  0.75000000  0.79166667  0.83333333  0.87500000  0.91666667  0.95833333
#>  [25]  1.00000000  1.07142857  1.14285714  1.21428571  1.28571429  1.35714286
#>  [31]  1.42857143  1.50000000  1.57142857  1.64285714  1.71428571  1.78571429
#>  [37]  1.85714286  1.92857143  2.00000000  2.10000000  2.20000000  2.30000000
#>  [43]  2.40000000  2.50000000  2.60000000  2.70000000  2.80000000  2.90000000
#>  [49]  3.00000000  3.14285714  3.28571429  3.42857143  3.57142857  3.71428571
#>  [55]  3.85714286  4.00000000  4.14285714  4.28571429  4.42857143  4.57142857
#>  [61]  4.71428571  4.85714286  5.00000000  5.22222222  5.44444444  5.66666667
#>  [67]  5.88888889  6.11111111  6.33333333  6.55555556  6.77777778  7.00000000
#>  [73]  7.33333333  7.66666667  8.00000000  8.33333333  8.66666667  9.00000000
#>  [79]  9.33333333  9.66666667 10.00000000 10.45454545 10.90909091 11.36363636
#>  [85] 11.81818182 12.27272727 12.72727273 13.18181818 13.63636364 14.09090909
#>  [91] 14.54545455 15.00000000 15.62500000 16.25000000 16.87500000 17.50000000
#>  [97] 18.12500000 18.75000000 19.37500000 20.00000000 20.50000000 21.00000000
#> [103] 22.80000000 24.60000000 26.40000000 28.20000000 30.00000000 32.00000000
#> [109] 34.00000000 36.00000000 38.00000000 40.00000000 43.33333333 46.66666667
#> [115] 50.00000000 53.33333333 56.66666667 60.00000000
# under half the resolution everywhere, same shape, over twice as fast
length(default_age_lower(n_group = 53))
#> [1] 53
```
