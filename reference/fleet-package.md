# fleet: a mean-field (ODE) twin of malariasimulation

A fast, deterministic mean-field (ODE) counterpart to the
malariasimulation individual-based model of *Plasmodium falciparum*
malaria, built on the odin2 / dust2 stack. It reproduces the age- and
biting-heterogeneity-structured human model (states S/D/A/U/Tr plus
treatment and chemoprevention prophylaxis, and the six immunity
functions) coupled to the compartmental mosquito model, and honours the
same intervention modules. It accepts an unmodified
[`malariasimulation::get_parameters()`](https://rdrr.io/pkg/malariasimulation/man/get_parameters.html)
list, is seeded at the malariaEquilibrium fixed point, and returns
malariasimulation-style outputs that feed straight into postie. A
multi-decade run completes in seconds, independent of population size.

## Getting started

- [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md):
  run the model; returns a wide, daily count table (P. falciparum only).

- [`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md):
  the default graded age grid.

- [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md):
  solver and discretisation settings, passed as `tuning =`.

The count table is malariasimulation-shaped, so post-process it with
postie exactly as you would an IBM run: `postie::get_rates(out)` and
`postie::get_prevalence(out, diagnostic = "lm")`. See
[`vignette("fleet")`](https://pwinskill.github.io/fleet/articles/fleet.md)
for a worked tour,
[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
for where and how much the ODE departs from the IBM and what to do about
it, and
[`vignette("comparison")`](https://pwinskill.github.io/fleet/articles/comparison.md)
for the measured agreement behind those claims.

## See also

Useful links:

- <https://pwinskill.github.io/fleet>

- <https://github.com/pwinskill/fleet>

- Report bugs at <https://github.com/pwinskill/fleet/issues>

## Author

**Maintainer**: Peter Winskill <p.winskill@imperial.ac.uk> \[copyright
holder\]
