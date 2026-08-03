# blink: a mean-field (ODE) twin of malariasimulation

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

- [`run_simulation_ode()`](https://pwinskill.github.io/blink/reference/run_simulation_ode.md)
  — run the model; returns a wide, daily count table (P. falciparum
  only).

- [`get_epi_outputs()`](https://pwinskill.github.io/blink/reference/get_epi_outputs.md)
  — post-process a run into postie `rates` (long) and `prevalence`
  (wide).

- [`default_age_lower()`](https://pwinskill.github.io/blink/reference/default_age_lower.md)
  — the default graded age grid.

See
[`vignette("blink")`](https://pwinskill.github.io/blink/articles/blink.md)
for a worked tour, and the README's *Mean-field approximations* section
for where and how much the ODE departs from the IBM.

## See also

Useful links:

- <https://pwinskill.github.io/blink>

- <https://github.com/pwinskill/blink>

- Report bugs at <https://github.com/pwinskill/blink/issues>

## Author

**Maintainer**: Peter Winskill <p.winskill@imperial.ac.uk>
