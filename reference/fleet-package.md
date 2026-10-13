# fleet: mean-field counterpart to malariasimulation

A deterministic mean-field counterpart to the malariasimulation
individual-based model (IBM) of *Plasmodium falciparum* and *P. vivax*
malaria, advanced one day at a time on the IBM's daily clock and built
on odin2 and dust2. It reproduces the human model, structured by age and
biting heterogeneity and, for *P. vivax*, by hypnozoite batch, coupled
to the compartmental mosquito model, with the same intervention modules.
It takes an unmodified
[`malariasimulation::get_parameters()`](https://rdrr.io/pkg/malariasimulation/man/get_parameters.html)
list, is seeded at the malariaEquilibrium (or malariaEquilibriumVivax)
solution, and returns `malariasimulation`'s output table, which postie
post-processes as it would an IBM run. Run time does not depend on
population size.

## Functions

- [`run_simulation_ode()`](https://pwinskill.github.io/fleet/reference/run_simulation_ode.md)
  runs the model.

- [`ode_tuning()`](https://pwinskill.github.io/fleet/reference/ode_tuning.md)
  holds the discretisation settings, passed as `tuning`.

- [`default_age_lower()`](https://pwinskill.github.io/fleet/reference/default_age_lower.md)
  gives the default age grid.

The `ode` in two of the names is kept for compatibility: the model is a
daily update, not an ODE.

[`vignette("fleet")`](https://pwinskill.github.io/fleet/articles/fleet.md)
is a worked introduction, and
[`vignette("using")`](https://pwinskill.github.io/fleet/articles/using.md)
covers where the mean field departs from the IBM and what to do about
it. Agreement with the IBM is measured in a separate project,
[fleetcheck](https://pwinskill.github.io/fleetcheck/).

## See also

Useful links:

- <https://pwinskill.github.io/fleet/>

- <https://github.com/pwinskill/fleet>

- Report bugs at <https://github.com/pwinskill/fleet/issues>

## Author

**Maintainer**: Peter Winskill <p.winskill@imperial.ac.uk> \[copyright
holder\]
