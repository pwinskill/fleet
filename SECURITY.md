# Security policy

`blink` is a research modelling package. It runs simulations from a parameter list you
supply; it opens no network connections, reads no credentials, and writes only where
you tell it to. The realistic risk surface is small, but reports are welcome.

## Reporting a vulnerability

Please report privately rather than opening a public issue:

- Use GitHub's [private vulnerability reporting](https://github.com/pwinskill/blink/security/advisories/new), or
- email <p.winskill@imperial.ac.uk>.

Please include the `blink` commit, your `sessionInfo()`, and a reproducible example.
We will acknowledge within a week.

## Scope

In scope:

- Anything in this repository that executes attacker-controlled input, for example a
  parameter list or a file path that leads to arbitrary code execution.
- A GitHub Actions workflow in `.github/workflows/` that could be made to leak a token
  or run untrusted code with elevated permissions.

Out of scope, and better raised as a normal issue:

- **Wrong numbers.** Model output being incorrect is a correctness bug, not a
  vulnerability, and it is the thing we most want to hear about. Open an issue.
- Vulnerabilities in dependencies (`odin2`, `dust2`, `malariasimulation`, and so on).
  Report those to their own maintainers; tell us if `blink` needs a version bump.

## Supported versions

There is no released version yet. Only the current `main` branch is supported. See the
warning at the top of the [README](README.md) before using this package for anything
that matters.
