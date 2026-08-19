# Sparlectra workshop notebooks

These notebooks run on Google Colab, no local installation required. The
first cell installs Sparlectra into a fresh temporary environment (takes a
few minutes); the isolation keeps Colab's preinstalled package stack out of
the install, and a commented line in the same cell shows how to test an
unreleased development version instead of the latest release.

| Notebook | What it covers | Open |
|---|---|---|
| [workshop_tour.ipynb](workshop_tour.ipynb) | The whole workshop in one session: install once, warm up once, then power flow, slack types and short circuit, OLTC tap control, Q(U) control, remote voltage control, a steerable HVDC link, and state estimation | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb) |
| [workshop_intro.ipynb](workshop_intro.ipynb) | First contact: build a 7-bus network from scratch, validate it, solve the power flow with Newton-Raphson, read the classical result tables | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_intro.ipynb) |
| [workshop_slack_short_circuit.ipynb](workshop_slack_short_circuit.ipynb) | Three ways to model the grid connection (ideal slack, non-ideal external-grid source, distributed slack) plus an IEC 60909-0 short-circuit sweep from the declared feeder data | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb) |
| [workshop_state_estimation.ipynb](workshop_state_estimation.ipynb) | State estimation: derive a noisy measurement set from a reference power flow, check observability, run the WLS estimator, and see when observability breaks down | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb) |
| [workshop_series_compensation.ipynb](workshop_series_compensation.ipynb) | Series compensation (TCSC): why flow follows reactance, and how the series-reactance controller steers a loop-network flow split onto a branch target, including the honest at_limit case | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb) |

The same content is rendered as documentation pages in the
[Notebooks section](https://welthulk.github.io/Sparlectra.jl/generated/workshop_intro/)
of the hosted docs.

**Do not edit the `.ipynb` files directly**: they are generated. Edit the
Literate.jl source in [`docs/lit/`](../docs/lit/) and regenerate with
`julia --project=docs docs/generate_notebooks.jl`; the CI Documentation
workflow also regenerates them on every push to `main`.
