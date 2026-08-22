# Sparlectra workshop notebooks

These notebooks run on Google Colab, no local installation required. The
first cell installs Sparlectra from GitHub (branch `main`) into a fresh
temporary environment (takes a few minutes); the isolation keeps Colab's
preinstalled package stack out of the install, and commented lines in the
same cell switch to another branch or to the latest registered release.

| Notebook | What it covers | Open |
|---|---|---|
| [workshop_tour.ipynb](workshop_tour.ipynb) | The full workshop in one session, five tiers from Newcomer to Beyond: first network step by step, model editing and Q-limits, slack types and short circuit, OLTC tap control, Q(U) control, remote voltage control, a steerable HVDC link, state estimation, parallel sweeps on Julia threads | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_tour.ipynb) |
| [workshop_slack_short_circuit.ipynb](workshop_slack_short_circuit.ipynb) | Three ways to model the grid connection (ideal slack, non-ideal external-grid source, distributed slack) plus an IEC 60909-0 short-circuit sweep from the declared feeder data | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_slack_short_circuit.ipynb) |
| [workshop_distributed_slack.ipynb](workshop_distributed_slack.ipynb) | Distributed slack: how the participation weights are determined (schedule, capacity, headroom, imported APF/normalPF, explicit), how they are normalized, and the fallback when no valid participant exists | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_distributed_slack.ipynb) |
| [workshop_transformers.ipynb](workshop_transformers.ipynb) | Transformer taps: ratio taps (OLTC) move voltage, phase taps (PST/Schrägregler) move power; device tap math, the closed control loop with the X(α) characteristic, and the 3WT star equivalent | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_transformers.ipynb) |
| [workshop_series_compensation.ipynb](workshop_series_compensation.ipynb) | Series compensation (TCSC): why flow follows reactance, and how the series-reactance controller steers a loop-network flow split onto a branch target, including the at_limit case | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_series_compensation.ipynb) |
| [workshop_state_estimation.ipynb](workshop_state_estimation.ipynb) | State estimation: derive a noisy measurement set from a reference power flow, check observability, run the WLS estimator, and see when observability breaks down | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_state_estimation.ipynb) |

The same content is rendered as documentation pages in the
[Notebooks section](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour/)
of the hosted docs.

**Do not edit the `.ipynb` files directly**: they are generated. Edit the
Literate.jl source in [`docs/lit/`](../docs/lit/) and regenerate with
`julia --project=docs docs/generate_notebooks.jl`; the CI Documentation
workflow also regenerates them on every push to `main`.
