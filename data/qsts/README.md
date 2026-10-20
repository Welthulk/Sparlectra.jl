# Shipped QSTS example

## sp_case14_day.qsts.zip (a day on sp_case14)

An SCF time series: `sp_case14.scf.json` unchanged as base case and 96
steps of 15 minutes on 2026-06-15 (UTC). The three 110 kV loads follow a
commercial shape, the five 20 kV loads a household shape, each with 2
percent AR(1) noise of its own; the two machines other than the reference
follow the load, the reference balances. `Ostheim_Neuwiese_1` is out in
steps 57 to 64 (14:00 to 15:45). The taps of `Moorau_110_20` and
`Wehrden_110_20` and the shunt at `Moorau_20` stay as in the base case.
Built by `tools/build_scf_series.jl` (seed 20261013; without `--write` it
checks this file). `sp_case14_day.control.yaml` holds the controllers for
runs with control: a tap controller for `Wehrden_110_20`, the control
range `tap_min: 1.0` for the Moorau tap controller of the case and the
shunt at `Moorau_20` as a switched bank behind it (target 1.01 pu). See
the QSTS page of the documentation.
