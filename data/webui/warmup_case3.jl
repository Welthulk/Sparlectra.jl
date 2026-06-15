# Original Sparlectra-owned synthetic three-bus MATPOWER-compatible warm-up case.
# This small case is intentionally independent of external MATPOWER case data.
(
  name = "sparlectra_webui_warmup_case3",
  baseMVA = 100.0,
  bus = [
    1.0 3.0  0.0  0.0 0.0 0.0 1.0 1.02 0.0 110.0 1.0 1.10 0.90;
    2.0 2.0 20.0 10.0 0.0 0.0 1.0 1.01 0.0 110.0 1.0 1.10 0.90;
    3.0 1.0 45.0 15.0 0.0 0.0 1.0 1.00 0.0 110.0 1.0 1.10 0.90;
  ],
  gen = [
    1.0 50.0  0.0 100.0 -100.0 1.02 100.0 1.0 150.0 0.0;
    2.0 25.0  0.0  80.0  -80.0 1.01 100.0 1.0 100.0 0.0;
  ],
  branch = [
    1.0 2.0 0.02 0.08 0.02 999.0 999.0 999.0 0.0 0.0 1.0 -360.0 360.0;
    1.0 3.0 0.03 0.12 0.02 999.0 999.0 999.0 0.0 0.0 1.0 -360.0 360.0;
    2.0 3.0 0.02 0.10 0.02 999.0 999.0 999.0 0.0 0.0 1.0 -360.0 360.0;
  ],
  gencost = nothing,
  bus_name = ["Warmup Slack", "Warmup PV", "Warmup Load"],
)
