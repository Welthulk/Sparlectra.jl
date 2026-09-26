#!/usr/bin/env python3
# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# file: tools/powsybl_dump.py
# purpose: produce a PowSyBl table bundle for Sparlectra from an IIDM file
#          or a pypowsybl example network: pypowsybl loads the network,
#          resolves the node-breaker topology and the tap positions, runs
#          OpenLoadFlow (tight tolerances, see below) for the reference, and
#          this script writes every table with all attributes as CSV plus
#          manifest.json. Sparlectra reads the bundle without Python.
#          Dependencies: pypowsybl and pandas only, Python 3.10 or later.
#          No subprocess is started.
#
# usage:   python tools/powsybl_dump.py <input> <outdir> [--case NAME]
#          <input> is an IIDM file (.xiidm, .xiidm.bz2, .xml) or
#          builtin:<factory>, where <factory> names a
#          pypowsybl.network.create_<factory> function (builtin:ieee14,
#          builtin:four_substations_node_breaker_network).
#
# CSV conventions (the Sparlectra reader relies on them): comma separated,
# one header line, UTF-8, LF; strings always double-quoted with embedded
# quotes doubled; floats in full precision with NaN, Inf, -Inf as literals;
# booleans true/false; integers plain, a missing integer (pandas NaN in an
# integer column) written as -1; index columns first.

import argparse
import bz2
import datetime
import json
import math
import os
import re
import sys

import pandas as pd
import pypowsybl as pp

TABLE_NAMES = [
    "substations", "voltage_levels", "buses", "bus_breaker_view_buses",
    "switches", "busbar_sections", "lines", "2_windings_transformers",
    "3_windings_transformers", "ratio_tap_changers", "ratio_tap_changer_steps",
    "phase_tap_changers", "phase_tap_changer_steps", "generators",
    "reactive_capability_curve_points", "loads", "shunt_compensators",
    "linear_shunt_compensator_sections", "dangling_lines", "tie_lines",
    "hvdc_lines", "vsc_converter_stations", "lcc_converter_stations",
    "static_var_compensators", "operational_limits",
]

# Columns the Sparlectra schema types as Int. pandas turns an integer column
# into floats as soon as it holds a NaN; these are written as integers
# again, a NaN becomes -1.
INT_COLUMNS = {
    "connected_component", "synchronous_component", "node", "node1", "node2", "node3",
    "tap", "low_tap", "high_tap", "step_count", "position", "num",
    "max_section_count", "section_count", "acceptable_duration",
    "ratio_tap_position1", "ratio_tap_position2", "ratio_tap_position3",
    "phase_tap_position1", "phase_tap_position2", "phase_tap_position3",
}


def csv_string(value):
    return '"' + str(value).replace('"', '""') + '"'


def csv_cell(value, column):
    if value is None:
        return '""'
    if isinstance(value, (bool,)) or type(value).__name__ == "bool_":
        return "true" if bool(value) else "false"
    if column in INT_COLUMNS:
        if isinstance(value, float) and math.isnan(value):
            return "-1"
        try:
            return str(int(value))
        except (TypeError, ValueError):
            return csv_string(value)
    if isinstance(value, (int,)) and not isinstance(value, bool):
        return str(value)
    if isinstance(value, float) or type(value).__name__ in ("float64", "float32"):
        f = float(value)
        if math.isnan(f):
            return "NaN"
        if math.isinf(f):
            return "Inf" if f > 0 else "-Inf"
        return repr(f)
    if hasattr(value, "item"):
        return csv_cell(value.item(), column)
    if isinstance(value, float) is False and pd.isna(value):
        return '""'
    return csv_string(value)


def write_table(path, df, index_names):
    df = df.reset_index()
    # index columns first, then the value columns in pypowsybl order
    columns = [c for c in index_names if c in df.columns] + [c for c in df.columns if c not in index_names]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(",".join(csv_string(c) for c in columns) + "\n")
        for row in df[columns].itertuples(index=False, name=None):
            fh.write(",".join(csv_cell(v, c) for v, c in zip(row, columns)) + "\n")
    return len(df)


def get_table(network, name):
    getter = getattr(network, "get_" + name)
    try:
        return getter(all_attributes=True)
    except TypeError:
        return getter()


def iidm_version_of(path):
    if path.endswith(".bz2"):
        with bz2.open(path, "rt", encoding="utf-8") as fh:
            head = fh.read(4096)
    else:
        with open(path, "r", encoding="utf-8") as fh:
            head = fh.read(4096)
    m = re.search(r"schema/iidm/([0-9_]+)", head)
    return m.group(1) if m else ""


def main(argv):
    parser = argparse.ArgumentParser(description="Write a PowSyBl table bundle for Sparlectra.")
    parser.add_argument("input")
    parser.add_argument("outdir")
    parser.add_argument("--case", default=None)
    args = parser.parse_args(argv)

    if args.input.startswith("builtin:"):
        factory = args.input[len("builtin:"):]
        creator = getattr(pp.network, "create_" + factory, None)
        if creator is None:
            print("powsybl_dump: no pypowsybl factory create_" + factory, file=sys.stderr)
            return 2
        network = creator()
        case = args.case or factory
        source = "pypowsybl.network.create_" + factory
        source_file = case + ".xiidm"
        write_xiidm = True
    else:
        network = pp.network.load(args.input)
        stem = os.path.basename(args.input)
        for ext in (".xiidm.bz2", ".xiidm", ".xml"):
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        case = args.case or stem
        source = args.input
        source_file = os.path.basename(args.input)
        write_xiidm = False

    os.makedirs(args.outdir, exist_ok=True)

    iidm_version = ""
    if write_xiidm:
        # the factory network is saved unsolved and read back: the bundle
        # then carries the tables in the row order a load of that .xiidm
        # produces, so a later dump of the file (this script or the
        # PythonCall extension) reproduces the bundle cell by cell; a
        # factory network keeps its own element order, which differs
        xiidm_path = os.path.join(args.outdir, source_file)
        network.save(xiidm_path, format="XIIDM")
        iidm_version = iidm_version_of(xiidm_path)
        network = pp.network.load(xiidm_path)
    else:
        try:
            iidm_version = iidm_version_of(args.input)
        except (OSError, UnicodeDecodeError):
            iidm_version = ""

    # the reference solution first: the tables then carry the OpenLoadFlow
    # state (bus voltages, branch and injection flows) that the tests use.
    # OpenLoadFlow's defaults leave up to 1 MW undistributed at the slack
    # bus (slackBusPMaxMismatch) and stop Newton-Raphson at 1e-4 pu per
    # equation; a reference judged at 1e-6 pu needs both tightened. Every
    # other parameter keeps its default (distributed slack proportional to
    # max_p, remote voltage control on, reactive limits on).
    parameters = pp.loadflow.Parameters(provider_parameters={"slackBusPMaxMismatch": "0.0001", "newtonRaphsonConvEpsPerEq": "1e-9"})
    results = pp.loadflow.run_ac(network, parameters)
    components = []
    for r in results:
        components.append({
            "num": int(r.connected_component_num),
            "synchronous_component_num": int(r.synchronous_component_num),
            "status": str(r.status).split(".")[-1],
            "iterations": int(r.iteration_count),
            "reference_bus_id": str(r.reference_bus_id),
        })

    tables = {}
    for name in TABLE_NAMES:
        try:
            df = get_table(network, name)
        except Exception as err:  # the task asks for the table name on failure
            print("powsybl_dump: table " + name + " failed: " + repr(err), file=sys.stderr)
            return 1
        index_names = [n for n in df.index.names if n is not None]
        file_name = name + ".csv"
        rows = write_table(os.path.join(args.outdir, file_name), df, index_names)
        tables[name] = {"file": file_name, "rows": rows, "index": index_names}

    # reference_buses.csv: one row per bus-breaker bus, joined with the
    # bus-view bus and the voltage level after the load flow
    bb = network.get_bus_breaker_view_buses().reset_index()
    vls = network.get_voltage_levels()
    ref_rows = []
    for row in bb.itertuples(index=False):
        nominal_v = float(vls.loc[row.voltage_level_id, "nominal_v"])
        v_mag = float(row.v_mag)
        ref_rows.append((row.bus_id, row.id, v_mag, float(row.v_angle), v_mag / nominal_v if nominal_v > 0 else float("nan"), row.synchronous_component, nominal_v))
    ref_columns = ["id", "bus_breaker_id", "v_mag_kv", "v_angle_deg", "v_pu", "synchronous_component", "nominal_v"]
    with open(os.path.join(args.outdir, "reference_buses.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(",".join(csv_string(c) for c in ref_columns) + "\n")
        for row in ref_rows:
            fh.write(",".join(csv_cell(v, c) for v, c in zip(row, ref_columns)) + "\n")

    manifest = {
        "format": "powsybl_tables",
        "format_version": 1,
        "case": case,
        "source": source,
        "source_file": source_file,
        "pypowsybl_version": pp.__version__,
        "iidm_version": iidm_version,
        "all_attributes": True,
        "written": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tables": tables,
        "reference": {
            "file": "reference_buses.csv",
            "loadflow": "OpenLoadFlow via pypowsybl.loadflow.run_ac, default parameters except slackBusPMaxMismatch=0.0001 and newtonRaphsonConvEpsPerEq=1e-9",
            "components": components,
        },
    }
    with open(os.path.join(args.outdir, "manifest.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print("powsybl_dump: wrote " + args.outdir + " (" + case + ", " + str(sum(t["rows"] for t in tables.values())) + " rows in " + str(len(tables)) + " tables)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
