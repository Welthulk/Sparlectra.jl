#!/bin/sh
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

# file: tools/gen_powsybl_fixtures.sh
# purpose: regenerate the five PowSyBl fixture bundles under
#          test/fixtures/powsybl from the pypowsybl example networks.
#          Needs a Python with pypowsybl and pandas; activate it BEFORE
#          running this script, for example:
#            python3 -m venv ~/.venv-powsybl
#            ~/.venv-powsybl/bin/pip install pypowsybl pandas
#            PYTHON=~/.venv-powsybl/bin/python sh tools/gen_powsybl_fixtures.sh
#          One dump per line, no command chaining. Run from the repository
#          root.

PYTHON=${PYTHON:-python3}
"$PYTHON" tools/powsybl_dump.py builtin:ieee14 test/fixtures/powsybl/ieee14.powsybl --case ieee14
"$PYTHON" tools/powsybl_dump.py builtin:ieee57 test/fixtures/powsybl/ieee57.powsybl --case ieee57
"$PYTHON" tools/powsybl_dump.py builtin:four_substations_node_breaker_network test/fixtures/powsybl/four_substations.powsybl --case four_substations
"$PYTHON" tools/powsybl_dump.py builtin:micro_grid_be_network test/fixtures/powsybl/micro_grid_be.powsybl --case micro_grid_be
"$PYTHON" tools/powsybl_dump.py builtin:eurostag_tutorial_example1_with_tie_lines_and_areas test/fixtures/powsybl/eurostag_tie_lines.powsybl --case eurostag_tie_lines
