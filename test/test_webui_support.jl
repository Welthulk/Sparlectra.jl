# Copyright 2023–2026 Udo Schmitz
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

_dtf_network_fixture() = "P\nT\nN\nX\nY\n110\n2 1 0 0 SLACK\nL1A  PV       SLACK   1.0 2.0 0.0 0.0\n11   PV      110 0 0 0 10 2 -5 5\n21   SLACK   110 0 0 0 20 3 -10 10\n"
_dtf_network_with_outage_fixture() = _dtf_network_fixture() * "AUSFALL\nL1A\nENDE\n"
_dtf_outage_fixture() = "AUSFALL\nL1A\nENDE\n"
_dtf_reference_fixture() = "FOR002 reference result report\nBUS VOLTAGE RESULTS\n"
