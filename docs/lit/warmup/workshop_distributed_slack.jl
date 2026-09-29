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

# file: docs/lit/warmup/workshop_distributed_slack.jl
# purpose: compile warm-up of the distributed-slack workshop
#          (docs/lit/workshop_distributed_slack.jl): runs every solver path
#          the chapters use once on the study case, so the chapters do not
#          stall on first-call compilation. Included by the first code cell
#          of the notebook, which defines the helpers used here (load_case,
#          print_beyond_schedule, print_participation). Not part of the
#          library.

"""
    warmup()

Run every compile-heavy path of the distributed-slack workshop once and
print how long the first solves and the rest took. The output of the
paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  t_first = @elapsed runpf!(load_case(), 30, 1e-8, 0)
  t_dslack = @elapsed runpf!(load_case(), 30, 1e-8, 0; distributed_slack_enabled = true)
  println("warm: classical ", round(t_first; digits = 2), " s, distributed ", round(t_dslack; digits = 2), " s (first calls compile)")

  ## every further path of the chapters once, output discarded: the weight
  ## modes with their keyword forms and both print helpers
  t_paths = @elapsed redirect_stdout(devnull) do
    wn = load_case()
    runpf!(wn, 30, 1e-8, 0)
    print_beyond_schedule(wn)
    for mode in (:pg_weighted, :imported, :headroom_weighted)
      wn = load_case()
      runpf!(wn, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = mode)
      print_participation(wn)
      print_beyond_schedule(wn)
    end
    wn = load_case()
    runpf!(wn, 30, 1e-8, 0; distributed_slack_enabled = true, distributed_slack_p_mode = :explicit, distributed_slack_weights = Dict("2" => 1.0))
    print_participation(wn)
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (weight modes, participation tables); everything warm")
  return nothing
end
