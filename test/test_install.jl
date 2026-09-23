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
# file: test/test_install.jl
# purpose: the installation path of a fresh checkout: start_webui.jl sets up
#          the library and the application environments and compiles both on
#          the first start, and does nothing on the second. Its own profile
#          (install): a cold compile of both packages in a copy of the
#          checkout costs about a minute and belongs before a release, not
#          into every Web UI run.

using Sparlectra
using Test

function run_install_tests()
    @testset "installation path" begin
    @testset "fresh checkout: the first start sets up both environments, the second does nothing" begin (function ()
        # a checkout without any Manifest.toml, as git clone or a release
        # download leaves it: the first start resolves, installs and
        # compiles the library and the application and says so, the
        # second start finds both up to date and compiles nothing. The
        # copy carries the library manifest so the resolution stays
        # offline and deterministic; the application manifest is what a
        # fresh checkout lacks. Two child processes, about three minutes
        # on the development machine (a cold compile of both packages).
        mktempdir() do tmp
            root = Sparlectra.SPARLECTRA_ROOT
            copy_dir = joinpath(tmp, "checkout")
            for item in ("Project.toml", "Manifest.toml", "start_webui.jl", "src", "tools", joinpath("app", "Project.toml"), joinpath("app", "src"))
                src = joinpath(root, item)
                dst = joinpath(copy_dir, item)
                mkpath(dirname(dst))
                cp(src, dst)
            end
            @test !isfile(joinpath(copy_dir, "app", "Manifest.toml"))
            exe = joinpath(Sys.BINDIR, Base.julia_exename())
            env = copy(ENV)
            env["JULIA_PKG_OFFLINE"] = "true"          # resolve against the depot, no registry traffic
            env["SPARLECTRA_PRECOMPILE_WORKLOAD"] = "off"
            start(label) = begin
                io = IOBuffer()
                cmd = Cmd(`$(exe) --startup-file=no --project=$(copy_dir) $(joinpath(copy_dir, "start_webui.jl")) --env-only`; env = env, dir = copy_dir)
                ok = success(pipeline(cmd; stdout = io, stderr = io))
                text = String(take!(io))
                ok || println("      ", label, " start failed:\n", text)
                @test ok
                text
            end
            first_start = start("first")
            @test occursin("First start: setting up the application environment", first_start)
            @test occursin("The application environment is ready after", first_start)
            @test occursin("Environment check finished.", first_start)
            @test isfile(joinpath(copy_dir, "app", "Manifest.toml"))
            second_start = start("second")
            @test occursin("Environments up to date: library and application.", second_start)
            @test occursin("Packages already compiled.", second_start)
            @test occursin("Environment check finished.", second_start)
            @test !occursin("First start", second_start)
        end
    end)() end
    end
end
