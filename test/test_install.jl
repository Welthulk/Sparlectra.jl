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
            # ext/ holds the package extensions; since 0.30.1 the application
            # loads KLU, so its extension must be in the checkout
            for item in ("Project.toml", "Manifest.toml", "start_webui.jl", "src", "ext", "tools", joinpath("app", "Project.toml"), joinpath("app", "src"))
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
    @testset "sysimage build environment takes a checkout with another version" begin (function ()
        # The build environment is shared by all checkouts of a user. It
        # holds both packages developed from checkout A (0.40.1, application
        # compat Sparlectra 0.40); the build from checkout B (0.42.0) must
        # replace both. Developing them one at a time failed with a
        # ResolverError on the library. The two checkouts are stand-ins with
        # the real names and UUIDs and only a Project.toml and an empty
        # module each, so the child resolves offline in a second or two.
        # The child runs the real develop step of tools/build_sysimage.jl
        # (dry-run include, then `_develop_checkout`) in a temporary
        # environment; the shared build environment is never touched.
        mktempdir() do tmp
            fake_checkout(dir, version, compat) = begin
                mkpath(joinpath(dir, "src"))
                mkpath(joinpath(dir, "app", "src"))
                write(joinpath(dir, "src", "Sparlectra.jl"), "module Sparlectra end\n")
                write(joinpath(dir, "app", "src", "SparlectraApp.jl"), "module SparlectraApp end\n")
                write(joinpath(dir, "Project.toml"), """
                    name = "Sparlectra"
                    uuid = "31ce9bba-fd9d-44a1-b005-f5f509afda38"
                    version = "$(version)"
                    """)
                write(joinpath(dir, "app", "Project.toml"), """
                    name = "SparlectraApp"
                    uuid = "2e39659a-4f4c-4326-b34e-e3df865cea74"
                    version = "$(version)"

                    [deps]
                    Sparlectra = "31ce9bba-fd9d-44a1-b005-f5f509afda38"

                    [sources]
                    Sparlectra = {path = ".."}

                    [compat]
                    Sparlectra = "$(compat)"
                    """)
                dir
            end
            a = fake_checkout(joinpath(tmp, "a"), "0.40.1", "0.40")
            b = fake_checkout(joinpath(tmp, "b"), "0.42.0", "0.42")
            envdir = mkpath(joinpath(tmp, "env"))
            script = joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "build_sysimage.jl")
            child = """
                using Pkg
                sandbox = Module()
                redirect_stdout(devnull) do
                    Base.include(sandbox, $(repr(script)))
                end
                Pkg.activate($(repr(envdir)); io = devnull)
                Base.invokelatest(sandbox._develop_checkout, Pkg, $(repr(a)))
                Base.invokelatest(sandbox._develop_checkout, Pkg, $(repr(b)))
                for (_, d) in Pkg.dependencies()
                    println("dep ", d.name, " ", d.version, " ", d.source)
                end
                """
            exe = joinpath(Sys.BINDIR, Base.julia_exename())
            env = copy(ENV)
            env["JULIA_PKG_OFFLINE"] = "true"
            env["JULIA_PKG_PRECOMPILE_AUTO"] = "0"     # never compile the stand-ins into the depot
            env["JULIA_LOAD_PATH"] = "@:@stdlib"
            env["SPARLECTRA_SYSIMAGE_DRY_RUN"] = "1"
            env["XDG_STATE_HOME"] = joinpath(tmp, "state")
            io = IOBuffer()
            cmd = Cmd(`$(exe) --startup-file=no -e $(child)`; env = env, dir = tmp)
            ok = success(pipeline(cmd; stdout = io, stderr = io))
            text = String(take!(io))
            ok || println("      develop child failed:\n", text)
            @test ok
            for name in ("Sparlectra", "SparlectraApp")
                src = name == "Sparlectra" ? b : joinpath(b, "app")
                @test occursin("dep $(name) 0.42.0 $(src)", text)
            end
        end
    end)() end
    @testset "workshop installer: kernel per folder, backups, restore on a failed install" begin (function ()
        # tools/install_workshops.jl, the parts without network and Jupyter.
        # (1) Two folders get two kernels, the same folder the same one.
        # (2) A notebook the installer wrote and nobody touched is replaced
        # without a backup when the release brings a new one; a notebook the
        # user changed is kept with a date stamp, and a second backup does
        # not overwrite the first. (3) A failed install of Sparlectra puts
        # Project.toml and Manifest.toml back (child process, offline, on a
        # stand-in checkout): before, Pkg.rm ran first and a failed add left
        # the folder without Sparlectra.
        script = joinpath(Sparlectra.SPARLECTRA_ROOT, "tools", "install_workshops.jl")
        ws = Module(:WorkshopInstaller)
        Base.include(ws, script)
        call(name, args...; kw...) = Base.invokelatest(Base.invokelatest(getfield, ws, name), args...; kw...)
        mktempdir() do tmp
            k1 = call(:kernel_names, joinpath(tmp, "sparlectra-workshops"))
            k2 = call(:kernel_names, joinpath(tmp, "other", "sparlectra-workshops"))
            @test k1.spec != k2.spec && k1 == call(:kernel_names, joinpath(tmp, "sparlectra-workshops"))
            @test occursin(r"^sparlectra-sparlectra-workshops-[0-9a-f]{8}$", k1.spec)
            @test occursin(r"^Sparlectra workshops \(sparlectra-workshops [0-9a-f]{8}\)$", k1.display)
            src = mkpath(joinpath(tmp, "src"))
            dst = mkpath(joinpath(tmp, "dst"))
            nb(v) = "{\"metadata\": {\"kernelspec\": {\"name\": \"julia\"}}, \"v\": $(v)}"
            files = ["workshop_a.ipynb", "workshop_b.ipynb"]
            foreach(f -> write(joinpath(src, f), nb(1)), files)
            c1 = call(:copy_notebooks, src, dst, files, k1.spec, k1.display; stamp = "s1")
            @test (c1.updated, c1.backups) == (0, String[])
            @test occursin(k1.spec, read(joinpath(dst, "workshop_a.ipynb"), String))
            # a new release; the user ran notebook b
            foreach(f -> write(joinpath(src, f), nb(2)), files)
            write(joinpath(dst, "workshop_b.ipynb"), "ran")
            c2 = call(:copy_notebooks, src, dst, files, k1.spec, k1.display; stamp = "s2")
            @test (c2.updated, c2.backups) == (2, ["workshop_b.ipynb.s2.bak"])
            @test read(joinpath(dst, "workshop_b.ipynb.s2.bak"), String) == "ran"
            # the same release again: nothing changes, nothing is backed up
            c3 = call(:copy_notebooks, src, dst, files, k1.spec, k1.display; stamp = "s3")
            @test (c3.updated, c3.backups) == (0, String[])
            write(joinpath(dst, "workshop_b.ipynb"), "ran again")
            c4 = call(:copy_notebooks, src, dst, files, k1.spec, k1.display; stamp = "s4")
            @test c4.backups == ["workshop_b.ipynb.s4.bak"] && read(joinpath(dst, "workshop_b.ipynb.s2.bak"), String) == "ran"

            fake = joinpath(tmp, "fake")
            mkpath(joinpath(fake, "src"))
            write(joinpath(fake, "src", "Sparlectra.jl"), "module Sparlectra end\n")
            write(joinpath(fake, "Project.toml"), "name = \"Sparlectra\"\nuuid = \"31ce9bba-fd9d-44a1-b005-f5f509afda38\"\nversion = \"0.1.0\"\n")
            envdir = mkpath(joinpath(tmp, "ws"))
            child = """
                using Pkg
                sandbox = Module()
                Base.include(sandbox, $(repr(script)))
                Pkg.activate($(repr(envdir)); io = devnull)
                Pkg.develop(path = $(repr(fake)); io = devnull)
                before = read.(joinpath.($(repr(envdir)), ["Project.toml", "Manifest.toml"]), String)
                o = (dir = $(repr(envdir)), path = $(repr(joinpath(tmp, "missing"))), rev = nothing, url = "")
                failed = try
                    redirect_stdout(devnull) do
                        Base.invokelatest(sandbox.install_sparlectra, o)
                    end
                    false
                catch err
                    # the expected failure: Pkg refuses the missing checkout
                    err isa Pkg.Types.PkgError || rethrow()
                    true
                end
                after = read.(joinpath.($(repr(envdir)), ["Project.toml", "Manifest.toml"]), String)
                println("failed=", failed, " restored=", before == after, " has=", occursin("Sparlectra", after[1]))
                """
            env = copy(ENV)
            env["JULIA_PKG_OFFLINE"] = "true"
            env["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
            env["JULIA_LOAD_PATH"] = "@:@stdlib"
            io = IOBuffer()
            exe = joinpath(Sys.BINDIR, Base.julia_exename())
            ok = success(pipeline(Cmd(`$(exe) --startup-file=no -e $(child)`; env = env, dir = tmp); stdout = io, stderr = io))
            text = String(take!(io))
            ok || println("      workshop installer child failed:\n", text)
            @test ok && occursin("failed=true restored=true has=true", text)
        end
    end)() end
    end
end
