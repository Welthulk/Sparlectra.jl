# SPDX-License-Identifier: MIT

"""
    checker.jl

Utility script to scan all Julia source files in the `src` directory
and identify function parameters that are declared but never used in the function body.

This is a code-quality diagnostic tool that helps maintainers find dead parameters
that may indicate incomplete refactoring, API changes, or unintended overloads.

Unused parameters starting with underscore (`_`) are excluded by design, as they
follow Julia convention for intentionally unused arguments.

# Execution

```julia
julia --project=. checker.jl
```

# Output

Reports a list of function definitions with their unused parameters, grouped by file path.
If no unused parameters are found, reports success.
"""

root = joinpath(@__DIR__, "src")

if !isdir(root)
  error("Directory not found: $root")
end

# Remove :where clauses to extract the core function signature
unwrap_where(sig) = (sig isa Expr && sig.head == :where) ? unwrap_where(sig.args[1]) : sig

"""
    arg_symbol(a)

Extract the parameter name from a function argument expression.

Handles type annotations (`::`), defaults (`=`), keyword arguments (`kw`),
and variadic arguments (`...`).

Returns a Symbol for named parameters, or nothing if the argument has no extractable name.
"""
function arg_symbol(a)
  a isa Symbol && return a
  a isa Expr || return nothing
  if a.head == :(::) || a.head == :(=) || a.head == :kw || a.head == :(...)
    return arg_symbol(a.args[1])
  end
  return nothing
end

"""
    func_name(sig)

Extract the function name from a function signature expression.

Handles plain function names, method calls (for parameterized or type-constraint signatures),
and returns `"<anonymous>"` for lambda or other unnamed expressions.
"""
function func_name(sig)
  s = unwrap_where(sig)
  if s isa Symbol
    return String(s)
  elseif s isa Expr && s.head == :call
    return string(s.args[1])
  else
    return "<anonymous>"
  end
end

"""
    collect_params(sig)

Collect all named parameter symbols from a function signature.

Extracts positional arguments, keyword arguments, and variadic parameters.
Excludes parameters starting with underscore (Julia convention for intentionally unused args).

Returns a Vector of Symbol.
"""
function collect_params(sig)
  s = unwrap_where(sig)
  s isa Expr && s.head == :call || return Symbol[]
  out = Symbol[]
  for a in s.args[2:end]
    if a isa Expr && a.head == :parameters
      # Keyword arguments are stored in a :parameters Expr
      for kw in a.args
        sym = arg_symbol(kw)
        sym isa Symbol && !startswith(String(sym), "_") && push!(out, sym)
      end
    else
      # Positional arguments
      sym = arg_symbol(a)
      sym isa Symbol && !startswith(String(sym), "_") && push!(out, sym)
    end
  end
  return out
end

"""
    uses_symbol(ex, sym::Symbol)

Recursively check whether a symbol is used (referenced) in an expression.

Traverses the expression tree deeply; returns true if the symbol is found anywhere.
"""
function uses_symbol(ex, sym::Symbol)
  ex === sym && return true
  ex isa Expr || return false
  for a in ex.args
    uses_symbol(a, sym) && return true
  end
  return false
end

"""
    scan_expr!(hits, ex, file)

Recursively scan an expression for function definitions with unused parameters.

When a `:function` expression is found, collects all named parameters and checks
which ones do not appear in the function body. Adds a NamedTuple to `hits` if
any unused parameters are found.

Continues recursively through nested expressions to find all function definitions.
"""
function scan_expr!(hits, ex, file)
  if ex isa Expr && ex.head == :function
    sig, body = ex.args[1], ex.args[2]
    params = collect_params(sig)
    unused = [p for p in params if !uses_symbol(body, p)]
    !isempty(unused) && push!(hits, (file = file, fn = func_name(sig), unused = unused))
  end
  if ex isa Expr
    for a in ex.args
      scan_expr!(hits, a, file)
    end
  end
end

begin
  println("Starting unused-parameter scan in: ", root)

  hits = NamedTuple[]
  local nfiles = 0

  # Walk the source tree and scan each .jl file
  for (dir, _, files) in walkdir(root)
    for f in files
      endswith(f, ".jl") || continue
      nfiles += 1
      path = joinpath(dir, f)
      println("  Scanning: ", path)

      # Parse the file incrementally to handle parse errors gracefully
      src = read(path, String)
      pos = firstindex(src)
      while true
        ex, newpos = Meta.parse(src, pos; greedy = false, raise = false, filename = path)
        ex === nothing && break

        # Guard against parser stalls (e.g. parse errors with unchanged position)
        if newpos <= pos
          pos == lastindex(src) && break
          pos = nextind(src, pos)
          continue
        end

        # Ignore parse-error nodes, keep scanning
        if !(ex isa Expr && ex.head == :error)
          scan_expr!(hits, ex, path)
        end
        pos = newpos
      end
    end
  end

  println("\nScan complete. Processed $nfiles .jl files.\n")

  # Report results
  if nfiles == 0
    println("No .jl files found under: ", root)
  elseif isempty(hits)
    println("No unused-parameter candidates found.")
  else
    println("Unused-parameter candidates:")
    for h in hits
      println("- ", h.file, " :: ", h.fn, " -> ", join(string.(h.unused), ", "))
    end
  end
end