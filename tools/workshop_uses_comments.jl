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

# file: tools/workshop_uses_comments.jl
# purpose: maintain the "## uses:" line of every code cell of the Literate
#          workshops (docs/lit/workshop_*.jl): the names a cell takes from
#          another cell (functions, variables) and where they come from, so
#          a participant who skips a cell knows what to run first. Splits a
#          file into cells the way Literate does, reads definitions and uses
#          with the Julia parser, and rewrites the first line of each code
#          cell. `--check` only reports cells whose line is missing or stale
#          (the workshops test runs it). Usage:
#          julia tools/workshop_uses_comments.jl [--check] docs/lit/workshop_tour.jl ...

const USES_PREFIX = "## uses: "

# the guard line under a uses line: a cell that takes a name from another
# cell than the warm-up stops with a message naming that cell, instead of
# an UndefVarError, when a participant skipped it
_is_guard(l::AbstractString) = occursin(r"^\(?@isdefined\(.*\|\| error\(\"Run ", l)

# Literate: a line starting with "# " or equal to "#" is markdown; "#src"
# lines exist only in the script, "#nb" lines only in the notebook. Returns
# chunks (is_code, first_line_index, lines).
function literate_chunks(lines::Vector{String})
  chunks = Tuple{Bool,Int,Vector{String}}[]
  kind(l) = (startswith(l, "# ") || l == "#") ? :md : startswith(l, "#nb") ? :nb : :code
  for (i, l) in enumerate(lines)
    k = kind(l) === :md ? false : true
    if isempty(chunks) || first(chunks[end]) != k
      push!(chunks, (k, i, String[]))
    end
    push!(chunks[end][3], l)
  end
  return chunks
end

# the code a cell runs in the notebook: without #src lines, #nb lines with
# their marker removed, without the uses line itself
function cell_code(lines::Vector{String})::String
  kept = String[]
  for l in lines
    occursin("#src", l) && continue
    startswith(l, USES_PREFIX) && continue
    _is_guard(l) && continue
    push!(kept, startswith(l, "#nb ") ? l[5:end] : l)
  end
  return join(kept, "\n")
end

function _defs!(out::Set{Symbol}, ex)
  ex isa Expr || return out
  if ex.head === :function || ((ex.head === :(=)) && ex.args[1] isa Expr && ex.args[1].head === :call)
    sig = ex.args[1]
    while sig isa Expr && sig.head in (:where, :(::))
      sig = sig.args[1]
    end
    sig isa Expr && sig.head === :call && sig.args[1] isa Symbol && push!(out, sig.args[1])
  elseif ex.head === :(=)
    lhs = ex.args[1]
    lhs isa Symbol && push!(out, lhs)
    lhs isa Expr && lhs.head === :tuple && foreach(a -> a isa Symbol && push!(out, a), lhs.args)
    _defs!(out, ex.args[2])
  elseif ex.head in (:block, :toplevel, :macrocall, :const, :global)
    foreach(a -> _defs!(out, a), ex.args)
  end
  return out
end

# every symbol an expression reads; the name of a keyword argument
# (`f(net = pnet)`) is not a read of `net`, only its value is
function _syms!(out::Set{Symbol}, ex)
  if ex isa Symbol
    push!(out, ex)
  elseif ex isa Expr
    if ex.head === :kw
      _syms!(out, ex.args[2])
    elseif _is_function_def(ex)
      # a function body reads its parameters and its own locals, not the
      # globals of the same name (a builder's local `net` is not the
      # chapter's `net`)
      sig, body = ex.args[1], ex.args[2]
      own = _param_names!(Set{Symbol}(), sig)
      _local_assignments!(own, body)
      _signature_reads!(out, sig)
      inner = _syms!(Set{Symbol}(), body)
      union!(out, setdiff(inner, own))
    else
      foreach(a -> _syms!(out, a), ex.args)
    end
  end
  return out
end

_is_function_def(ex::Expr) = ex.head === :function || (ex.head === :(=) && ex.args[1] isa Expr && (ex.args[1].head === :call || (ex.args[1].head in (:where, :(::)) && ex.args[1].args[1] isa Expr && ex.args[1].args[1].head === :call)))

# the names a signature binds: positional, keyword, typed, defaulted
function _param_names!(out::Set{Symbol}, sig)
  while sig isa Expr && sig.head in (:where, :(::))
    sig = sig.args[1]
  end
  sig isa Expr && sig.head === :call || return out
  for a in sig.args[2:end]
    _bind_name!(out, a)
  end
  return out
end

function _bind_name!(out::Set{Symbol}, a)
  if a isa Symbol
    push!(out, a)
  elseif a isa Expr && a.head === :parameters
    foreach(b -> _bind_name!(out, b), a.args)
  elseif a isa Expr && a.head in (:kw, :(::), :(...)) && !isempty(a.args)
    _bind_name!(out, a.args[1])
  end
  return out
end

# what a signature reads: default values (a default may name a global)
function _signature_reads!(out::Set{Symbol}, sig)
  sig isa Expr || return out
  if sig.head === :kw
    _syms!(out, sig.args[2])
  else
    foreach(a -> _signature_reads!(out, a), sig.args)
  end
  return out
end

function _local_assignments!(out::Set{Symbol}, ex)
  ex isa Expr || return out
  if ex.head === :(=) && ex.args[1] isa Symbol
    push!(out, ex.args[1])
  elseif ex.head === :(=) && ex.args[1] isa Expr && ex.args[1].head === :tuple
    foreach(a -> a isa Symbol && push!(out, a), ex.args[1].args)
  elseif ex.head === :for
    _local_assignments!(out, ex.args[1])
  end
  foreach(a -> _local_assignments!(out, a), ex.args)
  return out
end

# where a name comes from, as a reader finds it: the warm-up cell, or the
# nearest "**Example x.y" before the defining cell
function origin_label(lines::Vector{String}, line_idx::Int, warmup_line::Int)::String
  line_idx == warmup_line && return "warm-up cell"
  for i in line_idx:-1:1
    m = match(r"\*\*Example ([0-9]+(\.[0-9]+)*)", lines[i])
    m === nothing || return "Example " * m.captures[1]
    m = match(r"^# #+ (.+)$", lines[i])
    if m !== nothing
      heading = String(strip(m.captures[1]))
      # "Example 7: what current magnitudes can do" is Example 7
      short = match(r"^(Example [0-9]+(\.[0-9]+)*)", heading)
      return short === nothing ? heading : String(short.captures[1])
    end
  end
  return "above"
end

# the names a cell changes without assigning them: the first argument or
# the `net =` keyword of a call to a `!` function (`addVmMeasurement!(net,
# ...)`), a field or an index set (`t.tap_min = 0.9`). A later cell that
# uses the name needs these cells as well, not only the one that built it.
function _mods!(out::Set{Symbol}, ex)
  ex isa Expr || return out
  _is_function_def(ex) && return out
  if ex.head === :call && ex.args[1] isa Symbol && endswith(String(ex.args[1]), "!")
    for a in ex.args[2:end]
      if a isa Symbol
        push!(out, a)
        break
      elseif a isa Expr && a.head === :parameters
        foreach(k -> (k isa Expr && k.head === :kw && k.args[1] === :net && k.args[2] isa Symbol) && push!(out, k.args[2]), a.args)
      elseif a isa Expr && a.head === :kw && a.args[1] === :net && a.args[2] isa Symbol
        push!(out, a.args[2])
      end
    end
  elseif ex.head === :(=) && ex.args[1] isa Expr && ex.args[1].head in (:., :ref) && ex.args[1].args[1] isa Symbol
    push!(out, ex.args[1].args[1])
  end
  foreach(a -> _mods!(out, a), ex.args)
  return out
end

function uses_lines(lines::Vector{String})
  chunks = literate_chunks(lines)
  code_chunks = [(start, ls) for (is_code, start, ls) in chunks if is_code && !all(l -> startswith(l, "#nb") || isempty(strip(l)), ls)]
  warmup = findfirst(c -> occursin("using Sparlectra", cell_code(c[2])), code_chunks)
  warmup_line = warmup === nothing ? 0 : code_chunks[warmup][1]
  # every cell that built or changed a name, in order
  touched = Dict{Symbol,Vector{Int}}()
  result = Dict{Int,Vector{String}}()   # first line of a code chunk => uses line and guards (empty for none)
  for (start, ls) in code_chunks
    ex = Meta.parseall(cell_code(ls))
    mine = _defs!(Set{Symbol}(), ex)
    used = _syms!(Set{Symbol}(), ex)
    changed = _mods!(Set{Symbol}(), ex)
    foreign = sort!([s for s in used if haskey(touched, s) && !(s in mine)])
    # the places a name comes from: the cell that built it and every cell
    # that changed it since; names with the same places share a group
    groups = Dict{Vector{String},Vector{String}}()
    order = Vector{String}[]
    for s in foreign
      labels = unique([origin_label(lines, c, warmup_line) for c in touched[s]])
      haskey(groups, labels) || (push!(order, labels); groups[labels] = String[])
      push!(groups[labels], string(s))
    end
    sort!(order; by = l -> (l != ["warm-up cell"], findfirst(==(l), order)))
    head = isempty(order) ? String[] : [USES_PREFIX * join((string(join(groups[l], ", "), " (", join(l, ", "), ")") for l in order), "; ")]
    guards = String[]
    for l in order
      places = filter(!=("warm-up cell"), l)
      isempty(places) && continue
      ns = groups[l]
      cond = length(ns) == 1 ? "@isdefined($(ns[1]))" : string("(", join(("@isdefined($(n))" for n in ns), " && "), ")")
      where_text = join((startswith(p, "Example") ? p : string("the section \\\"", p, "\\\"") for p in places), ", ")
      push!(guards, string(cond, " || error(\"Run ", where_text, " first: ", length(places) == 1 ? "it sets up " : "they set up ", join(ns, ", "), ".\")"))
    end
    result[start] = vcat(head, guards)
    for s in union(mine, changed)
      (haskey(touched, s) || s in mine) || continue
      start in get!(touched, s, Int[]) || push!(touched[s], start)
    end
  end
  return result
end

# the file with every code cell's uses line set; the line goes first in the
# cell (after leading blank lines), replacing an older one
function rewrite(lines::Vector{String})::Vector{String}
  wanted = uses_lines(lines)
  out = String[]
  starts = Set(keys(wanted))
  i = 1
  # every pass moves i forward by at least one line; a pass count beyond
  # the line count is a bug of this function, stopped here instead of
  # filling memory (a stuck i once grew the output to 26 GB)
  passes = 0
  while i <= length(lines)
    passes += 1
    passes > length(lines) + 1 && error("workshop_uses_comments: rewrite made no progress at line $(i) ($(repr(lines[i]))); stopped after $(passes - 1) passes")
    if i in starts
      # drop an existing uses line at the top of the cell
      j = i
      block = String[]
      while j <= length(lines) && isempty(strip(lines[j]))
        push!(block, lines[j]); j += 1
      end
      j <= length(lines) && startswith(lines[j], USES_PREFIX) && (j += 1)
      while j <= length(lines) && _is_guard(lines[j])
        j += 1
      end
      append!(out, block)
      append!(out, wanted[i])
      # a cell right under a markdown line (no blank line, no old uses
      # line) leaves j == i: copy that line here, or the loop never moves
      if j == i
        push!(out, lines[i])
        j += 1
      end
      i = j
      continue
    end
    push!(out, lines[i])
    i += 1
  end
  return out
end

function main(args)
  check = "--check" in args
  files = filter(a -> a != "--check", args)
  stale = String[]
  for f in files
    lines = readlines(f)
    new = rewrite(lines)
    new == lines && continue
    check ? push!(stale, f) : write(f, join(new, "\n") * "\n")
    check || println("updated ", f)
  end
  if check && !isempty(stale)
    println("stale uses lines in: ", join(stale, ", "))
    exit(1)
  end
  return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && Base.invokelatest(main, ARGS)
