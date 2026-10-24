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
#
# file: app/src/webui/help_math.jl
# purpose: formulas of the help texts as readable HTML. The documentation
#          writes math as `$...$` and ```` ```math ```` blocks, which
#          Documenter renders on the published site; the Web UI has no
#          math renderer, so the formulas are converted here: Greek letters
#          and relations as Unicode, indices as <sub>/<sup>, \sqrt and
#          \frac as text. Only the commands listed below are converted; any
#          other command stays visible as `\name`, and the excerpt check of
#          tools/generate_webui_help_excerpts.jl fails on it.

# Commands that stand for one character or word.
const _WEBUI_TEX_SYMBOLS = Dict{String,String}(
  # Greek, lowercase (TeX shapes: \phi is the closed form, \varphi the open one)
  "alpha" => "α", "beta" => "β", "gamma" => "γ", "delta" => "δ", "epsilon" => "ϵ", "varepsilon" => "ε",
  "zeta" => "ζ", "eta" => "η", "theta" => "θ", "vartheta" => "ϑ", "iota" => "ι", "kappa" => "κ",
  "lambda" => "λ", "mu" => "μ", "nu" => "ν", "xi" => "ξ", "pi" => "π", "varpi" => "ϖ", "rho" => "ρ",
  "varrho" => "ϱ", "sigma" => "σ", "varsigma" => "ς", "tau" => "τ", "upsilon" => "υ", "phi" => "ϕ",
  "varphi" => "φ", "chi" => "χ", "psi" => "ψ", "omega" => "ω",
  # Greek, uppercase
  "Gamma" => "Γ", "Delta" => "Δ", "Theta" => "Θ", "Lambda" => "Λ", "Xi" => "Ξ", "Pi" => "Π",
  "Sigma" => "Σ", "Upsilon" => "Υ", "Phi" => "Φ", "Psi" => "Ψ", "Omega" => "Ω",
  # relations and operators
  "le" => "≤", "leq" => "≤", "ge" => "≥", "geq" => "≥", "ne" => "≠", "neq" => "≠", "approx" => "≈",
  "sim" => "∼", "simeq" => "≃", "equiv" => "≡", "propto" => "∝", "ll" => "≪", "gg" => "≫",
  "cdot" => "·", "times" => "×", "div" => "÷", "pm" => "±", "mp" => "∓", "ast" => "∗", "star" => "⋆",
  "circ" => "∘", "bullet" => "•", "in" => "∈", "notin" => "∉", "subset" => "⊂", "subseteq" => "⊆",
  "cup" => "∪", "cap" => "∩", "setminus" => "∖", "wedge" => "∧", "vee" => "∨", "neg" => "¬",
  "forall" => "∀", "exists" => "∃", "emptyset" => "∅", "infty" => "∞", "partial" => "∂", "nabla" => "∇",
  "sum" => "∑", "prod" => "∏", "int" => "∫", "to" => "→", "rightarrow" => "→", "leftarrow" => "←",
  "Rightarrow" => "⇒", "Leftarrow" => "⇐", "Leftrightarrow" => "⇔", "leftrightarrow" => "↔",
  "mapsto" => "↦", "angle" => "∠", "perp" => "⊥", "parallel" => "∥", "top" => "⊤", "dagger" => "†",
  "prime" => "′", "degree" => "°", "ell" => "ℓ", "hbar" => "ℏ", "Re" => "ℜ", "Im" => "ℑ",
  # delimiters and dots
  "vert" => "|", "lvert" => "|", "rvert" => "|", "mid" => "|", "Vert" => "‖", "lVert" => "‖",
  "rVert" => "‖", "langle" => "⟨", "rangle" => "⟩", "lfloor" => "⌊", "rfloor" => "⌋", "lceil" => "⌈",
  "rceil" => "⌉", "lbrace" => "{", "rbrace" => "}", "ldots" => "…", "dots" => "…", "cdots" => "⋯",
  "vdots" => "⋮", "colon" => ":",
  # operator names, written upright as words
  "max" => "max", "min" => "min", "log" => "log", "ln" => "ln", "exp" => "exp", "sin" => "sin",
  "cos" => "cos", "tan" => "tan", "arcsin" => "arcsin", "arccos" => "arccos", "arctan" => "arctan",
  "sinh" => "sinh", "cosh" => "cosh", "tanh" => "tanh", "arg" => "arg", "det" => "det", "dim" => "dim",
  "lim" => "lim", "sup" => "sup", "inf" => "inf", "deg" => "deg", "ker" => "ker", "bmod" => "mod",
  # spacing
  "quad" => " ", "qquad" => " ",
)

# Commands written as a backslash and one non-letter character.
const _WEBUI_TEX_ESCAPES = Dict{Char,String}(',' => " ", ':' => " ", ';' => " ", ' ' => " ", '!' => "", '{' => "{", '}' => "}", '%' => "%", '_' => "_", '#' => "#", '&' => "&", '|' => "‖")

# Accents: the combining character set after a single letter.
const _WEBUI_TEX_ACCENTS = Dict{String,Char}("hat" => '̂', "bar" => '̄', "overline" => '̅', "tilde" => '̃', "dot" => '̇', "ddot" => '̈', "vec" => '⃗', "check" => '̌')

# Size and pairing commands that only scale a delimiter: dropped.
const _WEBUI_TEX_DROPPED = Set(["left", "right", "big", "Big", "bigg", "Bigg", "bigl", "bigr", "Bigl", "Bigr", "biggl", "biggr", "displaystyle", "textstyle"])

# Greek letters: a following space before a letter or digit is dropped,
# as TeX does ("\Delta x" reads as "Δx").
const _WEBUI_TEX_LETTERLIKE = Set(["alpha", "beta", "gamma", "delta", "epsilon", "varepsilon", "zeta", "eta", "theta", "vartheta", "iota", "kappa", "lambda", "mu", "nu", "xi", "pi", "varpi", "rho", "varrho", "sigma", "varsigma", "tau", "upsilon", "phi", "varphi", "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Upsilon", "Phi", "Psi", "Omega", "partial", "nabla", "ell", "hbar"])

# One symbol, optionally indexed: needs no parentheses under a root or in
# a fraction ("√x", "a/b", but "√(a + b)").
const _WEBUI_TEX_ATOM = r"^\s*(?:[A-Za-z]|[0-9]+(?:\.[0-9]+)?|\\[A-Za-z]+)\s*(?:[_^]\s*(?:\{[^{}]*\}|\\[A-Za-z]+|[A-Za-z0-9])\s*)*$"

mutable struct _WebUITexReader
  chars::Vector{Char}
  pos::Int
end

_webui_tex_eof(r::_WebUITexReader)::Bool = r.pos > length(r.chars)
_webui_tex_peek(r::_WebUITexReader)::Char = r.chars[r.pos]

function _webui_tex_skip_space!(r::_WebUITexReader)
  while !_webui_tex_eof(r) && isspace(_webui_tex_peek(r))
    r.pos += 1
  end
end

# The text up to the matching closing brace; the reader stands after the
# opening one and ends after the closing one.
function _webui_tex_group_source!(r::_WebUITexReader)::String
  depth = 1
  start = r.pos
  while !_webui_tex_eof(r)
    c = _webui_tex_peek(r)
    if c == '\\' && r.pos < length(r.chars)
      r.pos += 2
      continue
    end
    c == '{' && (depth += 1)
    c == '}' && (depth -= 1)
    r.pos += 1
    depth == 0 && return String(r.chars[start:(r.pos - 2)])
  end
  return String(r.chars[start:end])
end

# A command name after the backslash: a run of letters, or one other character.
function _webui_tex_command_name!(r::_WebUITexReader)::String
  _webui_tex_eof(r) && return ""
  start = r.pos
  if isletter(_webui_tex_peek(r)) && isascii(_webui_tex_peek(r))
    while !_webui_tex_eof(r) && isascii(_webui_tex_peek(r)) && isletter(_webui_tex_peek(r))
      r.pos += 1
    end
  else
    r.pos += 1
  end
  return String(r.chars[start:(r.pos - 1)])
end

# One argument of a command or an index: a brace group, a command or one
# character. Returns the HTML and the TeX source of the argument.
function _webui_tex_argument!(r::_WebUITexReader)
  _webui_tex_skip_space!(r)
  _webui_tex_eof(r) && return ("", "")
  c = _webui_tex_peek(r)
  if c == '{'
    r.pos += 1
    source = _webui_tex_group_source!(r)
    return (_webui_tex_html(source), source)
  elseif c == '\\'
    start = r.pos
    r.pos += 1
    html = _webui_tex_command!(r, _webui_tex_command_name!(r))
    return (html, String(r.chars[start:(r.pos - 1)]))
  end
  r.pos += 1
  return (_webui_escape(string(c)), string(c))
end

_webui_tex_wrap(html::AbstractString, source::AbstractString)::String = occursin(_WEBUI_TEX_ATOM, source) ? String(html) : string("(", html, ")")

# The HTML of one command; the reader stands after its name.
function _webui_tex_command!(r::_WebUITexReader, name::AbstractString)::String
  if length(name) == 1 && !isletter(name[1])
    return _webui_escape(get(_WEBUI_TEX_ESCAPES, name[1], "\\" * name))
  elseif haskey(_WEBUI_TEX_SYMBOLS, name)
    if name in _WEBUI_TEX_LETTERLIKE
      save = r.pos
      _webui_tex_skip_space!(r)
      (!_webui_tex_eof(r) && (isletter(_webui_tex_peek(r)) || isdigit(_webui_tex_peek(r)))) || (r.pos = save)
    end
    return _WEBUI_TEX_SYMBOLS[name]
  elseif name in _WEBUI_TEX_DROPPED
    # \left. and \right. are invisible delimiters
    _webui_tex_skip_space!(r)
    !_webui_tex_eof(r) && _webui_tex_peek(r) == '.' && (r.pos += 1)
    return ""
  elseif name == "sqrt"
    _webui_tex_skip_space!(r)
    degree = ""
    if !_webui_tex_eof(r) && _webui_tex_peek(r) == '['
      start = r.pos + 1
      while !_webui_tex_eof(r) && _webui_tex_peek(r) != ']'
        r.pos += 1
      end
      degree = "<sup>" * _webui_tex_html(String(r.chars[start:(r.pos - 1)])) * "</sup>"
      r.pos += 1
    end
    html, source = _webui_tex_argument!(r)
    return string(degree, "√", _webui_tex_wrap(html, source))
  elseif name in ("frac", "dfrac", "tfrac")
    num, num_source = _webui_tex_argument!(r)
    den, den_source = _webui_tex_argument!(r)
    return string(_webui_tex_wrap(num, num_source), "/", _webui_tex_wrap(den, den_source))
  elseif name in ("text", "textrm", "mbox", "textit")
    # text mode: the content as it stands
    _webui_tex_skip_space!(r)
    if !_webui_tex_eof(r) && _webui_tex_peek(r) == '{'
      r.pos += 1
      return _webui_escape(_webui_tex_group_source!(r))
    end
    return first(_webui_tex_argument!(r))
  elseif name in ("mathrm", "operatorname", "mathit", "mathsf", "mathcal")
    return first(_webui_tex_argument!(r))
  elseif name in ("mathbf", "boldsymbol", "bm")
    return string("<b>", first(_webui_tex_argument!(r)), "</b>")
  elseif name == "pmod"
    return string(" (mod ", first(_webui_tex_argument!(r)), ")")
  elseif haskey(_WEBUI_TEX_ACCENTS, name)
    save = r.pos
    html, source = _webui_tex_argument!(r)
    # a combining accent sits on one character only; anything longer stays
    # unconverted, so the excerpt check names it
    length(html) == 1 && return string(html, _WEBUI_TEX_ACCENTS[name])
    r.pos = save
  end
  # not covered: the command stays as written, with its brace argument,
  # so the output reads `\begin{bmatrix}` and the check names `\begin`
  if !_webui_tex_eof(r) && _webui_tex_peek(r) == '{'
    r.pos += 1
    return _webui_escape(string("\\", name, "{", _webui_tex_group_source!(r), "}"))
  end
  return _webui_escape("\\" * name)
end

"""
    _webui_tex_html(formula) -> String

The TeX source of a formula (without the `\$` delimiters) as HTML text:
Greek letters and relations as Unicode, `_`/`^` as `<sub>`/`<sup>`,
`\\sqrt{x}` as `√x` (with parentheses when the argument is more than one
symbol), `\\frac{a}{b}` as `a/b`, `\\text`/`\\mathrm` content as plain text.
Characters are HTML-escaped. A command not covered stays as `\\name`
(escaped), so a check of the output finds it; nothing is dropped silently.
"""
function _webui_tex_html(formula::AbstractString)::String
  r = _WebUITexReader(collect(String(formula)), 1)
  io = IOBuffer()
  while !_webui_tex_eof(r)
    c = _webui_tex_peek(r)
    r.pos += 1
    if c == '\\'
      print(io, _webui_tex_command!(r, _webui_tex_command_name!(r)))
    elseif c == '_' || c == '^'
      tag = c == '_' ? "sub" : "sup"
      print(io, "<", tag, ">", first(_webui_tex_argument!(r)), "</", tag, ">")
    elseif c == '{'
      print(io, _webui_tex_html(_webui_tex_group_source!(r)))
    elseif c == '\''
      print(io, "′")
    elseif isspace(c)
      _webui_tex_skip_space!(r)
      print(io, ' ')
    else
      print(io, _webui_escape(string(c)))
    end
  end
  return String(take!(io))
end

"""
    _webui_tex_plain(formula) -> String

The formula as plain text for the search index: the HTML of
`_webui_tex_html` with `<sub>` written as `_` and `<sup>` as `^`, the
remaining tags removed and the entities resolved.
"""
function _webui_tex_plain(formula::AbstractString)::String
  text = replace(_webui_tex_html(formula), "<sub>" => "_", "<sup>" => "^")
  text = replace(text, r"<[^>]+>" => "")
  return replace(text, "&lt;" => "<", "&gt;" => ">", "&quot;" => "\"", "&#39;" => "'", "&amp;" => "&")
end

# A formula node in a parsed Markdown tree: written as converted HTML by the
# Markdown standard library's HTML output (show for text/html).
struct _WebUIMathNode
  html::String
  display::Bool
end

Base.show(io::IO, ::MIME"text/html", m::_WebUIMathNode) = print(io, m.display ? "<p class=\"math\">" : "<span class=\"math\">", m.html, m.display ? "</p>" : "</span>")

# Replace every Markdown.LaTeX node of a parsed Markdown tree by its
# converted HTML. A formula inside a paragraph, heading, emphasis, link or
# table cell is inline, one at block level (```` ```math ````) a display.
function _webui_math_nodes!(x, inline::Bool)
  x isa Markdown.LaTeX && return _WebUIMathNode(_webui_tex_html(x.formula), !inline)
  if x isa Vector{Any}
    for i in eachindex(x)
      x[i] = _webui_math_nodes!(x[i], inline)
    end
    return x
  elseif x isa AbstractVector
    foreach(item -> _webui_math_nodes!(item, inline), x)
    return x
  end
  parentmodule(typeof(x)) === Markdown || return x
  child_inline = x isa Union{Markdown.Paragraph,Markdown.Header,Markdown.Bold,Markdown.Italic,Markdown.Link,Markdown.Table}
  for field in (:content, :text, :items, :rows)
    hasfield(typeof(x), field) || continue
    value = getfield(x, field)
    value isa AbstractVector && _webui_math_nodes!(value, child_inline)
  end
  return x
end
