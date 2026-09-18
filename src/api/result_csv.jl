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

# file: src/api/result_csv.jl
# purpose: the one writer for result CSV artifacts (issue #386): every run
#          type keeps its own columns and its own function, the delimiter,
#          decimal separator, thousands separator and quoting come from
#          output.csv_format and from nowhere else. Data formats that
#          Sparlectra reads back with a fixed layout (the measurement CSV
#          "# sparlectra-measurements v1") are NOT result artifacts and do
#          not use this writer.

"""
    result_csv_format() -> String

The CSV format name of result artifacts, `output.csv_format` of the active
configuration (`technical`, `excel_de`, `excel_us`). Writers that are called
without an explicit `format` take it from here, so a run writes every
artifact in one format.
"""
result_csv_format()::String = String(output_config().csv_format)

"""
    write_result_csv(path, header, rows; format = result_csv_format(), comments = String[]) -> path
    write_result_csv(io::IO, header, rows; format = result_csv_format(), comments = String[])

Write one result CSV: optional `# comment` lines, the header, one line per
row. `rows` is any iterable of iterables (tuples, vectors, generators);
every cell goes through the shared formatting (`_csv_field`): numbers get
the decimal and thousands separators of `format`, strings are quoted when
they contain the delimiter, a quote or a line break, `nothing`/`missing`
become empty cells. `format` is a format name (`technical`, `excel_de`,
`excel_us`) or an already resolved format; the default is the central
setting `output.csv_format` of the active configuration.
"""
function write_result_csv(io::IO, header, rows; format = result_csv_format(), comments = String[])
  resolved = format isa NamedTuple ? format : _resolve_detailed_csv_format(format)
  delimiter = resolved.delimiter
  for c in comments
    println(io, "# ", c)
  end
  println(io, join((String(h) for h in header), delimiter))
  for row in rows
    println(io, join((_csv_field(v, delimiter, resolved) for v in row), delimiter))
  end
  return nothing
end

function write_result_csv(path::AbstractString, header, rows; format = result_csv_format(), comments = String[])
  open(path, "w") do io
    write_result_csv(io, header, rows; format = format, comments = comments)
  end
  return path
end

"""
    result_csv_delimiter(line) -> Char

The delimiter a result CSV header line was written with (`;` or `,`), for
readers that accept both formats of `output.csv_format`.
"""
result_csv_delimiter(line::AbstractString)::Char = occursin(';', line) ? ';' : ','

# parse a numeric cell written by `write_result_csv` in any of the three
# formats: thousands separators are stripped, a decimal comma becomes a dot
function _parse_result_csv_number(text::AbstractString, delimiter::Char)
  t = strip(text)
  if delimiter == ';'
    t = replace(replace(t, "." => ""), "," => ".")
  else
    t = replace(t, "," => "")
  end
  return tryparse(Float64, t)
end
