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

# file: test/test_webui_support.jl
# purpose: shared fixtures for the Web UI tests: small inline DTF network,
#          outage, and FOR002 reference-report strings, and a form
#          serializer that posts a rendered form the way a browser does
_dtf_network_fixture() = "P\nT\nN\nX\nY\n110\n2 1 0 0 SLACK\nL1A  PV       SLACK   1.0 2.0 0.0 0.0\n11   PV      110 0 0 0 10 2 -5 5\n21   SLACK   110 0 0 0 20 3 -10 10\n"
_dtf_network_with_outage_fixture() = _dtf_network_fixture() * "AUSFALL\nL1A\nENDE\n"
_dtf_outage_fixture() = "AUSFALL\nL1A\nENDE\n"
_dtf_reference_fixture() = "FOR002 reference result report\nBUS VOLTAGE RESULTS\n"

# The entries a browser submits for the form `form_id` in `html`: text,
# number and hidden inputs with their value, checkboxes and radios only when
# checked, a select with its selected (else first) option, disabled controls
# dropped. A repeated name becomes a vector, as the route's form parser
# delivers it. A form without an id is found by a piece of its markup
# (`containing`), the first form that holds it.
function _webui_browser_form(html::AbstractString, form_id::AbstractString; containing::AbstractString = "")
    body = ""
    for m in eachmatch(r"<form\b", html)
        stop = findnext("</form>", html, m.offset)
        candidate = html[m.offset:last(stop)]
        hit = isempty(containing) ? startswith(candidate, "<form id=\"$(form_id)\"") : occursin(containing, candidate)
        hit && (body = candidate; break)
    end
    isempty(body) && error("form $(isempty(containing) ? form_id : containing) not found")
    unescape(s) = replace(String(s), "&quot;" => "\"", "&#39;" => "'", "&lt;" => "<", "&gt;" => ">", "&amp;" => "&")
    attr(tag, a) = (m = match(Regex("\\s$(a)=\"([^\"]*)\""), tag); m === nothing ? nothing : unescape(m.captures[1]))
    flag(tag, a) = occursin(Regex("\\s$(a)(\\s|>|=|/)"), tag)
    form = Dict{String,Any}()
    add!(name, value) = (form[name] = haskey(form, name) ? vcat(form[name], value) : value)
    for m in eachmatch(r"<(input|select|textarea)\b[^>]*>"s, body)
        tag = m.match
        name = attr(tag, "name")
        (name === nothing || flag(tag, "disabled")) && continue
        if m.captures[1] == "input"
            type = something(attr(tag, "type"), "text")
            type in ("submit", "button", "reset", "image", "file") && continue
            type in ("checkbox", "radio") && !flag(tag, "checked") && continue
            add!(name, something(attr(tag, "value"), type in ("checkbox", "radio") ? "on" : ""))
        elseif m.captures[1] == "select"
            options = collect(eachmatch(r"<option\b[^>]*>(?:([^<]*))?"s, body[m.offset:first(findnext("</select>", body, m.offset))]))
            isempty(options) && continue
            chosen = something(findfirst(o -> flag(o.match, "selected"), options), 1)
            add!(name, something(attr(options[chosen].match, "value"), strip(something(options[chosen].captures[1], ""))))
        else
            add!(name, unescape(body[m.offset+length(tag):first(findnext("</textarea>", body, m.offset))-1]))
        end
    end
    return form
end
