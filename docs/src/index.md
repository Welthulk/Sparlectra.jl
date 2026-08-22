# Sparlectra

```@meta
EditURL = "../../README.md"
```

```@eval
using Markdown

readme_path = normpath(joinpath(@__DIR__, "..", "..", "README.md"))
readme_text = read(readme_path, String)

# Keep docs-internal links working after reusing README content. The README
# addresses images and pages relative to the repository root; on the built
# site the same targets live relative to this page. Markdown links AND raw
# HTML img/anchor attributes must be rewritten (the screenshots and the logo
# are plain <img> tags, invisible to the Markdown-link regex).
readme_text = replace(readme_text, r"\(docs/src/([^)]+)\)" => s"(\1)")
readme_text = replace(readme_text, r"src=\"docs/src/([^\"]+)\"" => s"src=\"\1\"")
readme_text = replace(readme_text, r"href=\"docs/src/([^\"]+)\"" => s"href=\"\1\"")

Markdown.parse(readme_text)
```

## Documentation quick links

* [Feature Matrix](feature_matrix.md)
* [Central Configuration](configuration.md)
* [N-1 Contingency Analysis](contingency.md)
* [Q-limit Switching Strategy](q_limit_switching_strategy.md)
* [Synthetic Tiled Grids](synthetic_grids.md)
* [Branch Model](branchmodel.md)
* [Examples Overview](examples_overview.md)
