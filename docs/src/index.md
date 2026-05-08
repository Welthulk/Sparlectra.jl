# Sparlectra

```@meta
EditURL = "../../README.md"
```

```@eval
using Markdown

readme_path = normpath(joinpath(@__DIR__, "..", "..", "README.md"))
readme_text = read(readme_path, String)

# Keep docs-internal links working after reusing README content.
readme_text = replace(readme_text, r"\(docs/src/([^)]+)\)" => s"(\1)")

Markdown.parse(readme_text)
```

## Documentation quick links

* [Feature Matrix](feature_matrix.md)
* [Synthetic Tiled Grids](synthetic_grids.md)
* [Transformer Control](transformer_control.md)
* [Examples Overview](examples_overview.md)
