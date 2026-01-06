# GameTheory.jl Documentation

This directory contains the documentation for GameTheory.jl, built using [Documenter.jl](https://documenter.juliadocs.org/).

## Documentation Structure

The documentation uses a static structure with manually maintained pages:

- `src/index.md` - Main homepage with installation, usage examples, and library outline
- `src/lib/` - Library documentation pages organized by topic

## Local Build Instructions

To build the documentation locally from the repository root:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

This will build the docs with the version that `docs/Manifest.toml` currently resolves
(e.g. a released version of GameTheory.jl).

If you are editing the source code and want the docs to reflect your local checkout:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'
```

To revert this (go back to the registry-resolved version in the `docs` environment):

```bash
julia --project=docs -e 'using Pkg; Pkg.free("GameTheory"); Pkg.instantiate()'
```

The generated documentation will be available in `docs/build/index.html`.

## Contributing to Documentation

### Adding New Content

To add documentation for new functionality:

1. **For new library functions/types**: Add docstrings to the source code and ensure they appear in the appropriate `docs/src/lib/*.md` file by including the source file in the `@autodocs` block.

2. **For new major features**: Create a new page in `docs/src/lib/` and add it to the `pages` array in `make.jl`.

3. **For examples or tutorials**: Add to the main `docs/src/index.md` or create dedicated pages as needed.

### Editing Existing Pages

- **Homepage**: Edit `docs/src/index.md` for installation instructions, basic usage, and overview
- **Library pages**: Edit files in `docs/src/lib/` to modify section organization or add manual content
- **Structure**: Modify the `pages` array in `make.jl` to change page ordering or add new sections

### Documentation Standards

- Use `@autodocs` blocks to automatically include docstrings from source files
- Maintain consistent section structure across library pages (Exported/Internal)
- Include practical examples in docstrings where helpful
- Keep the library outline in `index.md` synchronized with actual pages
