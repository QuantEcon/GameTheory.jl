# GameTheory.jl Documentation

This directory contains the documentation for GameTheory.jl, built using [Documenter.jl](https://documenter.juliadocs.org/).

## Documentation Structure

The documentation uses a static structure with manually maintained pages:

- `src/index.md` - Main homepage with installation, usage examples, and library outline
- `src/lib/` - Library documentation pages organized by topic:
  - `base_types_and_methods.md` - Core types and basic operations
  - `game_generators.md` - Tools for generating games
  - `computing_nash_equilibria.md` - Nash equilibrium algorithms  
  - `learning_algorithms.md` - Learning and evolutionary dynamics
  - `repeated_games.md` - Repeated games functionality
  - `util.md` - Utility functions
  - `index.md` - Library index page

## Local Build Instructions

To build the documentation locally:

1. Navigate to the docs directory:
   ```bash
   cd docs/
   ```

2. Install documentation dependencies:
   ```bash
   julia --project=. -e "import Pkg; Pkg.instantiate()"
   ```

3. Build the documentation:
   ```bash
   julia --project=. make.jl
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
