## Generate docs automatically

To generate documentation for `Games.jl`, run

```
$julia make.jl
```

The generated docs can be accessed in `docs/build/index.html`.

There are three parts to generate the documentation.

### auto_doc_gen.jl

Main part of documentation generation. It reads `src/Games.jl` to find
out all the source files, and then write markdown files for each with
the structure instructed by `docs/build/Structure`, which will be utilized
by "Documenter.jl" later.

Note that the pattern of each page can be modified by changing the
strings in `auto_doc_gen.jl` called `file_page` and `section_page`.

### manually written files

Two files should be written manually:

#### structure

You may modify `Structure` to arrange the structure of the website to be generated.
For example:

```
Order: Base Types and Methods, Computing Nash Equilibria, Repeated Games
Base Types and Methods: normal_form_game, random
Computing Nash Equilibria: pure_nash
Repeated Games: repeated_game_util, repeated_game
```

The first line sets the order of pages. You can add a file name if you want it to
generate a single page, or the name of the section you want to create. In the following
you need to declare which files are included in each section. Note that you should not
put a comma or dot at the end of each line.

It is not necessary to declare the structure of every file included in the Module.
A separate single page will be generated for each file not mentioned in `Structure` automatically.

#### Homepage

The Homepage is written in `/docs/src/index.md`. Please keep the last part in Homepage,
called `Library Outline`, empty. It will be filled by `auto_doc_gen.jl` automatically.

### make.jl

This file includes information for "Documenter.jl" to generate documentation.
Details on how to use it can be found [here](https://juliadocs.github.io/Documenter.jl/stable/man/guide.html#Output-formats-1).
