#=
Generate markdown files automatically, for generating
docs using `Documenter.jl`.

@author: Zejin Shi

=#

# find all files used in Games.jl
re = r"include\(\"(.*)\.jl\"\)"
files = String[]
path = Pkg.dir("Games")

open(joinpath(path, "src/Games.jl")) do f
    for match in eachmatch(re, readstring(f))
        push!(files, match.captures[1])
    end
end

# create PAGES for makedocs()
PAGES = ["Home" => "index.md",
         "Library" => push!(
                        ["lib/$file.md" for file in files],
                         "lib/index.md"
                        )
         ]

# generate paths
if !ispath(joinpath(path, "docs/src/lib"))
    mkpath(joinpath(path, "docs/src/lib"))
end

# write index.md as Homepage
home_page = """
# Games.jl

"""

open(joinpath(path, "docs/src/index.md"), "w") do f
    write(f, home_page)
    for file in files
        file_name = 
            replace(join(map(ucfirst, split(file, "_")), " "),
                    "Util",
                    "Utilities"
                )
        write(f, "* [$file_name](@ref)\n")
    end
end

# write .md for each file in Library
for file in files
    file_name = 
            replace(join(map(ucfirst, split(file, "_")), " "),
                    "Util",
                    "Utilities"
                )

    file_page = """
# $file_name

This is documentation for `$file.jl`.

## Exported

```@autodocs
Modules = [Games]
Pages   = ["$file.jl"]
Private = false
```

## Internal

```@autodocs
Modules = [Games]
Pages   = ["$file.jl"]
Public = false
```
"""
    open(joinpath(path, "docs/src/lib/$file.md"), "w") do f
        write(f, file_page)
    end
end

# write index page for Liabrary
index = """
# Index

```@index
modules = [Games]
```
"""

open(joinpath(path, "docs/src/lib/index.md"), "w") do f
    write(f, index)
end
