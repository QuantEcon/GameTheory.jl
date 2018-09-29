#=
Generate markdown files automatically, for generating
docs using `Documenter.jl`.

@author: Zejin Shi

=#

path = replace(pathof(Games), "src/Games.jl" => "")

# read the basic structures
order = SubString{String}[]
sections_names = SubString{String}[]
sections = Dict{SubString{String}, Vector{SubString{String}}}()
re = r"(.*): (.*)\n"
open(joinpath(path, "docs/Structure")) do f
    for match in eachmatch(re, String(read(f)))
        list = map(i -> strip(i), split(match.captures[2], ","))
        if list ==[""]
            continue
        end

        if match.captures[1] == "Order"
            order = list
            global order
        else
            section_name = match.captures[1]
            push!(sections_names, section_name)
            sections[section_name] = list
        end
    end
end

modules = String["Games"]
# find all files used in Games.jl
re = r"include\(\"(.*)\.jl\"\)"
files = String[]

open(joinpath(path, "src/Games.jl")) do f
    for match in eachmatch(re, String(read(f)))
        # check if it is a submodule
        if occursin("/", match.captures[1])
            submodule_path = match.captures[1]
            submodule_dir, submodule_name = split(submodule_path, "/")
            # add submodule name to modules
            push!(modules, "Games.$submodule_name")
            # find all files used in submodules
            open(joinpath(path, "src/$submodule_path.jl")) do f_submodule
                for match_sub in eachmatch(re, String(read(f_submodule)))
                    file = match_sub.captures[1]
                    push!(files, "$submodule_dir/$file")
                end
            end
        else
            push!(files, match.captures[1])
        end
    end
end

module_names = join(modules, ", ")

files_names = Dict{String, String}(
                file => replace(join(map(uppercasefirst, split(file, "_")), " "),
                                "Util" => "Utilities")
                for file in files
                )

# generate paths
if !ispath(joinpath(path, "docs/src/lib"))
    mkpath(joinpath(path, "docs/src/lib"))
end

# write .md for files in main modules not in section
for file in files
    if file in vcat(values(sections)...)
        continue
    end

    if !(file in order)
        push!(order, file)
    end
    file_name = files_names[file]
    file_page = """
# [$file_name](@id $file)

This is documentation for `$file.jl`.

## Exported

```@autodocs
Modules = [$module_names]
Pages   = ["$file.jl"]
Private = false
```

## Internal

```@autodocs
Modules = [$module_names]
Pages   = ["$file.jl"]
Public = false
```
"""
    open(joinpath(path, "docs/src/lib/$file.md"), "w") do f
        write(f, file_page)
    end
end

# write .md for sections
for section_name in sections_names
    if !(section_name in order)
        push!(order, section_name)
    end

    section_files = sections[section_name]
    section_name_lower = replace(lowercase(section_name), " " => "_")
    section_file_list = join(map(i -> string("\"", i, ".jl\""),
                             section_files), ", ")
    # include "Games.jl" in "Base Types and Methods"
    if section_name == "Base Types and Methods"
        section_file_list = string("\"Games.jl\", ", section_file_list)
    end
    section_page = """
# [$section_name](@id $section_name_lower)

## Exported
```@autodocs
Modules = [$module_names]
Pages   = [$section_file_list]
Private = false
```

## Internal
```@autodocs
Modules = [$module_names]
Pages   = [$section_file_list]
Public = false
```
"""
    open(joinpath(path, "docs/src/lib/$section_name_lower.md"), "w") do f
        write(f, section_page)
    end
end

# write index page for Liabrary
index = """
# Index

```@index
Modules = [$module_names]
```
"""

open(joinpath(path, "docs/src/lib/index.md"), "w") do f
    write(f, index)
end

# add Library Outline to Homepage index.md
open(joinpath(path, "docs/src/index.md"), "a") do f
    for page in order
        if page in keys(sections)
            section_name_lower = replace(lowercase(page), " " => "_")
            write(f, "* [$page](@ref $section_name_lower)\n\n")
        else
            file = page
            file_name = files_names[file]
            write(f, "* [$file_name](@ref $file)\n\n")
        end
    end
end

pages_names = map(i -> replace(lowercase(i), " " => "_"), order)

PAGES = ["Home" => "index.md",
         "Library" => push!(
                        ["lib/$page_name.md" 
                         for page_name in pages_names
                        ],
                         "lib/index.md"
                    )
    ]
