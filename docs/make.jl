using Documenter, Games

include("auto_doc_gen.jl")

makedocs(
    modules = [Games],
    format = Documenter.HTML(prettyurls = false),
    sitename = "Games.jl",
    pages = PAGES,
)

deploydocs(
    repo = "github.com/QuantEcon/Games.jl.git",
    branch = "gh-pages",
    target = "build",
    deps = nothing,
    make = nothing,
)
