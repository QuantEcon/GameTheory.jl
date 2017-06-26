using Documenter, Games

include("auto_doc_gen.jl")

makedocs(
    modules = [Games],
    format = :html,
    sitename = "Games.jl",
    pages = PAGES,
)

deploydocs(
    repo = "github.com/QuantEcon/Games.jl.git",
    branch = "gh-pages",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = nothing,
    make = nothing,
)
