using Documenter, GameTheory

include("auto_doc_gen.jl")

makedocs(
    modules = [GameTheory],
    format = Documenter.HTML(prettyurls = false),
    sitename = "GameTheory.jl",
    pages = PAGES,
)

deploydocs(
    repo = "github.com/QuantEcon/GameTheory.jl.git",
    branch = "gh-pages",
    target = "build",
    deps = nothing,
    make = nothing,
)
