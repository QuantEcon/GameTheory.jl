using Documenter, Games

include("auto_doc_gen.jl")

makedocs(
    modules = [Games],
    format = :html,
    sitename = "Games.jl",
    pages = PAGES,
)
