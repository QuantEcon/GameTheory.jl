using Documenter, GameTheory

makedocs(
    modules = [GameTheory],
    format = Documenter.HTML(prettyurls = false),
    sitename = "GameTheory.jl",
    pages = [
        "Home" => "index.md",
        "Library" => [
            "lib/base_types_and_methods.md",
            "lib/game_generators.md", 
            "lib/computing_nash_equilibria.md",
            "lib/learning_algorithms.md",
            "lib/repeated_games.md",
            "lib/util.md",
            "lib/index.md"
        ]
    ],
)

deploydocs(
    repo = "github.com/QuantEcon/GameTheory.jl.git",
    branch = "gh-pages",
    target = "build",
    deps = nothing,
    make = nothing,
)
