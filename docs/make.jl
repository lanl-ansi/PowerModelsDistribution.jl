using Documenter, ThreePhasePowerModels

makedocs(
    modules = [ThreePhasePowerModels],
    format = :html,
    sitename = "ThreePhasePowerModels",
    authors = "Carleton Coffrin, David Fobes, and contributors.",
    analytics = "",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickguide.md",
        ],
        "Library" => "library.md",
        "Developer" => "developer.md",
    ]
)

deploydocs(
    deps = nothing,
    make = nothing,
    target = "build",
    repo = "github.com/lanl-ansi/ThreePhasePowerModels.jl.git",
    julia = "0.6"
)
