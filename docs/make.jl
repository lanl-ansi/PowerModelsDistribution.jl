using Documenter, ThreePhasePowerModels

makedocs(
    modules = [ThreePhasePowerModels],
    format = :html,
    sitename = "ThreePhasePowerModels",
    authors = "Carleton Coffrin, David Fobes, Frederik Geth and contributors.",
    analytics = "",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickguide.md",
            "Mathematical Model" => "math-model.md",
        ],
        "Library" => [
            "Network Formulations" => "formulations.md",
            "Modeling Components" => "library.md",
        ],
        "Developer" => [
            "Developer" => "developer.md",
            "Formulation Details" => "formulation-details.md",
        ],
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/ThreePhasePowerModels.jl.git",
)
