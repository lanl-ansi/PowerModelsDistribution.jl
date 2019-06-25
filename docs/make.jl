using Documenter, ThreePhasePowerModels

makedocs(
    modules = [ThreePhasePowerModels],
    format = Documenter.HTML(analytics = ""),
    sitename = "ThreePhasePowerModels",
    authors = "Carleton Coffrin, David Fobes, Frederik Geth and contributors.",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickguide.md",
            "Mathematical Model" => "math-model.md",
        ],
        "Library" => [
            "Network Formulations" => "formulations.md",
            "Problem Specifications" => "specifications.md",
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
