using Documenter, PowerModelsDistribution

makedocs(
    modules = [PowerModelsDistribution],
    format = Documenter.HTML(analytics = "", mathengine = Documenter.MathJax()),
    sitename = "PowerModelsDistribution",
    authors = "David Fobes, Carleton Coffrin, Frederik Geth, Sander Claeys, and contributors.",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "quickguide.md",
            "Connecting Components" => "connections.md",
            "Mathematical Model" => "math-model.md",
            "Engineering Data Model" => "eng-data-model.md",
            "Enums in Engineering Model" => "enums.md",
            "Conversion to Mathematical Model" => "eng2math.md",
            "External Data Formats" => "external-data-formats.md",
            "Examples" => "engineering_model.md",
        ],
        "Library" => [
            "Network Formulations" => "formulations.md",
            "Problem Specifications" => "specifications.md",
            "Modeling Components" => "library.md",
        ],
        "Basics" => [
            "Load Models" => "load-model.md",
        ],
        "Developer" => [
            "Developer" => "developer.md",
            "Formulation Details" => "formulation-details.md",
        ],
    ]
)

deploydocs(
    repo = "github.com/lanl-ansi/PowerModelsDistribution.jl.git",
)
