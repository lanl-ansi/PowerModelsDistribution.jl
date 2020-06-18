using Documenter, PowerModelsDistribution
import Weave

Weave.set_chunk_defaults!(Dict{Symbol, Any}(:line_width => 120))

examples = []
cd("examples")
for file in readdir(".")
    if endswith(file, "ipynb")
        doc_filepath = Weave.weave("$file"; out_path="../docs/src", fig_path="../docs/src/assets", doctype="github", mod=Main)
        doc_file = basename(doc_filepath)

        push!(examples, doc_file)
    end
end
cd("..")

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
            "Examples" => [titlecase(replace(file[1:end-3], "_" => " ")) => file for file in examples],
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
