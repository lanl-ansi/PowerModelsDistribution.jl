using Documenter
using PowerModelsDistribution

import Pluto
import Gumbo

const _FAST = findfirst(isequal("--fast"), ARGS) !== nothing

makedocs(
    modules = [PowerModelsDistribution],
    format = Documenter.HTML(
        analytics = "",
        mathengine = Documenter.MathJax(),
        prettyurls=false,
        collapselevel=1,
    ),
    strict=false,
    sitename = "PowerModelsDistribution",
    authors = "David M Fobes, Carleton Coffrin, Frederik Geth, Sander Claeys, and contributors",
    pages = [
        "Introduction" => "index.md",
        "installation.md",
        "Manual" => [
            "Getting Started" => "manual/quickguide.md",
            "External Data Formats" => "manual/external-data-formats.md",
            "Engineering Data Model" => "manual/eng-data-model.md",
            "Enums in Engineering Model" => "manual/enums.md",
            "Mathematical Model" => "manual/math-model.md",
            "Conversion to Mathematical Model" => "manual/eng2math.md",
            "Unbalanced Formulations" => "manual/formulations.md",
            "Problem Specifications" => "manual/specifications.md",
            "Load Models" => "manual/load-model.md",
            "Connecting Components" => "manual/connections.md",
        ],
        "Tutorials" => [
            "Beginners Guide" => "tutorials/Beginners Guide.md",
            "The Engineering Data Model" => "tutorials/The Engineering Model.md",
            "Engineering Model: Helper Functions" => "tutorials/Engineering Model - Helper Functions.md",
            "Basics" => "tutorials/basic.md",
            "Explicit Neutral Models" => "tutorials/Explicit Neutral Models.md",
        ],
        "API Reference" => [
            "Base" => "reference/base.md",
            "Logging" => "reference/logging.md",
            "Data Models" => "reference/data_models.md",
            "Enums" => "reference/enums.md",
            "Formulations" => "reference/formulations.md",
            "Problems" => "reference/problems.md",
            "Variables" => "reference/variables.md",
            "Constraints" => "reference/constraints.md",
            "Objectives" => "reference/objectives.md",
            "Constants" => "reference/constants.md",
        ],
        "Developer Docs" => [
            "Extensions" => "developer/extensions.md",
            "Contributing" => "developer/contributing.md",
            "Style Guide" => "developer/style.md",
            "Roadmap" => "developer/roadmap.md",
        ],
    ]
)

# Insert HTML rendered from Pluto.jl into tutorial stubs as iframes
if !_FAST
    ss = Pluto.ServerSession()
    client = Pluto.ClientSession(Symbol("client", rand(UInt16)), nothing)
    ss.connected_clients[client.id] = client
    for file in readdir("examples", join=true)
        if endswith(file, ".jl")
            nb = Pluto.load_notebook_nobackup(file)
            client.connected_notebook = nb
            Pluto.update_run!(ss, nb, nb.cells)
            html = Pluto.generate_html(nb)

            fileout = "docs/build/tutorials/$(basename(file)).html"
            open(fileout, "w") do io
                write(io, html)
            end

            doc = open("docs/build/tutorials/$(replace(basename(file), ".jl" => ".html"))", "r") do io
                Gumbo.parsehtml(read(io, String))
            end

            # add style for full height iframe
            style = Gumbo.HTMLElement(:style)
            style.children = Gumbo.HTMLNode[Gumbo.HTMLText("iframe { height: 100vh; width: 100%; }")]
            push!(doc.root[1], style)

            # create iframe containing Pluto.jl rendered HTML
            iframe = Gumbo.HTMLElement(:iframe)
            iframe.attributes = Dict{AbstractString,AbstractString}(
                "src" => "$(basename(file)).html",
            )

            # edit existing html to replace :article with :iframe
            doc.root[2][1][2][2] = iframe

            # Overwrite HTML
            open("docs/build/tutorials/$(replace(basename(file), ".jl" => ".html"))", "w") do io
                Gumbo.prettyprint(io, doc)
            end
        end
    end
end


deploydocs(
    repo = "github.com/lanl-ansi/PowerModelsDistribution.jl.git",
    push_preview = false,
    devbranch = "main",
)
