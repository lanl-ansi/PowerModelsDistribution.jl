using Documenter, PowerModelsDistribution
import Pluto
import Literate

# TODO convert Pluto to Markdown
# ss = Pluto.ServerSession()
# client = Pluto.ClientSession(Symbol("client", rand(UInt16)), nothing)
# ss.connected_clients[client.id] = client

examples = []
for file in readdir("examples", join=true)
    if endswith(file, ".jl")
        # nb = Pluto.load_notebook_nobackup("examples/$file")
        # client.connected_notebook = nb
        # Pluto.update_run!(ss, nb, nb.cells)
        # html = Pluto.generate_html(nb)

        Literate.markdown(file, "docs/src"; documenter=true)
        push!(examples, replace(basename(file), ".jl" => ".md"))
        # fileout = "docs/src/$file"
        # open(fileout, "w") do io
        #     write(io, html)
        # end
    end
end

makedocs(
    modules = [PowerModelsDistribution],
    format = Documenter.HTML(analytics = "", mathengine = Documenter.MathJax()),
    sitename = "PowerModelsDistribution",
    authors = "David M Fobes, Carleton Coffrin, Frederik Geth, Sander Claeys, and contributors.",
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
