using Documenter
using Qumulants

pages = [
        "index.md",
        "tutorial.md",
        "api.md",
        "Examples" => [
            "examples/many-atom-laser.md"
            ]
    ]

makedocs(
    sitename = "Qumulants.jl",
    modules = [Qumulants],
    pages = pages,
    checkdocs=:exports,
    format = Documenter.HTML(mathengine=MathJax())
    )

deploydocs(
    repo = "github.com/david-pl/Qumulants.jl.git",
    )
