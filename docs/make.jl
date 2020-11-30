using Documenter
using Qumulants

pages = [
        "index.md",
        "tutorial.md",
        "correlation.md",
        "api.md",
        "Examples" => [
            "examples/single-atom-laser-spectrum.md"
            "examples/mollow.md"
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
