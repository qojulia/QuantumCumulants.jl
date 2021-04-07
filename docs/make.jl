using Documenter
using QuantumCumulants

pages = [
        "index.md",
        # "implementation.md",
        "theory.md",
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
    sitename = "QuantumCumulants.jl",
    modules = [QuantumCumulants],
    pages = pages,
    checkdocs=:exports,
    format = Documenter.HTML(mathengine=MathJax())
    )

deploydocs(
    repo = "github.com/qojulia/QuantumCumulants.jl.git",
    )
