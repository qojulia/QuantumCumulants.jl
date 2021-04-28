using Documenter
using QuantumCumulants

pages = [
        "index.md",
        "theory.md",
        "tutorial.md",
        "correlation.md",
        "implementation.md",
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
    format = Documenter.HTML(
                            mathengine=MathJax(),
                            footer="[**Back to GitHub**](https://github.com/qojulia/QuantumCumulants.jl)"
                            )
    )

deploydocs(
    repo = "github.com/qojulia/QuantumCumulants.jl.git",
    )
