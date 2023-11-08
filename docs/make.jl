using Documenter
using QuantumCumulants

ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

pages = [
        "index.md",
        "theory.md",
        "tutorial.md",
        "correlation.md",
        "symbolic_sums.md",
        "implementation.md",
        "api.md",
        "Examples" => [
            "examples/single-atom-laser-spectrum.md"
            "examples/mollow.md"
            "examples/many-atom-laser.md"
            "examples/optomechanical-cooling.md"
            "examples/excitation-transport-chain.md"
            "examples/ramsey_spectroscopy.md"
            "examples/superradiant_laser_indexed.md"
            "examples/cavity_antiresonance_indexed.md"
            "examples/filter-cavity_indexed.md"
            "examples/unique_squeezing.md"
            "examples/superradiant-laser.md"
            #"examples/heterodyne_detection.md"
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
    push_preview = false,
    )
