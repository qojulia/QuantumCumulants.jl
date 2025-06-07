using Documenter
using QuantumCumulants, SecondQuantizedAlgebra

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
            "examples/waveguide.md"
            "examples/superradiant-laser.md"
            # "examples/heterodyne_detection.md"
            ]
    ]

using Pkg
status = sprint(io -> Pkg.status("SecondQuantizedAlgebra"; io = io))
version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
gh_moi = Documenter.Remotes.GitHub("qojulia", "SecondQuantizedAlgebra.jl")
remotes = Dict(pkgdir(SecondQuantizedAlgebra) => (gh_moi, version))

makedocs(
    sitename = "QuantumCumulants.jl",
    modules = [QuantumCumulants, SecondQuantizedAlgebra],
    pages = pages,
    remotes = remotes,
    checkdocs=:exports,
    format = Documenter.HTML(
                            mathengine=MathJax(),
                            footer="[**Back to GitHub**](https://github.com/qojulia/QuantumCumulants.jl)",
                            example_size_threshold = 800 * 2^10,
                            size_threshold_warn = 400 * 2^10,
                            size_threshold = 600 * 2^10,
                            )
    )

deploydocs(
    repo = "github.com/qojulia/QuantumCumulants.jl.git",
    push_preview = false,
    )
