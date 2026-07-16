using Documenter
using QuantumCumulants, SecondQuantizedAlgebra
using SymbolicUtils, ModelingToolkitBase

ENV["JULIA_DEBUG"] = "Documenter"
ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

include("make_md_examples.jl")

# The repo-root CHANGELOG.md is the single source of truth; expose it in the docs as the
# migration guide / changelog page. Write only when the content changed so a LiveServer
# build does not loop on its own output.
let src = normpath(@__FILE__, "../../CHANGELOG.md"),
        dst = normpath(@__FILE__, "../src/changelog.md"),
        new = read(src, String)

    if !isfile(dst) || read(dst, String) != new
        write(dst, new)
    end
end

# The repo-root README.md is the single source of truth for the landing page; expose it as
# the docs index. Same write-only-on-change guard as the changelog above.
let src = normpath(@__FILE__, "../../README.md"),
        dst = normpath(@__FILE__, "../src/index.md"),
        new = read(src, String)

    if !isfile(dst) || read(dst, String) != new
        write(dst, new)
    end
end

pages = [
    "index.md",
    "theory.md",
    "tutorial.md",
    "implementation.md",
    "symbolic_sums.md",
    "correlation.md",
    "noise.md",
    "ode_backends.md",
    "api.md",
    "Examples" => [
        "examples/single-atom-laser-spectrum.md"
        "examples/mollow.md"
        "examples/pulsed_mollow_spectrum.md"
        "examples/many-atom-laser.md"
        "examples/optomechanical-cooling.md"
        "examples/excitation-transport-chain.md"
        "examples/driven_dissipative_ising.md"
        "examples/ramsey_spectroscopy.md"
        "examples/superradiant_laser_indexed.md"
        "examples/weighted_phase_invariance.md"
        "examples/cavity_antiresonance_indexed.md"
        "examples/filter-cavity_indexed.md"
        "examples/unique_squeezing.md"
        "examples/waveguide.md"
        "examples/heterodyne_detection.md"
        # "examples/retrodiction_homodyne.md"

    ],
    "changelog.md",
    "devdocs.md",
]

using Pkg
status = sprint(io -> Pkg.status("SecondQuantizedAlgebra"; io = io))
version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
gh_moi = Documenter.Remotes.GitHub("qojulia", "SecondQuantizedAlgebra.jl")
remotes = Dict(pkgdir(SecondQuantizedAlgebra) => (gh_moi, version))

DocMeta.setdocmeta!(
    SecondQuantizedAlgebra,
    :DocTestSetup,
    :(using SecondQuantizedAlgebra);
    recursive = true,
)
DocMeta.setdocmeta!(
    QuantumCumulants,
    :DocTestSetup,
    :(using QuantumCumulants);
    recursive = true,
)

makedocs(
    sitename = "QuantumCumulants.jl",
    modules = [QuantumCumulants, SecondQuantizedAlgebra],
    pages = pages,
    remotes = remotes,
    checkdocs = :exports,
    doctest = false,
    # cross_references stay non-fatal: the only unresolved refs come from
    # SecondQuantizedAlgebra docstrings (e.g. `(i::Index)(k)`) and must be fixed upstream.
    warnonly = [:cross_references],
    format = Documenter.HTML(
        mathengine = MathJax(),
        footer = "[**Back to GitHub**](https://github.com/qojulia/QuantumCumulants.jl)",
        example_size_threshold = 800 * 2^10,
        size_threshold_warn = 400 * 2^10,
        size_threshold = 600 * 2^10,
    ),
)

deploydocs(repo = "github.com/qojulia/QuantumCumulants.jl.git", push_preview = true)
