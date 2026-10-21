using Documenter
using DocumenterCodeBlocks
using DocumenterInterLinks
using DocumenterLandingPage
using QuantumCumulants
using SymbolicUtils, ModelingToolkitBase

ENV["JULIA_DEBUG"] = "Documenter"
ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

# DOCS_DRAFT=true skips the Literate examples and code execution for fast local previews.
DRAFT = get(ENV, "DOCS_DRAFT", "") == "true"

DRAFT || include("make_md_examples.jl")

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

include("pages.jl")

existing_pages(page::Pair) = page.second isa AbstractVector ?
    page.first => filter(!isnothing, map(existing_pages, page.second)) :
    (isfile(joinpath(@__DIR__, "src", page.second)) ? page : nothing)
existing_pages(page::AbstractString) = isfile(joinpath(@__DIR__, "src", page)) ? page : nothing
DRAFT && (pages = filter(!isnothing, map(existing_pages, pages)))

# Names re-exported from SecondQuantizedAlgebra are documented there. Every `@ref` that does
# not resolve locally falls back to SQA's inventory, so docstrings and pages keep plain refs.
links = InterLinks(
    "SecondQuantizedAlgebra" => "https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/",
)
fallbacks = ExternalFallbacks(; automatic = true)

DocMeta.setdocmeta!(
    QuantumCumulants,
    :DocTestSetup,
    :(using QuantumCumulants);
    recursive = true,
)

makedocs(
    sitename = "QuantumCumulants.jl",
    modules = [QuantumCumulants],
    pages = pages,
    plugins = [CodeBlocks(), LandingPage(), links, fallbacks],
    checkdocs = :exports,
    doctest = false,
    draft = DRAFT,
    format = Documenter.HTML(
        mathengine = KaTeX(),
        assets = ["assets/favicon.ico"],
        canonical = "https://qojulia.github.io/QuantumCumulants.jl/stable/",
        footer = "[**Back to GitHub**](https://github.com/qojulia/QuantumCumulants.jl)",
        example_size_threshold = 800 * 2^10,
        size_threshold_warn = 400 * 2^10,
        size_threshold = 600 * 2^10,
    ),
)

deploydocs(repo = "github.com/qojulia/QuantumCumulants.jl.git", push_preview = true)
