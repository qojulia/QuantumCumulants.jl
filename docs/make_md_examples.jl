using Literate

### Process examples
# Always rerun examples
const EXAMPLES_IN = joinpath(@__DIR__, "..", "examples")
const OUTPUT_MD_DIR = joinpath(@__DIR__, "src", "examples")

examples = filter!(file -> file[(end - 2):end] == ".jl", readdir(EXAMPLES_IN; join=true))
filter!(file -> !contains(file, "make_nb_examples"), examples)

if isempty(get(ENV, "CI", ""))
    # only needed when building docs locally; set automatically when built under CI
    # https://fredrikekre.github.io/Literate.jl/v2/outputformats/#Configuration
    extra_literate_config = Dict(
        "repo_root_path" => abspath(joinpath(@__DIR__, "..")),
        "repo_root_url" => "file://" * abspath(joinpath(@__DIR__, "..")),
    )
else
    extra_literate_config = Dict()
end

function preprocess(content)
    sub = SubstitutionString("""
                             """)
    content = replace(content, r"^# # [^\n]*"m => sub; count=1)

    # remove VSCode `##` block delimiter lines
    content = replace(content, r"^##$."ms => "")
    return content
end

for example in examples
    Literate.markdown(
        example,
        OUTPUT_MD_DIR;
        flavor=Literate.DocumenterFlavor(),
        config=extra_literate_config,
        # preprocess=preprocess,
    )
end
