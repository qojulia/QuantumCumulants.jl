using Literate

using QuantumCumulants, Plots
default(; fmt = :png)

### Process examples
# Always rerun examples
const EXAMPLES_IN = @__DIR__
const OUTPUT_NB_DIR = @__DIR__

examples = filter!(file -> file[(end-2):end] == ".jl", readdir(EXAMPLES_IN; join = true))
filter!(
    file -> !contains(file, "make_nb_examples") && !contains(file, "heterodyne_detection"),
    examples,
)

for example in examples
    Literate.notebook(example, OUTPUT_NB_DIR; documenter = false, execute = true)
end
