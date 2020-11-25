names = [
    "test_fock.jl"
    "test_nlevel.jl"
    "test_embed.jl"
    "test_heisenberg.jl"
    "test_parameters.jl"
    "test_average.jl"
    "test_diffeq.jl"
    "test_v-level.jl"
    "test_mixed-order.jl"
    "test_correlation.jl"
    # "test_multisystems.jl"
    # "test_higher-order.jl"
]

detected_tests = filter(
    name->startswith(name, "test_") && endswith(name, ".jl"),
    readdir("."))

unused_tests = setdiff(detected_tests, names)
if length(unused_tests) != 0
    @warn string("The following tests are not used:\n", join(unused_tests, "\n"))
end

unavailable_tests = setdiff(names, detected_tests)
if length(unavailable_tests) != 0
    error("The following tests could not be found:\n", join(unavailable_tests, "\n"))
end

for name=names
    if startswith(name, "test_") && endswith(name, ".jl")
        include(name)
    end
end
