using QuantumCumulants
using ParallelTestRunner: ParallelTestRunner

# Start with autodiscovered tests
testsuite = ParallelTestRunner.find_tests(@__DIR__)

# `pending/` holds ports of master tests that depend on v1 features still
# being implemented (see TODO.md). The files are kept in the tree as
# source-of-truth so a future port is a one-line uncomment.
for name in collect(keys(testsuite))
    startswith(name, "pending/") && delete!(testsuite, name)
end

# Parse arguments
args = ParallelTestRunner.parse_args(ARGS)

if ParallelTestRunner.filter_tests!(testsuite, args)
    delete!(testsuite, "quality/JET")
end

ParallelTestRunner.runtests(QuantumCumulants, args; testsuite)
