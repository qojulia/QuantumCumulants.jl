using QuantumCumulants
using BenchmarkTools

include("models.jl")

_meanfield(ops, H, J, rates, order) = meanfield(ops, H, J; rates = rates, order = order)
_complete(eqs, ::Nothing) = complete(eqs)
_complete(eqs, filt) = complete(eqs; filter_func = filt)
_find_missing(eqs, ::Nothing) = find_missing(eqs)
_find_missing(eqs, filt) = find_missing(eqs; filter_func = filt)
_corr(op1, op2, eqs, ::Nothing, ss) = CorrelationFunction(op1, op2, eqs; steady_state = ss)
_corr(op1, op2, eqs, filt, ss) =
    CorrelationFunction(op1, op2, eqs; steady_state = ss, filter_func = filt)

const SUITE = BenchmarkGroup()
for f in ("meanfield", "cumulant_expansion", "complete", "find_missing", "scale", "correlation")
    SUITE[f] = BenchmarkGroup()
end

sub!(g, k) = haskey(g, k) ? g[k] : (g[k] = BenchmarkGroup())

# evals=1: the symbolic passes are deterministic, so repeated evals add no signal.
const EVALS = 1
const SECONDS = 30

for m in models()
    name = m.name
    raw = m.pipeline ? meanfield(m.ops, m.H, m.J; rates = m.rates) : nothing
    avg = Dict(o => meanfield(m.ops, m.H, m.J; rates = m.rates, order = o) for o in m.orders)
    completed = Dict(o => _complete(avg[o], m.filt) for o in m.orders)

    if m.pipeline
        sub!(SUITE["meanfield"], name)
        sub!(SUITE["cumulant_expansion"], name)
        sub!(SUITE["find_missing"], name)
        SUITE["meanfield"][name]["order 2"] = @benchmarkable(
            _meanfield($(m.ops), $(m.H), $(m.J), $(m.rates), 2), evals = EVALS, seconds = SECONDS
        )
        SUITE["cumulant_expansion"][name]["order 2"] = @benchmarkable(
            cumulant_expansion($raw, 2), evals = EVALS, seconds = SECONDS
        )
        # find_missing on the open system; on a closed one it returns empty.
        SUITE["find_missing"][name]["order 2"] = @benchmarkable(
            _find_missing($(avg[2]), $(m.filt)), evals = EVALS, seconds = SECONDS
        )

        sub!(SUITE["complete"], name)
        for o in m.orders
            SUITE["complete"][name]["order $o"] = @benchmarkable(
                _complete($(avg[o]), $(m.filt)), evals = EVALS, seconds = SECONDS
            )
        end
    end

    if m.indexed
        sub!(SUITE["scale"], name)
        for o in m.orders
            SUITE["scale"][name]["order $o"] = @benchmarkable(
                scale($(completed[o])), evals = EVALS, seconds = SECONDS
            )
        end
    end

    if m.corr
        sub!(SUITE["correlation"], name)
        SUITE["correlation"][name]["order 2"] = @benchmarkable(
            _corr($(m.op1), $(m.op2), $(completed[2]), $(m.filt), $(m.steady)),
            evals = EVALS, seconds = SECONDS
        )
    end
end


results = BenchmarkTools.run(SUITE; verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
