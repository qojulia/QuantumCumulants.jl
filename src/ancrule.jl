using SymbolicUtils: AbstractRule, Term, symtype, @timer, operation, arguments
import SymbolicUtils: Rule, getdepth

#### Associative Noncommutative Rules
struct ANCRule{F,R} <: AbstractRule
    sets::F
    rule::R
    arity::Int
end

Rule(acr::ANCRule)   = acr.rule
getdepth(r::ANCRule) = getdepth(r.rule)

Base.show(io::IO, acr::ANCRule) = print(io, "ANCRule(", acr.rule, ")")

function (acr::ANCRule)(term)
    r = Rule(acr)
    if !(term isa Term)
        r(term)
    else
        f =  operation(term)
        # Assume that the matcher was formed by closing over a term
        if f != operation(r.lhs) # Maybe offer a fallback if m.term errors.
            return nothing
        end

        T = symtype(term)
        args = arguments(term)

        itr = acr.sets(eachindex(args), acr.arity)

        for inds in itr
            result = r(Term{T}(f, @views args[inds]))
            if !isnothing(result)
                # Keep order
                args_l = args[1:inds[1]-1]
                args_r = args[inds[end]+1:end]
                return @timer "arule" Term{T}(f, [args_l..., result, args_r...])
            end
        end
    end
end
