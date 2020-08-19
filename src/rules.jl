for f in [:isdestroy,:iscreate,:istransition,:issigmax,:issigmay,:issigmaz,
    :isspinop,:is_symspin,:isanyspin,:is_symmetrized_sigmax,:is_symmetrized_sigmay,:is_symmetrized_sigmaz,]
    @eval $(f)(a) = false
end
isspinN(N) = x->isspinN(N,_to_qumulants(x))
_isnotflat(f) = x -> _isnotflat(f,x)
_isnotflat(f,x) = SymbolicUtils.is_operation(f)(x) && SymbolicUtils.isnotflat(f)(x)

let
    NC_TIMES_RULES = [
        SymbolicUtils.@rule(~x::_isnotflat(*) => SymbolicUtils.flatten_term(*, ~~x))
        SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)

        SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~y)^(~n) * ~y => (~y)^(~n+1)), 3)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m)), 3)

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._isone  * ~x) => ~x), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero *  ~x) => ~z), 2)
        SymbolicUtils.@rule(*(~x) => ~x)
    ]

    COMMUTATOR_RULES = [
        # Fock space rules
        SymbolicUtils.@rule(*(~~a, ~x::isdestroy, ~y::iscreate, ~~b) => apply_commutator(commute_bosonic, ~~a, ~~b, ~x, ~y))

        # NLevel rules
        SymbolicUtils.@rule(*(~~a, ~x::istransition, ~y::istransition, ~~b) => apply_commutator(merge_transitions, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(~x::istransition => rewrite_gs(~x))

        # Spin space
        # Symmetrize
        SymbolicUtils.@rule(*(~~a, ~x::isspinop, ~y::is_symspin, ~~b) => *((~~a)..., symmetrize(~x), ~y, (~~b)...))
        SymbolicUtils.@rule(*(~~a, ~x::is_symspin, ~y::isspinop, ~~b) => *((~~a)..., ~x, symmetrize(~y), (~~b)...))
        # Replace with commutation relations
        SymbolicUtils.@rule(*(~~a, ~x::issigmaz, ~y::issigmax, ~~b) => apply_commutator(commute_spin, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(*(~~a, ~x::issigmay, ~y::issigmax, ~~b) => apply_commutator(commute_spin, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(*(~~a, ~x::issigmaz, ~y::issigmay, ~~b) => apply_commutator(commute_spin, ~~a, ~~b, ~x, ~y))
        # SymbolicUtils.@rule(*(~~a, ~x::isanyspin, ~y::isspinop, ~~b) => apply_commutator(commute_spin, ~~a, ~~b, ~x, ~y))

        # Special rules for spin-1/2 particles
        # SymbolicUtils.@rule(*(~~a, ~x::isspinN(1//2), ~y::isspinN(1//2), ~~b) => rewrite_spinhalf(~~a,~~b,~x,~y))
    ]

    EXPAND_RULES = [
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)

        # Expand powers
        SymbolicUtils.@rule(^(~x::SymbolicUtils.sym_isa(AbstractOperator),~y::SymbolicUtils.isliteral(Integer)) => *((~x for i=1:~y)...))
        SymbolicUtils.@rule(~x::_isnotflat(*) => SymbolicUtils.flatten_term(*, ~~x))
        SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))

        # Expand sums
        SymbolicUtils.@rule(*(~~a, +(~~b), ~~c) => +(map(b -> *((~~a)..., b, (~~c)...), ~~b)...))
        SymbolicUtils.@rule(~x::_isnotflat(+) => SymbolicUtils.flatten_term(+,~~x))
    ]


    # Copied directly from SymbolicUtils
    PLUS_RULES = [
        SymbolicUtils.@rule(~x::_isnotflat(+) => SymbolicUtils.flatten_term(+, ~~x))
        SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting(+) => SymbolicUtils.sort_args(+, ~~x))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...)), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(~x + *(~β, ~x) => *(1 + ~β, ~x)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α::SymbolicUtils.isnumber, ~x) + ~x => *(~α + 1, ~x)), 2)
        SymbolicUtils.@rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero + ~x) => ~x), 2)
        SymbolicUtils.@rule(+(~x) => ~x)
    ]

    POW_RULES = [
        SymbolicUtils.@rule(^(*(~~x), ~y::SymbolicUtils.isliteral(Integer)) => *(map(a->SymbolicUtils.pow(a, ~y), ~~x)...))
        SymbolicUtils.@rule((((~x)^(~p::SymbolicUtils.isliteral(Integer)))^(~q::SymbolicUtils.isliteral(Integer))) => (~x)^((~p)*(~q)))
        SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._iszero) => 1)
        SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._isone) => ~x)
    ]

    ASSORTED_RULES = [
        SymbolicUtils.@rule(identity(~x) => ~x)
        SymbolicUtils.@rule(-(~x) => -1*~x)
        SymbolicUtils.@rule(-(~x, ~y) => ~x + -1(~y))
        SymbolicUtils.@rule(~x / ~y => ~x * SymbolicUtils.pow(~y, -1))
        SymbolicUtils.@rule(one(~x) => one(SymbolicUtils.symtype(~x)))
        SymbolicUtils.@rule(zero(~x) => zero(SymbolicUtils.symtype(~x)))
        # SymbolicUtils.@rule(SymbolicUtils.cond(~x::SymbolicUtils.isnumber, ~y, ~z) => ~x ? ~y : ~z)
    ]

    UNSYMMETRIZE_RULES = [
        # After symmetrically ordered commutation relations have been applied, rewrite to original operators
        SymbolicUtils.@rule(~x::is_symspin => unsymmetrize(~x))
    ]

    ASSORTED_EXPAND_RULES = [EXPAND_RULES;ASSORTED_RULES]

    # Rewriter functions
    global default_commutator_simplifier
    global commutator_simplifier
    global default_expand_simplifier
    global default_noncommutative_simplifier
    global noncommutative_simplifier

    function default_commutator_simplifier(; kwargs...)
        SymbolicUtils.IfElse(
            SymbolicUtils.sym_isa(AbstractOperator), commutator_simplifier(),
            SymbolicUtils.default_simplifier(; kwargs...)
        )
    end

    function commutator_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_EXPAND_RULES)),
                    SymbolicUtils.Chain(COMMUTATOR_RULES)
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function default_noncommutative_simplifier(; kwargs...)
            SymbolicUtils.IfElse(
                SymbolicUtils.sym_isa(AbstractOperator), noncommutative_simplifier(),
                SymbolicUtils.default_simplifier(; kwargs...)
            )
    end

    function noncommutative_simplifier()
        rule_tree = [SymbolicUtils.Chain(ASSORTED_RULES),
                    SymbolicUtils.Chain(UNSYMMETRIZE_RULES),
                    SymbolicUtils.If(SymbolicUtils.is_operation(+), SymbolicUtils.Chain(PLUS_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(NC_TIMES_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(POW_RULES))
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function default_expand_simplifier()
        rule_tree = SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_EXPAND_RULES))
        SymbolicUtils.Postwalk(rule_tree)
    end
end
