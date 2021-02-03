
let
    NC_TIMES_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~~x))
        SymbolicUtils.@rule(~x::needs_sorting_nc => sort_args_nc(~x))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.is_literal_number * ~b::SymbolicUtils.is_literal_number => ~a * ~b), 2)

        SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~y)^(~n) * ~y => (~y)^(~n+1)), 3)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m)), 3)

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._isone  * ~x) => ~x), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero *  ~x) => ~z), 2)
        SymbolicUtils.@rule(*(~x) => ~x)
    ]

    COMMUTATOR_RULES = [
        # Fock space rules
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(Destroy), ~y::SymbolicUtils.sym_isa(Create), ~~b) => apply_commutator(commute_bosonic, ~~a, ~~b, ~x, ~y))

        # NLevel rules
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(Transition), ~y::SymbolicUtils.sym_isa(Transition), ~~b) => apply_commutator(merge_transitions, ~~a, ~~b, ~x, ~y))
    ]

    EXPAND_TIMES_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~x))
        SymbolicUtils.@rule(~x::needs_sorting_nc => sort_args_nc(~x))
        SymbolicUtils.@rule(*(~~a, +(~~b), ~~c) => +(map(b -> *((~~a)..., b, (~~c)...), ~~b)...))
    ]

    EXPAND_TIMES_COMMUTATOR_RULES = vcat(EXPAND_TIMES_RULES, COMMUTATOR_RULES)

    EXPAND_POW_RULES = [
        SymbolicUtils.@rule(^(~x::SymbolicUtils.sym_isa(QNumber),~y::SymbolicUtils.isliteral(Integer)) => *((~x for i=1:~y)...))
    ]

    EXPAND_PLUS_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(+) => SymbolicUtils.flatten_term(+, ~~x))
        SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting(+) => SymbolicUtils.sort_args(+, ~~x))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.is_literal_number + ~b::SymbolicUtils.is_literal_number => ~a + ~b), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α::SymbolicUtils.is_literal_number, ~x) + ~x => *(~α + 1, ~x)), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero + ~x) => ~x), 2)
        SymbolicUtils.@rule(+(~x) => ~x)
    ]

    CONJ_RULES = [
        SymbolicUtils.@rule(conj(conj(~x)) => ~x)
        SymbolicUtils.@rule(conj(average(~x)) => average(adjoint(~x)))
        SymbolicUtils.@rule(conj(*(~~x)) => *(map(conj, ~~x)...))
    ]

    # Copied directly from SymbolicUtils
    PLUS_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(+) => SymbolicUtils.flatten_term(+, ~~x))
        SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting(+) => SymbolicUtils.sort_args(+, ~~x))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.is_literal_number + ~b::SymbolicUtils.is_literal_number => ~a + ~b), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...)), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(~x + *(~β, ~x) => *(1 + ~β, ~x)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α::SymbolicUtils.is_literal_number, ~x) + ~x => *(~α + 1, ~x)), 2)
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
        # SymbolicUtils.@rule(SymbolicUtils.cond(~x::SymbolicUtils.is_literal_number, ~y, ~z) => ~x ? ~y : ~z)
    ]


    # Rewriter functions
    global operator_simplifier
    global commutator_simplifier
    global default_operator_simplifier
    global default_expand_simplifier
    global noncommutative_simplifier
    global conj_rewriter

    function default_operator_simplifier(; kwargs...)
        SymbolicUtils.IfElse(
            SymbolicUtils.sym_isa(QNumber), SymbolicUtils.Postwalk(operator_simplifier(); kwargs...),
            SymbolicUtils.Postwalk(number_simplifier(;kwargs...))
        )
    end

    # For whatever reason, making this global fixes an "Unreachable reached" error
    # on Julia 1.5; should be removed at some point
    global _conj_rewriter
    _conj_rewriter() = SymbolicUtils.Chain(CONJ_RULES)

    function number_simplifier(;kwargs...)
        rule_tree = [SymbolicUtils.If(SymbolicUtils.is_operation(conj), _conj_rewriter()),
                    SymbolicUtils.default_simplifier(;kwargs...)
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Postwalk(rule_tree)
    end

    function conj_rewriter()
        rw = [SymbolicUtils.If(SymbolicUtils.is_operation(conj), SymbolicUtils.Chain(CONJ_RULES))
            ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rw))
    end


    function operator_simplifier()
        rw_comms = commutator_simplifier()
        rw_nc = noncommutative_simplifier()
        rw = SymbolicUtils.Chain([rw_comms,rw_nc])
        return SymbolicUtils.Postwalk(rw)
    end

    function commutator_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(EXPAND_TIMES_COMMUTATOR_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(EXPAND_POW_RULES)),
                    SymbolicUtils.@rule(~x::SymbolicUtils.sym_isa(Transition) => rewrite_gs(~x))
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function noncommutative_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.is_operation(+), SymbolicUtils.Chain(PLUS_RULES)),
                     SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(NC_TIMES_RULES)),
                     SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(POW_RULES))
                     ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function default_expand_simplifier(;kwargs...)
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(EXPAND_POW_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(EXPAND_TIMES_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(+), SymbolicUtils.Chain(EXPAND_PLUS_RULES)),
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Postwalk(rule_tree;kwargs...)
    end

    global serial_q_simplifier
    global serial_c_simplifier
    global threaded_q_simplifier
    global threaded_c_simplifier
    global serial_c_polynorm
    global serial_expand_simplifier
    global threaded_expand_simplifier

    serial_q_simplifier = default_operator_simplifier()
    serial_c_simplifier = SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Fixpoint(number_simplifier()))

    threaded_q_simplifier(cutoff) = default_operator_simplifier(threaded=true,
                                                              thread_cutoff=cutoff)
    threaded_c_simplifier(cutoff) = SymbolicUtils.Fixpoint(number_simplifier(threaded=true,
                                                              thread_cutoff=cutoff))

    serial_c_polynorm = SymbolicUtils.If(SymbolicUtils.istree,
                            SymbolicUtils.Fixpoint(SymbolicUtils.Chain((SymbolicUtils.polynormalize,
                                            SymbolicUtils.Fixpoint(number_simplifier())))))

    serial_expand_simplifier = SymbolicUtils.If(SymbolicUtils.istree,
                            SymbolicUtils.Fixpoint(default_expand_simplifier()))
    threaded_expand_simplifier(cutoff) = SymbolicUtils.Fixpoint(default_expand_simplifier(threaded=true,
                                                              thread_cutoff=cutoff))

    serial_expand_polynorm = SymbolicUtils.If(SymbolicUtils.istree,
                          SymbolicUtils.Fixpoint(SymbolicUtils.Chain((SymbolicUtils.polynormalize,
                                          SymbolicUtils.Fixpoint(default_expand_simplifier())))))
end
