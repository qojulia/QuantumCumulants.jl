
let
    NC_TIMES_RULES = [
        @rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~~x))
        @rule(~x::needs_sorting_nc => sort_args_nc(~x))

        @ordered_acrule(~a::SymbolicUtils.is_literal_number * ~b::SymbolicUtils.is_literal_number => ~a * ~b)

        @rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
        @ordered_acrule((~y)^(~n) * ~y => (~y)^(~n+1))
        @ordered_acrule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m))

        @ordered_acrule((~z::SymbolicUtils._isone  * ~x) => ~x)
        @ordered_acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)
        @rule(*(~x) => ~x)
    ]

    COMMUTATOR_RULES = [
        # Fock space rules
        @rule(*(~~a, ~x::SymbolicUtils.sym_isa(Destroy), ~y::SymbolicUtils.sym_isa(Create), ~~b) => apply_commutator(commute_bosonic, ~~a, ~~b, ~x, ~y))

        # NLevel rules
        @rule(*(~~a, ~x::SymbolicUtils.sym_isa(Transition), ~y::SymbolicUtils.sym_isa(Transition), ~~b) => apply_commutator(merge_transitions, ~~a, ~~b, ~x, ~y))
    ]

    EXPAND_TIMES_RULES = [
        @rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~x))
        @rule(~x::needs_sorting_nc => sort_args_nc(~x))
        @rule(*(~~a, +(~~b), ~~c) => +(map(b -> *((~~a)..., b, (~~c)...), ~~b)...))
    ]

    EXPAND_TIMES_COMMUTATOR_RULES = vcat(EXPAND_TIMES_RULES, COMMUTATOR_RULES)

    EXPAND_POW_RULES = [
        @rule(^(~x::SymbolicUtils.sym_isa(QNumber),~y::SymbolicUtils.isliteral(Integer)) => *((~x for i=1:~y)...))
    ]

    EXPAND_PLUS_RULES = [
        @rule(~x::SymbolicUtils.isnotflat(+) => SymbolicUtils.flatten_term(+, ~~x))
        @rule(~x::SymbolicUtils.needs_sorting(+) => SymbolicUtils.sort_args(+, ~~x))
        @ordered_acrule(~a::SymbolicUtils.is_literal_number + ~b::SymbolicUtils.is_literal_number => ~a + ~b)
        @acrule(*(~α::SymbolicUtils.is_literal_number, ~x) + ~x => *(~α + 1, ~x))
        @ordered_acrule((~z::SymbolicUtils._iszero + ~x) => ~x)
        @rule(+(~x) => ~x)
    ]

    CONJ_RULES = [
        @rule(conj(conj(~x)) => ~x)
        @rule(conj(sym_average(~x)) => sym_average(adjoint(~x)))
        @rule(conj(*(~~x)) => *(map(conj, ~~x)...))
    ]

    # Copied directly from SymbolicUtils
    PLUS_RULES = [
        @rule(~x::SymbolicUtils.isnotflat(+) => SymbolicUtils.flatten_term(+, ~~x))
        @rule(~x::SymbolicUtils.needs_sorting(+) => SymbolicUtils.sort_args(+, ~~x))
        @ordered_acrule(~a::SymbolicUtils.is_literal_number + ~b::SymbolicUtils.is_literal_number => ~a + ~b)

        @acrule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...))
        @acrule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...))

        @acrule(~x + *(~β, ~x) => *(1 + ~β, ~x))
        @acrule(*(~α::SymbolicUtils.is_literal_number, ~x) + ~x => *(~α + 1, ~x))
        @rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))

        @ordered_acrule((~z::SymbolicUtils._iszero + ~x) => ~x)
        @rule(+(~x) => ~x)
    ]

    POW_RULES = [
        @rule(^(*(~~x), ~y::SymbolicUtils.isliteral(Integer)) => *(map(a->SymbolicUtils.pow(a, ~y), ~~x)...))
        @rule((((~x)^(~p::SymbolicUtils.isliteral(Integer)))^(~q::SymbolicUtils.isliteral(Integer))) => (~x)^((~p)*(~q)))
        @rule(^(~x, ~z::SymbolicUtils._iszero) => 1)
        @rule(^(~x, ~z::SymbolicUtils._isone) => ~x)
    ]

    ASSORTED_RULES = [
        @rule(identity(~x) => ~x)
        @rule(-(~x) => -1*~x)
        @rule(-(~x, ~y) => ~x + -1(~y))
        @rule(~x / ~y => ~x * SymbolicUtils.pow(~y, -1))
        @rule(one(~x) => one(SymbolicUtils.symtype(~x)))
        @rule(zero(~x) => zero(SymbolicUtils.symtype(~x)))
        # @rule(SymbolicUtils.cond(~x::SymbolicUtils.is_literal_number, ~y, ~z) => ~x ? ~y : ~z)
    ]


    # Rewriter functions
    global operator_simplifier
    global commutator_simplifier
    global default_operator_simplifier
    global default_expand_simplifier
    global noncommutative_simplifier

    function default_operator_simplifier(; kwargs...)
        SymbolicUtils.IfElse(
            SymbolicUtils.sym_isa(QNumber), SymbolicUtils.Postwalk(operator_simplifier(); kwargs...),
            SymbolicUtils.Postwalk(number_simplifier(;kwargs...))
        )
    end

    function number_simplifier(;kwargs...)
        rule_tree = [SymbolicUtils.If(SymbolicUtils.is_operation(conj),
                                        SymbolicUtils.Chain(CONJ_RULES)
                                      ),
                    SymbolicUtils.default_simplifier(;kwargs...)
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Postwalk(rule_tree)
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
                    @rule(~x::SymbolicUtils.sym_isa(Transition) => rewrite_gs(~x))
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
