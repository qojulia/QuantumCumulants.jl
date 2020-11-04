
let
    NC_TIMES_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~~x))
        SymbolicUtils.@rule(~x::needs_sorting_nc => sort_args_nc(~x))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)

        SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~y)^(~n) * ~y => (~y)^(~n+1)), 3)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m)), 3)

        SymbolicUtils.@rule(*(~~a, ~x!=~y, ~x==~y, ~~c) => false)
        SymbolicUtils.@rule(*(~~a, ~x==~y, ~x!=~y, ~~c) => false)
        SymbolicUtils.@rule(~x * !(~x) => false)

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._isone  * ~x) => ~x), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero *  ~x) => ~z), 2)
        SymbolicUtils.@rule(*(~x) => ~x)
    ]

    COMMUTATOR_RULES = [
        # Fock space rules
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(Destroy), ~y::SymbolicUtils.sym_isa(Create), ~~b) => apply_commutator(commute_bosonic, ~~a, ~~b, ~x, ~y))

        # NLevel rules
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(Transition), ~y::SymbolicUtils.sym_isa(Transition), ~~b) => apply_commutator(merge_transitions, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(~x::SymbolicUtils.sym_isa(Transition) => rewrite_gs(~x))

        # Indexed rules
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(IndexedDestroy), ~y::SymbolicUtils.sym_isa(IndexedCreate), ~~b) => apply_commutator(commute_bosonic_idx, ~~a, ~~b, ~x, ~y))

        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(IndexedTransition), ~y::SymbolicUtils.sym_isa(IndexedTransition), ~~b) => apply_commutator(merge_idx_transitions, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(IndexedTransition), ~y::SymbolicUtils.sym_isa(IndexedTransition), ~~b) => apply_commutator(merge_idx_transitions, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(~x::SymbolicUtils.sym_isa(IndexedTransition) => rewrite_gs(~x))

        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.is_operation(nip), ~y::SymbolicUtils.sym_isa(IndexedTransition), ~~b) => apply_commutator(merge_nip_idx_transition, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(IndexedTransition), ~y::SymbolicUtils.is_operation(nip), ~~b) => apply_commutator(merge_idx_transition_nip, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.is_operation(nip), ~y::SymbolicUtils.is_operation(nip), ~~b) => apply_commutator(merge_nips, ~~a, ~~b, ~x, ~y))

        # Sums
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(AbstractOperator), Sum(~y, ~~i), ~~b) => *((~~a)..., Sum(~x*~y, (~~i)...), (~~b)...))
        SymbolicUtils.@rule(*(~~a, Sum(~y, ~~i), ~x::SymbolicUtils.sym_isa(AbstractOperator), ~~b) => *((~~a)..., Sum(~y*~x, (~~i)...), (~~b)...))

    ]

    NIP_RULES = [
        # SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting(nip) => SymbolicUtils.sort_args(nip, ~x))
        SymbolicUtils.@rule(~x::needs_sorting_nip => sort_args_nip(~x))
        SymbolicUtils.@rule(nip(~~a, +(~~b), ~~c) => +(map(b -> nip((~~a)..., b, (~~c)...), ~~b)...))
        SymbolicUtils.@rule(nip(~~a, ~b::SymbolicUtils.isnumber, ~~d) => *(~b, nip((~~a)..., (~~d)...)))
        SymbolicUtils.@rule(nip(~~a, *(~b::SymbolicUtils.isnumber, ~c), ~~d) => *(~b, nip((~~a)..., ~c, (~~d)...)))
        SymbolicUtils.@rule(nip(~~a, *(~b, ~c::SymbolicUtils.isnumber), ~~d) => *(~c, nip((~~a)..., ~b, (~~d)...)))
        SymbolicUtils.@rule(nip(~x::SymbolicUtils.isnumber) => ~x)
        SymbolicUtils.@rule(nip(~x::SymbolicUtils.sym_isa(IndexedTransition)) => ~x)
    ]

    SUM_RULES = [
        SymbolicUtils.@rule(Sum(+(~~a), ~~i) => +(map(a -> Sum(a, (~~i)...), ~~a)...))
        SymbolicUtils.@rule(Sum(*(~~a, ~b::SymbolicUtils.isnumber, ~~c), ~~i) => *(~b, Sum(*((~~a)..., (~~c)...), (~~i)...)))
        SymbolicUtils.@rule(Sum(~x, ~~i::(!SymbolicUtils.issortedₑ)) => Sum(~x, sort_idx(~~i)...))

        SymbolicUtils.@rule(Sum(*(~~a, ~i==~j, ~~c), ~~k, ~i, ~~l) => Sum(swap_index(*((~~a)..., (~~c)...), ~i, ~j), (~~k)..., (~~l)...))
        SymbolicUtils.@rule(Sum(*(~~a, ~j==~i, ~~c), ~~k, ~i, ~~l) => Sum(swap_index(*((~~a)..., (~~c)...), ~i, ~j), (~~k)..., (~~l)...))
        SymbolicUtils.@rule(Sum(~x) => ~x)
        SymbolicUtils.@rule(Sum(~x::SymbolicUtils.isnumber, ~~idx) => _multiply_idxs_borders(~x, ~~idx))

        SymbolicUtils.@rule(~S::sum_has_const => sum_extract_const(~S))
    ]

    EXPAND_TIMES_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~x))
        SymbolicUtils.@rule(~x::needs_sorting_nc => sort_args_nc(~x))
        SymbolicUtils.@rule(*(~~a, +(~~b), ~~c) => +(map(b -> *((~~a)..., b, (~~c)...), ~~b)...))
    ]

    EXPAND_POW_RULES = [
        SymbolicUtils.@rule(^(~x::SymbolicUtils.sym_isa(AbstractOperator),~y::SymbolicUtils.isliteral(Integer)) => *((~x for i=1:~y)...))
    ]


    # Copied directly from SymbolicUtils
    PLUS_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(+) => SymbolicUtils.flatten_term(+, ~x))
        SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting(+) => SymbolicUtils.sort_args(+, ~x))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...)), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(~x + *(~β, ~x) => *(1 + ~β, ~x)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α::SymbolicUtils.isnumber, ~x) + ~x => *(~α + 1, ~x)), 2)
        SymbolicUtils.@rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))

        SymbolicUtils.@rule(+(~~a, ~x!=~y, ~x==~y, ~~c) => +((~~a)..., true, (~~c)...))
        SymbolicUtils.@rule(+(~~a, ~x==~y, ~x!=~y, ~~c) => +((~~a)..., true, (~~c)...))
        SymbolicUtils.@rule(~x + !(~x) => true)

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero + ~x) => ~x), 2)
        SymbolicUtils.@rule(+(~x) => ~x)
    ]

    POW_RULES = [
        SymbolicUtils.@rule(^(*(~~x), ~y::SymbolicUtils.isliteral(Integer)) => *(map(a->SymbolicUtils.pow(a, ~y), ~~x)...))
        SymbolicUtils.@rule((((~x)^(~p::SymbolicUtils.isliteral(Integer)))^(~q::SymbolicUtils.isliteral(Integer))) => (~x)^((~p)*(~q)))
        SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._iszero) => 1)
        SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._isone) => ~x)
        SymbolicUtils.@rule(^(~x == ~y, ~z) => (~x == ~y))
        SymbolicUtils.@rule(^(~x != ~y, ~z::SymbolicUtils._iszero) => true)
        SymbolicUtils.@rule(^(~x != ~y, ~z) => (~x != ~y))
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

    BOOLEAN_RULES = [
        SymbolicUtils.@rule((true | (~x)) => true)
        SymbolicUtils.@rule(((~x) | true) => true)
        SymbolicUtils.@rule((false | (~x)) => ~x)
        SymbolicUtils.@rule(((~x) | false) => ~x)
        SymbolicUtils.@rule((true & (~x)) => ~x)
        SymbolicUtils.@rule(((~x) & true) => ~x)
        SymbolicUtils.@rule((false & (~x)) => false)
        SymbolicUtils.@rule(((~x) & false) => false)

        SymbolicUtils.@rule(!(~x) & ~x => false)
        SymbolicUtils.@rule(~x & !(~x) => false)
        SymbolicUtils.@rule(!(~x) | ~x => true)
        SymbolicUtils.@rule(~x | !(~x) => true)
        SymbolicUtils.@rule(xor(~x, !(~x)) => true)
        SymbolicUtils.@rule(xor(~x, ~x) => false)

        SymbolicUtils.@rule(~x == ~x => true)
        SymbolicUtils.@rule(~x < ~x => false)
        SymbolicUtils.@rule(~x > ~x => false)
        SymbolicUtils.@rule(!(~x == ~y) => ~x != ~y)

        SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting((==)) => SymbolicUtils.sort_args((==), ~x))
        SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting((!=)) => SymbolicUtils.sort_args((!=), ~x))

        # simplify terms with no symbolic arguments
        # e.g. this simplifies term(isodd, 3, type=Bool)
        # or term(!, false)
        SymbolicUtils.@rule((~f)(~x::SymbolicUtils.isnumber) => (~f)(~x))
        # and this simplifies any binary comparison operator
        SymbolicUtils.@rule((~f)(~x::SymbolicUtils.isnumber, ~y::SymbolicUtils.isnumber) => (~f)(~x, ~y))
    ]



    # Rewriter functions
    global operator_simplifier
    global commutator_simplifier
    global default_operator_simplifier
    global default_expand_simplifier
    global noncommutative_simplifier

    function default_operator_simplifier(; kwargs...)
        SymbolicUtils.IfElse(
            SymbolicUtils.sym_isa(AbstractOperator), SymbolicUtils.Postwalk(operator_simplifier()),
            SymbolicUtils.default_simplifier(; kwargs...)
        )
    end

    function operator_simplifier()
        rw_comms = commutator_simplifier()
        rw_nc = noncommutative_simplifier()
        rw = SymbolicUtils.Chain([rw_comms,rw_nc])
        return SymbolicUtils.Postwalk(rw)
    end

    function commutator_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(EXPAND_TIMES_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(EXPAND_POW_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(nip), SymbolicUtils.Chain(NIP_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(Sum), SymbolicUtils.Chain(SUM_RULES)),
                    SymbolicUtils.Chain(COMMUTATOR_RULES)
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function noncommutative_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.is_operation(+), SymbolicUtils.Chain(PLUS_RULES)),
                     SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(NC_TIMES_RULES)),
                     SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(POW_RULES)),
                     SymbolicUtils.If(x->SymbolicUtils.symtype(x)<:Bool, SymbolicUtils.Chain(BOOLEAN_RULES))
                     ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function default_expand_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(EXPAND_TIMES_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(EXPAND_POW_RULES))
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end
end
