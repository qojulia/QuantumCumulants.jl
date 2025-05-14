
## Commutator methods
"""
    commutator(a,b)

Computes the commutator `a*b - b*a`.
"""
commutator(a,b) = _commutator(a,b)
_commutator(a, b) = a*b - b*a
commutator(a::QNumber,b::SNuN) = 0
commutator(a::SNuN,b::QNumber) = 0
commutator(::SNuN,::SNuN) = 0
function commutator(a::QSym,b::QSym)
    acts_on(a)==acts_on(b) || return 0
    isequal(a,b) && return 0
    return _commutator(a,b)
end
function commutator(a::QMul,b::QSym)
    aon = acts_on(b)
    idx = findfirst(x->isequal(acts_on(x),aon),a.args_nc)
    idx===nothing && return 0
    return _commutator(a,b)
end
function commutator(a::QSym,b::QMul)
    aon = acts_on(a)
    idx = findfirst(x->isequal(acts_on(x),aon),b.args_nc)
    idx===nothing && return 0
    return _commutator(a,b)
end
function commutator(a::QMul,b::QMul)
    # isequal(a.h, b.h) && return 0
    aon_a = map(acts_on, a.args_nc)
    aon_b = map(acts_on, b.args_nc)
    aon = intersect(aon_a,aon_b)
    isempty(aon) && return 0
    return _commutator(a,b)
end
function commutator(a::QAdd,b::QNumber)
    args = []
    for a_∈a.arguments
        c = commutator(a_,b)
        push_or_append_nz_args!(args, c)
    end
    isempty(args) && return 0
    return QAdd(args)
end
function commutator(a::QNumber,b::QAdd)
    args = []
    for b_∈b.arguments
        ## Commutator methods
        """
            commutator(a,b)

        Computes the commutator `a*b - b*a`.
        """
        commutator(a,b) = _commutator(a,b)
        _commutator(a, b) = a*b - b*a
        commutator(a::QNumber,b::SNuN) = 0
        commutator(a::SNuN,b::QNumber) = 0
        commutator(::SNuN,::SNuN) = 0
        function commutator(a::QSym,b::QSym)
            acts_on(a)==acts_on(b) || return 0
            isequal(a,b) && return 0
            return _commutator(a,b)
        end
        function commutator(a::QMul,b::QSym)
            aon = acts_on(b)
            idx = findfirst(x->isequal(acts_on(x),aon),a.args_nc)
            idx===nothing && return 0
            return _commutator(a,b)
        end
        function commutator(a::QSym,b::QMul)
            aon = acts_on(a)
            idx = findfirst(x->isequal(acts_on(x),aon),b.args_nc)
            idx===nothing && return 0
            return _commutator(a,b)
        end
        function commutator(a::QMul,b::QMul)
            # isequal(a.h, b.h) && return 0
            aon_a = map(acts_on, a.args_nc)
            aon_b = map(acts_on, b.args_nc)
            aon = intersect(aon_a,aon_b)
            isempty(aon) && return 0
            return _commutator(a,b)
        end
        function commutator(a::QAdd,b::QNumber)
            args = []
            for a_∈a.arguments
                c = commutator(a_,b)
                push_or_append_nz_args!(args, c)
            end
            isempty(args) && return 0
            return QAdd(args)
        end
        function commutator(a::QNumber,b::QAdd)
            args = []
            for b_∈b.arguments
                c = commutator(a,b_)
                push_or_append_nz_args!(args, c)
            end
            isempty(args) && return 0
            return QAdd(args)
        end
        function commutator(a::QAdd,b::QAdd)
            args = []
            for a_∈a.arguments, b_∈b.arguments
                c = commutator(a_,b_)
                push_or_append_nz_args!(args, c)
            end
            isempty(args) && return 0
            return QAdd(args)
        end

        function push_or_append_nz_args!(args,c)
            if !SymbolicUtils._iszero(c)
                push!(args, c)
            end
            return args
        end
        function push_or_append_nz_args!(args,c::QAdd)
            @inbounds for i=1:length(c.arguments)
                push_or_append_nz_args!(args, c.arguments[i])
            end
            return args
        end

        c = commutator(a,b_)
        push_or_append_nz_args!(args, c)
    end
    isempty(args) && return 0
    return QAdd(args)
end
function commutator(a::QAdd,b::QAdd)
    args = []
    for a_∈a.arguments, b_∈b.arguments
        c = commutator(a_,b_)
        push_or_append_nz_args!(args, c)
    end
    isempty(args) && return 0
    return QAdd(args)
end

function push_or_append_nz_args!(args,c)
    if !SymbolicUtils._iszero(c)
        push!(args, c)
    end
    return args
end
function push_or_append_nz_args!(args,c::QAdd)
    @inbounds for i=1:length(c.arguments)
        push_or_append_nz_args!(args, c.arguments[i])
    end
    return args
end
